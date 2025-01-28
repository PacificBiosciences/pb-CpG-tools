use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Arc;
use std::time::Duration;

use indicatif::{ProgressBar, ProgressStyle};
use log::info;
use thousands::Separable;

fn progress_status_report(log_info: &LogInfo) {
    let progress = log_info.progress.load(Ordering::SeqCst);
    let width = log_info.display_chars as usize;
    let percent = if log_info.total > 0 {
        (100 * progress) / log_info.total
    } else {
        0
    };
    info!(
        "{} {:>width$} of {:>width$} {} ({}%)",
        log_info.event_verb,
        progress.separate_with_commas(),
        log_info.total.separate_with_commas(),
        log_info.event_label,
        percent,
    );
}

fn progress_status_reporter(log_info: Arc<LogInfo>) {
    // Report interval in seconds:
    let progress_status_reporter_interval = 5 * 60;
    loop {
        std::thread::sleep(Duration::from_secs(progress_status_reporter_interval));
        if !log_info.run_reporter.load(Ordering::Relaxed) {
            return;
        }
        progress_status_report(&log_info);
    }
}

struct LogInfo {
    total: u64,
    event_verb: String,
    event_label: String,
    display_chars: u32,
    progress: AtomicU64,
    run_reporter: AtomicBool,
}

enum ProgressReporterType {
    Bar(ProgressBar),
    PeriodicLog(Arc<LogInfo>),
}

/// A generalized progress reporter covering both tty and non-tty output.
///
/// Progress reporter uses a single-line progress bar for tty contexts and a periodic interval status update for non-tty contexts.
///
pub struct ProgressReporter {
    pr_type: ProgressReporterType,
}

impl ProgressReporter {
    /// Create a new reporter instance
    ///
    /// The report format is: "{event_verb} {completed_count} of {event_count} {event label} ({percent}%)"
    /// For instance "Refined  1,696 of 42,056 breakpoint evidence clusters (4%)"
    ///
    /// The reporter will be formatted as either a single-line progress bar or as a periodic update.
    /// The progress bar will be selected if:
    ///   1. The terminal allows for this (the program is running on tty without output redirection, etc..)
    ///   2. the force_periodic flag is not true
    ///
    /// # Arguments
    /// * event_count - Total events to be completed
    /// * event_verb - Part of report phrase as described above
    /// * event_label - Part of report phrase as described above
    /// * force_periodic_updates - Force the progess output into the periodic format even if tty is available
    ///
    pub fn new(
        event_count: u64,
        event_verb: &str,
        event_label: &str,
        force_periodic_updates: bool,
    ) -> Self {
        let display_chars = {
            if event_count == 0 {
                1
            } else {
                let display_digits = event_count.ilog10() + 1;
                let display_commas = (display_digits - 1) / 3;
                display_digits + display_commas
            }
        };

        let template_string = format!(
            "[{{elapsed_precise}}] [{{bar:40}}] {event_verb} {{human_pos:>{display_chars}}} of {{human_len:{display_chars}}} {event_label} ({{percent}}%)");
        let progress_bar = ProgressBar::new(event_count).with_style(
            ProgressStyle::with_template(&template_string)
                .unwrap()
                .progress_chars("=> "),
        );

        use ProgressReporterType::*;
        // Determine whether this will be reported as a tty progress bar, or periodic log updates.
        //
        // For consistent tty determination we reuse the criteria from indicatif, ie. "bar.is_hidden()", so
        // that the periodic update is always provided when the progress bar is not.
        let use_periodic_updates = force_periodic_updates || progress_bar.is_hidden();
        if use_periodic_updates {
            let log_info = LogInfo {
                total: event_count,
                event_verb: event_verb.to_string(),
                event_label: event_label.to_string(),
                display_chars,
                progress: AtomicU64::new(0),
                run_reporter: AtomicBool::new(true),
            };
            let log_info = Arc::new(log_info);
            {
                let log_info = log_info.clone();
                std::thread::spawn(|| {
                    progress_status_reporter(log_info);
                });
            }

            Self {
                pr_type: PeriodicLog(log_info),
            }
        } else {
            progress_bar.tick();
            Self {
                pr_type: Bar(progress_bar),
            }
        }
    }

    /// Increment progress by `delta` completed events
    ///
    /// This is mutating the object, but the components are designed to internally
    /// manage this via locked access to the mutable state.
    ///
    pub fn inc(&self, delta: u64) {
        use ProgressReporterType::*;
        match &self.pr_type {
            Bar(progress_bar) => {
                progress_bar.inc(delta);
            }
            PeriodicLog(log_info) => {
                log_info.progress.fetch_add(delta, Ordering::SeqCst);
            }
        }
    }

    /// Clear progress bar from screen if in tty mode
    ///
    /// Shutdown updater thread if in non-tty mode
    ///
    pub fn clear(&self) {
        use ProgressReporterType::*;
        match &self.pr_type {
            Bar(progress_bar) => {
                progress_bar.finish_and_clear();
            }
            PeriodicLog(log_info) => {
                log_info.run_reporter.store(false, Ordering::Relaxed);
            }
        }
    }
}

impl Drop for ProgressReporter {
    fn drop(&mut self) {
        self.clear();
    }
}
