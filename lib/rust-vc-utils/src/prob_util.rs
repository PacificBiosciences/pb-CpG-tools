/// returns equiv of log(exp(x1)+exp(x2)), requires that x1 >= x2
///
fn log_sum_sorted(x1: f64, x2: f64) -> f64 {
    x1 + (1.0 + (x2 - x1).exp()).ln()
}

/// returns equiv of log(exp(x1)+exp(x2))
///
pub fn log_sum(x1: f64, x2: f64) -> f64 {
    if x1 < x2 {
        log_sum_sorted(x2, x1)
    } else {
        log_sum_sorted(x1, x2)
    }
}

/// Normalize log-transformed probability distro in place
///
/// Returns index of the most probable component
///
pub fn norm_ln_distro(c: &mut [f64]) -> Option<usize> {
    // scale and exp pprob values:
    let mut is_first = true;
    let mut max_val = 0.0;
    let mut max_element = 0;
    for (i, &v) in c.iter().enumerate() {
        if is_first || v > max_val {
            is_first = false;
            max_val = v;
            max_element = i;
        }
    }

    if is_first {
        return None;
    }

    let mut sum = 0_f64;
    for v in c.iter_mut() {
        *v = (*v - max_val).exp();
        sum += *v;
    }

    // normalize:
    let sum = 1.0 / sum;
    for v in c.iter_mut() {
        *v *= sum;
    }

    Some(max_element)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_log_sum() {
        approx::assert_ulps_eq!(log_sum(0.2_f64.ln(), 0.3_f64.ln()).exp(), 0.5, max_ulps = 4);
        approx::assert_ulps_eq!(log_sum(0.3_f64.ln(), 0.2_f64.ln()).exp(), 0.5, max_ulps = 4);
    }

    #[test]
    fn test_norm_ln_distro() {
        let mut distro = Vec::new();
        let res = norm_ln_distro(&mut distro);
        assert!(res.is_none());

        let mut distro = vec![0.1_f64.ln(), 0.1_f64.ln(), 0.1_f64.ln(), 0.2_f64.ln()];
        let res = norm_ln_distro(&mut distro);
        assert_eq!(res, Some(3));
        approx::assert_ulps_eq!(distro[3] as f32, 0.4_f32, max_ulps = 4);

        let mut distro = vec![
            0.0001_f64.ln(),
            0.0001_f64.ln(),
            0.0001_f64.ln(),
            0.0002_f64.ln(),
        ];
        let res = norm_ln_distro(&mut distro);
        assert_eq!(res, Some(3));
        approx::assert_ulps_eq!(distro[3] as f32, 0.4_f32, max_ulps = 4);
    }
}
