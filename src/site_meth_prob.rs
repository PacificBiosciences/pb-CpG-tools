use rust_vc_utils::norm_ln_distro;

/// Report the expected methylation frequency at a site from that site's base methylation
/// probabilities, given that only a specified set of methylation frequencies are possible at a site
///
/// Assume that the prior over possible site methylation frequencies is uniform.
///
/// possible_site_meth_frequencies: An array of values between 0.0 and 1.0. The model assumes that
/// the site methylation frequency must be one of these values
///
/// base_meth_probs: Observed methylation probabilities for each base at a site
///
/// other_count: Number of additional background observations not included in base_meth_probs
///
/// other_meth_prob: Methylation probability for additional background observations
///
/// Returns Expectation[possible_site_meth_frequencies]
///
pub fn get_site_meth_prob(
    possible_site_meth_frequencies: &[f32],
    base_meth_probs: &[f32],
    other_count: u32,
    other_meth_prob: f32,
) -> f32 {
    let mut state_ln_lhood = Vec::new();
    for freq in possible_site_meth_frequencies.iter() {
        let mut ln_lhood = 0.0f64;
        for base_meth_prob in base_meth_probs.iter() {
            ln_lhood +=
                ((base_meth_prob * freq) + ((1.0 - base_meth_prob) * (1.0 - freq))).ln() as f64;
        }
        ln_lhood += other_count as f64
            * ((other_meth_prob * freq) + ((1.0 - other_meth_prob) * (1.0 - freq))).ln() as f64;
        state_ln_lhood.push(ln_lhood);
    }

    norm_ln_distro(&mut state_ln_lhood);
    let state_post_prob = state_ln_lhood;

    // Get expected meth freq from the posterior distro over all possible frequencies:
    let mut expect_freq = 0.0f32;
    for (&post_prob, &freq) in state_post_prob
        .iter()
        .zip(possible_site_meth_frequencies.iter())
    {
        expect_freq += post_prob as f32 * freq;
    }
    expect_freq
}

const HAPLOID_FREQS: &[f32] = &[0.0, 1.0];

#[allow(dead_code)]
pub fn get_haploid_site_meth_prob(base_meth_probs: &[f32]) -> f32 {
    get_site_meth_prob(HAPLOID_FREQS, base_meth_probs, 0, 0.1)
}

#[allow(dead_code)]
pub fn get_haploid_site_meth_prob_with_other_count(
    base_meth_probs: &[f32],
    other_count: u32,
    other_meth_prob: f32,
) -> f32 {
    get_site_meth_prob(HAPLOID_FREQS, base_meth_probs, other_count, other_meth_prob)
}

const DIPLOID_FREQS: &[f32] = &[0.0, 0.5, 1.0];

#[allow(dead_code)]
pub fn get_diploid_site_meth_prob(base_meth_probs: &[f32]) -> f32 {
    get_site_meth_prob(DIPLOID_FREQS, base_meth_probs, 0, 0.1)
}

#[allow(dead_code)]
pub fn get_diploid_site_meth_prob_with_other_count(
    base_meth_probs: &[f32],
    other_count: u32,
    other_meth_prob: f32,
) -> f32 {
    get_site_meth_prob(DIPLOID_FREQS, base_meth_probs, other_count, other_meth_prob)
}

#[allow(dead_code)]
pub fn get_uniform_site_meth_prob_with_other_count(
    base_meth_probs: &[f32],
    other_count: u32,
    other_meth_prob: f32,
) -> f32 {
    let sum = base_meth_probs.iter().sum::<f32>() + other_meth_prob * other_count as f32;
    let total = base_meth_probs.len() + other_count as usize;
    sum / total as f32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_site_meth_prob() {
        let val = get_haploid_site_meth_prob(&[0.9, 0.9]);
        approx::assert_abs_diff_eq!(val, 0.9878_f32, epsilon = 0.001);

        let val = get_haploid_site_meth_prob(&[0.1, 0.1]);
        approx::assert_abs_diff_eq!(val, 0.0122_f32, epsilon = 0.001);

        let val1 = get_haploid_site_meth_prob_with_other_count(&[0.9, 0.9], 2, 0.1);
        let val2 = get_haploid_site_meth_prob(&[0.9, 0.9, 0.1, 0.1]);
        approx::assert_abs_diff_eq!(val1, val2, epsilon = 0.001);

        let val = get_diploid_site_meth_prob(&[0.9, 0.9]);
        approx::assert_abs_diff_eq!(val, 0.8738_f32, epsilon = 0.001);

        let val = get_diploid_site_meth_prob(&[0.1, 0.1]);
        approx::assert_abs_diff_eq!(val, 0.1262_f32, epsilon = 0.001);

        let val1 = get_diploid_site_meth_prob_with_other_count(&[0.9, 0.9], 2, 0.1);
        let val2 = get_diploid_site_meth_prob(&[0.9, 0.9, 0.1, 0.1]);
        approx::assert_abs_diff_eq!(val1, val2, epsilon = 0.001);
    }
}
