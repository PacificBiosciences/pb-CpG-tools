use std::collections::BTreeMap;
use std::fmt;

use clap::ValueEnum;
use error_stack::{Context, Report, ResultExt};
use rust_htslib::bam;
use rust_vc_utils::aux::get_optional_int_aux_tag;
use rust_vc_utils::basemod::{decode_cpg_meth_info, CpgMethInfo};
use rust_vc_utils::get_complete_read_to_ref_pos_map;
use rust_vc_utils::RingBuffer;
use serde_derive::{Deserialize, Serialize};
use unwrap::unwrap;

use crate::chrom_sites::ChromSiteData;
use crate::tflite_model::TFLiteModelData;

/// Number of contiguous CpGs used as input to the pileup model
const MODEL_CPG_WINDOW_COUNT: usize = 11;

/// 0-index of the CpG that is actually called within the CpG window
const CALLED_CPG_INDEX: usize = 6;

/// The number of bins used for the meth prob distribution
const MODEL_METH_PROB_BIN_COUNT: usize = 20;

/// The number of pileup model features per CpG
const MODEL_FEATURES_PER_CPG: usize = 21;

/// Get a histogram of probabilities, L2-normalized
///
/// * `zero_prob_count` - The number of additional observation counts with probability zero, not
///   already included in probs.
///
fn get_norm_prob_histogram(
    probs: &[f32],
    bin_count: usize,
    zero_prob_count: u32,
    prob_offset: f32,
) -> Vec<f32> {
    assert!(bin_count > 0);

    fn prob_to_bin(prob: f32, prob_offset: f32, bin_count: usize) -> usize {
        (((prob + prob_offset) * bin_count as f32) as i32).clamp(0, bin_count as i32 - 1) as usize
    }
    let mut hist = vec![0.0f64; bin_count];

    for &prob in probs.iter() {
        hist[prob_to_bin(prob, prob_offset, bin_count)] += 1.0;
    }

    hist[0] += zero_prob_count as f64;

    // L2 norm
    let l2n = hist.iter().map(|x| x * x).sum::<f64>().sqrt();
    hist.iter().map(|x| (x / l2n) as f32).collect::<Vec<_>>()
}

/// Meth summary features used for a given site in the reference
#[derive(Clone, Deserialize, Serialize)]
pub struct SiteMethTrainingData {
    /// Feature volume for the target site.
    /// Outer vector contains one entry for each position in the sequence window.
    /// Inner vector contains the feature set for each position.
    ///
    pub features: Vec<Vec<f32>>,

    /// Not currently used as a training feature but added here in case useful
    pub cov: u32,

    /// This is the expected methylation fraction as specified in the sites input file. It will be
    /// passed on to the training routine as an output label.
    pub meth_frac: f32,
}

/// Object to accumulate all processed site features for a single chromosome
///
/// Information in this object will be translated into final json feature output
///
#[derive(Clone, Default, Deserialize, Serialize)]
pub struct ChromSiteMethTrainingData {
    /// key is chrom position (0-indexed)
    /// value is feature output for this the pileup model at this position
    ///
    pub data: BTreeMap<i64, SiteMethTrainingData>,
}

impl ChromSiteMethTrainingData {
    pub fn new() -> Self {
        Self {
            data: BTreeMap::new(),
        }
    }

    /// Move site feature data from other into this object, but require that it is non-conflicting
    ///
    pub fn merge(&mut self, other: Self) {
        // Use the explicit loop here instead of extend so that key collisions can be checked:
        for (k, v) in other.data {
            let val = self.data.insert(k, v);
            assert!(
                val.is_none(),
                "Key collision between genome segments at ref position {}",
                k + 1
            );
        }
    }
}

/// Site methylation features for the whole genome
///
/// This will include any sites specified in the input sites file
///
#[derive(Deserialize, Serialize)]
pub struct GenomeSiteMethTrainingData {
    pub chroms: BTreeMap<String, ChromSiteMethTrainingData>,
}

impl GenomeSiteMethTrainingData {
    pub fn site_count(&self) -> usize {
        self.chroms.values().map(|x| x.data.len()).sum()
    }
}

/// Simulates python round(x, 1)
pub fn python_round1(x: f64) -> f64 {
    (x * 10.0).round() / 10.0
}

/// Meth summary information for each site in the reference
#[derive(Clone)]
pub enum SiteMethInfo {
    /// Info produced in 'model' pileup mode
    Model {
        /// Output of the deep learning model
        pileup_model_prob: f32,
        cov: u32,
    },
    /// Info produced in 'count' pileup mode
    Count {
        meth_cov: u32,
        unmeth_cov: u32,
        mean_meth_prob: f32,
        mean_unmeth_prob: f32,
    },
}

impl SiteMethInfo {
    /// Provide the summary methylation score for this site that matches previous pileup python
    /// script
    pub fn pileup_output_score(&self) -> f64 {
        match self {
            SiteMethInfo::Model {
                pileup_model_prob, ..
            } => python_round1(100.0 * *pileup_model_prob as f64),
            SiteMethInfo::Count {
                meth_cov,
                unmeth_cov,
                ..
            } => {
                if *meth_cov == 0 {
                    0.0
                } else {
                    let cov = *meth_cov + *unmeth_cov;
                    python_round1(100.0 * *meth_cov as f64 / cov as f64)
                }
            }
        }
    }

    /// Summary coverage that matches python script
    pub fn cov(&self) -> u32 {
        match self {
            SiteMethInfo::Model { cov, .. } => *cov,
            SiteMethInfo::Count {
                meth_cov,
                unmeth_cov,
                ..
            } => *meth_cov + *unmeth_cov,
        }
    }
}

pub type SiteMethProbTrack = BTreeMap<i64, SiteMethInfo>;

#[derive(Copy, Clone)]
pub enum MethInfoTrackType {
    Combined,
    Hap1,
    Hap2,
}

pub fn get_meth_info_track_label(tt: MethInfoTrackType) -> &'static str {
    use MethInfoTrackType::*;
    match tt {
        Combined => "combined",
        Hap1 => "hap1",
        Hap2 => "hap2",
    }
}

pub fn get_python_script_bed_labels(tt: MethInfoTrackType) -> &'static str {
    use MethInfoTrackType::*;
    match tt {
        Combined => "Total",
        Hap1 => "hap1",
        Hap2 => "hap2",
    }
}

pub fn skip_if_empty(tt: MethInfoTrackType) -> bool {
    use MethInfoTrackType::*;
    match tt {
        Combined => false,
        Hap1 => true,
        Hap2 => true,
    }
}

/// Object to accumulate all processed site data for a single chromosome
///
/// Information in this object will be translated into final output bed/bigwig files
///
#[derive(Clone, Default)]
pub struct ChromSiteMethInfo {
    /// Site probabilities from all reads, all hap1 reads, and all hap2 reads
    pub track_site_probs: Vec<SiteMethProbTrack>,
}

impl ChromSiteMethInfo {
    pub fn new() -> Self {
        Self {
            track_site_probs: vec![BTreeMap::new(); HAPLOTYPE_TRACK_COUNT],
        }
    }

    pub fn get_track(&self, tt: MethInfoTrackType) -> &SiteMethProbTrack {
        use MethInfoTrackType::*;
        match tt {
            Combined => &self.track_site_probs[0],
            Hap1 => &self.track_site_probs[1],
            Hap2 => &self.track_site_probs[2],
        }
    }

    /// Move site prob data from other into this object, but require that it is non-conflicting
    ///
    pub fn merge(&mut self, other: Self) {
        fn merge_track(self_site_prob: &mut SiteMethProbTrack, other_site_prob: SiteMethProbTrack) {
            // Use the explicit loop here instead of extend so that key collisions can be checked:
            for (k, v) in other_site_prob {
                let val = self_site_prob.insert(k, v);
                assert!(
                    val.is_none(),
                    "Key collision between genome segments at ref position {k}"
                );
            }
        }

        for (track_id, other_track) in other.track_site_probs.into_iter().enumerate() {
            merge_track(&mut self.track_site_probs[track_id], other_track);
        }
    }
}

const BASE_COUNT: usize = 5;

#[allow(dead_code)]
pub enum Bases {
    A,
    C,
    G,
    T,
    N,
}

/// Summarize all base observations at a site for a given haplotype category.
///
/// Note that "haplotype category" includes the unphased/unknown haplotype code '0'
///
#[derive(Clone)]
pub struct HaplotypeMethProbPileupColumn {
    /// List of all observed base methylation probabilities
    pub meth_probs: Vec<f32>,

    /// Count Bases to help determine the majority base
    pub base_counts: [u32; BASE_COUNT],
}

impl HaplotypeMethProbPileupColumn {
    pub fn new() -> Self {
        Self {
            meth_probs: Vec::new(),
            base_counts: [0; BASE_COUNT],
        }
    }

    pub fn clear(&mut self) {
        self.meth_probs.clear();
        self.base_counts = [0; BASE_COUNT];
    }

    pub fn is_empty(&self) -> bool {
        self.meth_probs.is_empty() && (self.cov() == 0)
    }

    pub fn add_base(&mut self, base: u8) {
        let index = match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 4,
        };
        self.base_counts[index] += 1;
    }

    /// Count of all non-filtered reads at this site corresponding to the given haplotype
    pub fn cov(&self) -> u32 {
        self.base_counts.iter().sum::<u32>()
    }

    /// Produce all summary stats used in the 'count' pileup mode
    ///
    /// Note that bases are given a meth prob of 0 if the value is missing, so these are added to
    /// the unmeth_cov count first.
    ///
    pub fn get_count_mode_summary_stats(&self) -> SiteMethInfo {
        let mut meth_cov = 0;
        let mut unmeth_cov = self.cov() - self.meth_probs.len() as u32;
        let mut mean_meth_prob = 0.0;
        let mut mean_unmeth_prob = 0.0;
        for &meth_prob in self.meth_probs.iter() {
            if meth_prob > 0.5 {
                meth_cov += 1;
                mean_meth_prob += meth_prob;
            } else {
                unmeth_cov += 1;
                mean_unmeth_prob += meth_prob;
            }
        }

        if meth_cov > 0 {
            mean_meth_prob /= meth_cov as f32;
        }
        if unmeth_cov > 0 {
            mean_unmeth_prob /= unmeth_cov as f32;
        }

        SiteMethInfo::Count {
            meth_cov,
            unmeth_cov,
            mean_meth_prob,
            mean_unmeth_prob,
        }
    }

    pub fn is_most_probable_base(&self, base: Bases) -> bool {
        let query_index = base as usize;
        for base_index in 0..BASE_COUNT {
            if base_index == query_index {
                continue;
            }
            if self.base_counts[base_index] >= self.base_counts[query_index] {
                return false;
            }
        }
        true
    }

    pub fn merge(&mut self, other: &Self) {
        self.meth_probs.extend(&other.meth_probs);
        for base_index in 0..BASE_COUNT {
            self.base_counts[base_index] += other.base_counts[base_index];
        }
    }
}

const HAPLOTYPE_TRACK_COUNT: usize = 3;

/// Summarize all base observations at a site
#[derive(Clone)]
pub struct MethProbPileupColumn {
    haps: Vec<HaplotypeMethProbPileupColumn>,
}

impl MethProbPileupColumn {
    pub fn new() -> Self {
        Self {
            haps: vec![HaplotypeMethProbPileupColumn::new(); HAPLOTYPE_TRACK_COUNT],
        }
    }

    pub fn clear(&mut self) {
        self.haps.iter_mut().for_each(|x| x.clear());
    }

    pub fn is_empty(&self) -> bool {
        self.haps.iter().all(|x| x.is_empty())
    }

    /*
    pub fn cov(&self) -> u32 {
        self.haps.iter().map(|x| x.cov()).sum()
    }
     */

    /// Get a pileup for the 'total' merge of all (disjoint) haplotypes in haps
    ///
    pub fn get_hap_total(&self) -> HaplotypeMethProbPileupColumn {
        let mut total = HaplotypeMethProbPileupColumn::new();
        for x in self.haps.iter() {
            total.merge(x);
        }
        total
    }
}

/// Model features derived from each CpG site on a given haplotype track
#[derive(Clone)]
struct TrackCpGFeatures {
    ref_pos: i64,
    cov: u32,
    features: Vec<f32>,
}

impl TrackCpGFeatures {
    fn new() -> Self {
        Self {
            ref_pos: 0,
            cov: 0,
            features: Vec::new(),
        }
    }

    fn blank() -> Self {
        let mut x = Self::new();
        x.features = vec![0f32; MODEL_FEATURES_PER_CPG];
        x
    }
}

/// A ring buffer holding data from the last N CpG sites to be used as input to the CpG pileup model
///
type FeatureRing = RingBuffer<TrackCpGFeatures>;

struct CpGFeatureInfo {
    track_rings: Vec<FeatureRing>,
}

impl CpGFeatureInfo {
    pub fn new(is_chrom_start: bool) -> Self {
        fn init_ring(is_chrom_start: bool) -> FeatureRing {
            let mut x = FeatureRing::new(MODEL_CPG_WINDOW_COUNT);

            // At the start of a chromosome, insert several columns of zero features to get the
            // buffer started
            if is_chrom_start {
                for _ in 0..CALLED_CPG_INDEX {
                    x.push_item(TrackCpGFeatures::blank());
                }
            }
            x
        }

        Self {
            track_rings: vec![init_ring(is_chrom_start); HAPLOTYPE_TRACK_COUNT],
        }
    }
}

struct PostPileupData<'a> {
    last_cpg_ref_pos: Option<i64>,

    /// CpG feature information used as input to the CNN-model if pileup_mode is Model
    ///
    cpg_feature_info: CpGFeatureInfo,

    /// Processed site methylation information. This data is exported out of the object when read
    /// processing is complete.
    site_meth_info: ChromSiteMethInfo,

    /// Fully processed site CpG features information. This data is exported out of the object when
    /// read processing is complete.
    chrom_site_meth_training_data: ChromSiteMethTrainingData,

    /// End of site processing range - use this position to ensure that this bin does not duplicate
    /// content printed by the next bin
    end_pos: i64,

    /// Methylation probabilities will only be generated and output for a site with at least this
    /// much coverage
    min_coverage: u32,

    /// CpG sites must have at least this much coverage to be used as feature inputs to the pileup
    /// model
    min_pileup_model_feature_coverage: u32,

    /// Select pileup processing mode
    pileup_mode: MethReadsProcessorPileupMode,

    modsites_mode: ModSitesMode,

    chrom_training_sites: Option<&'a ChromSiteData>,

    chrom_ref: Option<&'a Vec<u8>>,
}

impl<'a> PostPileupData<'a> {
    pub fn new(
        is_chrom_start: bool,
        end_pos: i64,
        min_coverage: u32,
        pileup_mode: MethReadsProcessorPileupMode,
        modsites_mode: ModSitesMode,
        chrom_training_sites: Option<&'a ChromSiteData>,
        chrom_ref: Option<&'a Vec<u8>>,
    ) -> Self {
        Self {
            last_cpg_ref_pos: None,
            cpg_feature_info: CpGFeatureInfo::new(is_chrom_start),
            site_meth_info: ChromSiteMethInfo::new(),
            chrom_site_meth_training_data: ChromSiteMethTrainingData::new(),
            end_pos,
            min_coverage,
            min_pileup_model_feature_coverage: 4,
            pileup_mode,
            modsites_mode,
            chrom_training_sites,
            chrom_ref,
        }
    }

    /// Return true if any of the cpg feature rings still contain content < end_pos
    pub fn unprocessed_cpgs(&self) -> bool {
        fn unprocessed_cpgs(ring: &FeatureRing, end_pos: i64) -> bool {
            !ring.is_empty() && ring.get_item(0).ref_pos < end_pos
        }

        for track_ring in self.cpg_feature_info.track_rings.iter() {
            if unprocessed_cpgs(track_ring, self.end_pos) {
                return true;
            }
        }
        false
    }

    /// Return true if the site represented by meth_prob_pileup_column is the beginning of a
    /// qualifying CpG sequence.
    ///
    /// To be a qualifying CpG sequence. 'C' and 'G', must both be the majority bases observed at
    /// the site and next site respectively in any one of the 3 tracks (total, hap1, hap2)
    ///
    /// Note that total_col and next_total_col could be found directly from the pileup_column inputs
    /// and are just taken here to (theoretically) save time.
    ///
    fn keep_cpg_site(
        total_col: &HaplotypeMethProbPileupColumn,
        next_total_col: &HaplotypeMethProbPileupColumn,
        meth_prob_pileup_column: &MethProbPileupColumn,
        next_meth_prob_pileup_column: &MethProbPileupColumn,
    ) -> bool {
        fn test_track(
            hap_col: &HaplotypeMethProbPileupColumn,
            next_hap_col: &HaplotypeMethProbPileupColumn,
        ) -> bool {
            hap_col.is_most_probable_base(Bases::C) && next_hap_col.is_most_probable_base(Bases::G)
        }

        test_track(total_col, next_total_col)
            || test_track(
                &meth_prob_pileup_column.haps[1],
                &next_meth_prob_pileup_column.haps[1],
            )
            || test_track(
                &meth_prob_pileup_column.haps[2],
                &next_meth_prob_pileup_column.haps[2],
            )
    }

    fn process_pileup_to_features(
        col: &HaplotypeMethProbPileupColumn,
        is_adjacent_cpg: bool,
    ) -> Vec<f32> {
        // For compatibility with previous python script there are two steps added below:
        // (1) ensure that any entries in cov() which did not have a methylation prob have a zero
        // entry added now:
        assert!(col.cov() >= col.meth_probs.len() as u32);
        let missing_prob_count = col.cov() - col.meth_probs.len() as u32;

        // (2) offset all ML probabilities to the beginning (instead of the center) of each bin:
        let prob_offset = (-1.0 / 512.0) as f32;

        let mut features = get_norm_prob_histogram(
            &col.meth_probs,
            MODEL_METH_PROB_BIN_COUNT,
            missing_prob_count,
            prob_offset,
        );
        features.push(is_adjacent_cpg as u32 as f32);

        assert_eq!(
            features.len(),
            MODEL_FEATURES_PER_CPG,
            "Computed feature count per CpG does not match expected value of {MODEL_FEATURES_PER_CPG}"
        );

        features
    }

    fn get_model_input_feature_volume(cpg_feature_ring: &FeatureRing) -> Vec<Vec<f32>> {
        assert_eq!(cpg_feature_ring.len(), MODEL_CPG_WINDOW_COUNT);

        let mut model_features = Vec::new();

        for cpg_index in 0..MODEL_CPG_WINDOW_COUNT {
            let cpg = cpg_feature_ring.get_item(cpg_index);
            assert_eq!(cpg.features.len(), MODEL_FEATURES_PER_CPG);
            model_features.push(cpg.features.clone());
        }
        model_features
    }

    #[cfg(not(feature = "tflite"))]
    fn get_model_prob(cpg_feature_ring: &FeatureRing, _: &mut TFLiteModelData) -> f32 {
        assert_eq!(cpg_feature_ring.len(), MODEL_CPG_WINDOW_COUNT);

        0f32
    }

    #[cfg(feature = "tflite")]
    fn get_model_prob(cpg_feature_ring: &FeatureRing, pileup_model: &mut TFLiteModelData) -> f32 {
        assert_eq!(cpg_feature_ring.len(), MODEL_CPG_WINDOW_COUNT);

        let input_tensor: &mut [f32] = pileup_model
            .interpreter
            .tensor_data_mut(pileup_model.input_index)
            .unwrap();

        let mut input_tensor_index = 0;
        for cpg_model_features in Self::get_model_input_feature_volume(cpg_feature_ring) {
            for feature_val in cpg_model_features {
                input_tensor[input_tensor_index] = feature_val;
                input_tensor_index += 1;
            }
        }

        pileup_model.interpreter.invoke().unwrap();

        let output_tensor: &[f32] = pileup_model
            .interpreter
            .tensor_data(pileup_model.output_index)
            .unwrap();

        output_tensor[0]
    }

    /// This continues the logic in `model_mode_process_track_pileup` for the portion after the CpG
    /// ring entry has been created. It is a separate method becuase this part is called separately
    /// from its original parent function
    ///
    fn model_mode_process_track_cpg_features(
        &mut self,
        track_id: usize,
        pileup_model: &mut Option<TFLiteModelData>,
        cpg_features: TrackCpGFeatures,
    ) {
        // Add the new CpG slot into the ring, discarding the last CpG entry in the window:
        let cpg_feature_ring = &mut self.cpg_feature_info.track_rings[track_id];

        cpg_feature_ring.push_item(cpg_features);

        // Check if we've acquired enough predictions yet to start making model calls
        if cpg_feature_ring.len() < MODEL_CPG_WINDOW_COUNT {
            return;
        }

        // Check if coverage of the call site is high enough to call the pileup prob
        let call_site_data = cpg_feature_ring.get_item(CALLED_CPG_INDEX);
        if call_site_data.cov < self.min_coverage {
            return;
        }

        // Check if the call is still part of this bin's region of responsibility
        if cpg_feature_ring.get_item(0).ref_pos >= self.end_pos {
            return;
        }

        if track_id == 0 {
            if let Some(chrom_training_sites) = self.chrom_training_sites {
                // Check if this is a target site:
                if let Some(meth_frac) = chrom_training_sites.sites.get(&call_site_data.ref_pos) {
                    // process cpg ring features into training model input feature volume
                    let features = Self::get_model_input_feature_volume(cpg_feature_ring);

                    let smi = SiteMethTrainingData {
                        features,
                        cov: call_site_data.cov,
                        meth_frac: *meth_frac as f32,
                    };

                    self.chrom_site_meth_training_data
                        .data
                        .insert(call_site_data.ref_pos, smi);
                }
            }
        }

        if let Some(pileup_model) = pileup_model {
            // process cpg ring feature volume into predicted pileup prob
            let pileup_prob = Self::get_model_prob(cpg_feature_ring, pileup_model);

            let smi = SiteMethInfo::Model {
                pileup_model_prob: pileup_prob,
                cov: call_site_data.cov,
            };

            self.site_meth_info.track_site_probs[track_id].insert(call_site_data.ref_pos, smi);
        }
    }

    /// Fully process single track pileup data for 'model' pileup mode
    ///
    /// Pileup data is processed into per-CPG model features and inserted into the cpg feature ring.
    ///
    /// As the CpG feature ring fills, each new complete CpG feature volume is processed to
    /// either/both of:
    /// 1. Make a new tflite model prediction, and write prediction to output bed track data
    /// 2. Writing model input features into output feature data
    ///
    #[allow(clippy::too_many_arguments)]
    fn model_mode_process_track_pileup(
        &mut self,
        ref_pos: i64,
        col: &HaplotypeMethProbPileupColumn,
        is_adjacent_cpg: bool,
        track_id: usize,
        pileup_model: &mut Option<TFLiteModelData>,
    ) {
        // Don't create CpG features for sites that don't have minimum model coverage
        //
        // Note this is separate from minimum output coverage, which we can apply just to the call
        // site
        if col.cov() < self.min_pileup_model_feature_coverage {
            return;
        }

        // Generate all info for the CpG site:
        let mut cpg = TrackCpGFeatures::new();
        cpg.ref_pos = ref_pos;
        cpg.cov = col.cov();

        cpg.features = Self::process_pileup_to_features(col, is_adjacent_cpg);

        self.model_mode_process_track_cpg_features(track_id, pileup_model, cpg);
    }

    fn count_mode_process_track_pileup(
        &mut self,
        ref_pos: i64,
        col: &HaplotypeMethProbPileupColumn,
        track_id: usize,
    ) {
        // Check if the call is still part of this bin's region of responsibility
        // TODO: is this needed for the count case?
        if ref_pos >= self.end_pos {
            return;
        }

        // Check if coverage of the call site is high enough to call the pileup prob
        if col.cov() < self.min_coverage {
            return;
        }

        self.site_meth_info.track_site_probs[track_id]
            .insert(ref_pos, col.get_count_mode_summary_stats());
    }

    /// Update site data
    ///
    /// Site data will be meth prob tracks and/or site meth training data depending on whether
    /// the method is running meth model or reporting training features.
    ///
    fn update_site_data(
        &mut self,
        ref_pos: i64,
        meth_prob_pileup_column: &MethProbPileupColumn,
        next_meth_prob_pileup_column: Option<&MethProbPileupColumn>,
        pileup_model: &mut Option<TFLiteModelData>,
    ) {
        if meth_prob_pileup_column.is_empty() {
            return;
        }

        let total_col;

        // Determine if this qualifies as a CpG site:
        match self.modsites_mode {
            ModSitesMode::denovo => {
                // Next column is required to determine if the next base is the 'G' of a CpG:
                let next_meth_prob_pileup_column = match next_meth_prob_pileup_column {
                    None => {
                        return;
                    }
                    Some(x) => x,
                };

                total_col = meth_prob_pileup_column.get_hap_total();
                let next_total_col = next_meth_prob_pileup_column.get_hap_total();

                if !Self::keep_cpg_site(
                    &total_col,
                    &next_total_col,
                    meth_prob_pileup_column,
                    next_meth_prob_pileup_column,
                ) {
                    return;
                }
            }
            ModSitesMode::reference => {
                if ref_pos < 0 {
                    return;
                }
                let ref_pos = ref_pos as usize;
                let chrom_ref = unwrap!(
                    self.chrom_ref,
                    "ref fasta must be defined in reference modsites mode"
                );
                if ref_pos + 1 >= chrom_ref.len()
                    || chrom_ref[ref_pos] != b'C'
                    || chrom_ref[ref_pos + 1] != b'G'
                {
                    return;
                }
                total_col = meth_prob_pileup_column.get_hap_total();
            }
        };

        match self.pileup_mode {
            MethReadsProcessorPileupMode::model => {
                // Determine if this CpG is adjacent to the previous CpG
                let is_adjacent_cpg = match self.last_cpg_ref_pos {
                    None => false,
                    Some(pos) => ref_pos - pos == 2,
                };
                self.last_cpg_ref_pos = Some(ref_pos);

                // process new CpG site into model features for each track:
                self.model_mode_process_track_pileup(
                    ref_pos,
                    &total_col,
                    is_adjacent_cpg,
                    0,
                    pileup_model,
                );
                self.model_mode_process_track_pileup(
                    ref_pos,
                    &meth_prob_pileup_column.haps[1],
                    is_adjacent_cpg,
                    1,
                    pileup_model,
                );
                self.model_mode_process_track_pileup(
                    ref_pos,
                    &meth_prob_pileup_column.haps[2],
                    is_adjacent_cpg,
                    2,
                    pileup_model,
                );
            }
            MethReadsProcessorPileupMode::count => {
                self.count_mode_process_track_pileup(ref_pos, &total_col, 0);
                self.count_mode_process_track_pileup(ref_pos, &meth_prob_pileup_column.haps[1], 1);
                self.count_mode_process_track_pileup(ref_pos, &meth_prob_pileup_column.haps[2], 2);
            }
        }
    }
}

/// Buffer to accumulate position specific basecount and methylation data
///
/// Ring buffer scheme allows array-like position access while keeping total memory demand and small
/// chunk allocations down
///
struct MethProbPileupRingBuffer {
    data: Vec<MethProbPileupColumn>,
    head_index: usize,

    /// The lowest position indexed in the buffer, this value is in head_index in the array
    head_pos: i64,

    /// The highest position with data stored in the buffer
    max_pos: i64,
}

impl MethProbPileupRingBuffer {
    pub fn new(start_pos: i64) -> Self {
        Self {
            data: vec![MethProbPileupColumn::new(); 100_000],
            head_index: 0,
            head_pos: start_pos,
            max_pos: start_pos - 1,
        }
    }

    pub fn head_pos(&self) -> i64 {
        self.head_pos
    }

    pub fn max_pos(&self) -> i64 {
        self.max_pos
    }

    fn get_index(&self, pos: i64) -> usize {
        assert!(pos >= self.head_pos);
        let pos_offset = (pos - self.head_pos) as usize;
        let dsize = self.data.len();
        assert!(pos_offset < dsize);
        (self.head_index + pos_offset) % dsize
    }

    pub fn get(&self, pos: i64) -> &MethProbPileupColumn {
        let index = self.get_index(pos);
        &self.data[index]
    }

    pub fn get_mut(&mut self, pos: i64) -> &mut MethProbPileupColumn {
        assert!(pos >= self.head_pos);
        let pos_offset = (pos - self.head_pos) as usize;
        loop {
            let dsize = self.data.len();
            if pos_offset < dsize {
                break;
            }
            self.increase_capacity();
        }
        let dsize = self.data.len();
        let index = (self.head_index + pos_offset) % dsize;
        if pos > self.max_pos {
            self.max_pos = pos;
        }
        &mut self.data[index]
    }

    pub fn reset_head_pos(&mut self, new_head_pos: i64) {
        if new_head_pos <= self.head_pos {
            return;
        }
        let pos_offset = (new_head_pos - self.head_pos) as usize;
        let dsize = self.data.len();

        if pos_offset >= dsize {
            // reset the whole container:
            self.data.iter_mut().for_each(|x| x.clear());
            self.head_index = 0;
        } else {
            for index_offset in 0..pos_offset {
                let index = (self.head_index + index_offset) % dsize;
                self.data[index].clear();
            }
            self.head_index = (self.head_index + pos_offset) % dsize;
        }
        self.head_pos = new_head_pos;
    }

    // double the capacity of the buffer
    fn increase_capacity(&mut self) {
        let old_dsize = self.data.len();
        self.data
            .append(&mut vec![MethProbPileupColumn::new(); old_dsize]);
        for old_index in 0..self.head_index {
            self.data.swap(old_index, old_index + old_dsize);
        }
    }
}

#[derive(Debug)]
pub struct ProcessBamRecordContext;

impl fmt::Display for ProcessBamRecordContext {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt.write_str("Bam record parsing error")
    }
}

impl Context for ProcessBamRecordContext {}

/// Pileup mode options
///
/// Used to select which method the MRP uses to summarize each pileup column as a site probability
///
/// Note that enum values and their doc comments appear in the cmdline help.
///
#[derive(Copy, Clone, Default, PartialEq, Eq, PartialOrd, Ord, strum::Display, ValueEnum)]
#[allow(non_camel_case_types)]
pub enum MethReadsProcessorPileupMode {
    /// Deep learning model-based approach
    #[default]
    model,

    /// Simple count-based approach
    count,
}

/// Modsites mode options
///
/// Used to select which CpG sites are reported in the output.
///
/// Note that enum values and their doc comments appear in the cmdline help.
///
#[derive(Copy, Clone, Default, PartialEq, Eq, PartialOrd, Ord, strum::Display, ValueEnum)]
#[allow(non_camel_case_types)]
pub enum ModSitesMode {
    /// Pick sites where "CG" is present in the majority of the reads
    #[default]
    denovo,

    /// Pick sites where "CG" is present in the reference sequence
    reference,
}

#[derive(Clone)]
pub struct MethReadsProcessorOptions {
    /// Methylation probabilities will only be generated and output for a site with at least this
    /// much coverage
    pub min_coverage: u32,

    /// Minimum mapq for reads to be used in the pileup
    pub min_mapq: u32,

    /// 2-letter SAM aux tag used for per-read haplotype ids
    pub hap_tag: [u8; 2],

    pub pileup_mode: MethReadsProcessorPileupMode,

    pub modsites_mode: ModSitesMode,

    /// This defines how large of a chromosome segment should be processed by a single thread
    pub segment_size: u64,
}

#[derive(Debug)]
pub struct ProcessMethInfoContext;

impl fmt::Display for ProcessMethInfoContext {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt.write_str("Error processing bam record methylation data")
    }
}

impl Context for ProcessMethInfoContext {}

type ReadToRefWithMeth = Vec<(Option<i64>, Option<f32>)>;

/// Enhance the read_to_ref map with methylation data at corresponding reference positions
///
fn get_read_to_ref_with_meth(
    read_to_ref: &[Option<i64>],
    cpg_meth_info: &Result<CpgMethInfo, i32>,
) -> error_stack::Result<ReadToRefWithMeth, ProcessMethInfoContext> {
    let mut result = Vec::new();
    let mut r2riter = read_to_ref.iter();
    if let Ok(cpg_meth_info) = cpg_meth_info {
        for (meth_index, (&read_pos, &prob)) in cpg_meth_info.pos_prob.iter().enumerate() {
            let mut get_next_ref_pos =
                || -> error_stack::Result<Option<i64>, ProcessMethInfoContext> {
                    match r2riter.next() {
                        Some(&x) => Ok(x),
                        None => {
                            Err(Report::new(ProcessMethInfoContext).attach_printable(format!(
                            "Failed to match methylation index {meth_index} corresponding to read pos {read_pos} with meth prob {prob}."
                        )))
                        }
                    }
                };
            while result.len() < read_pos {
                result.push((get_next_ref_pos()?, None));
            }
            result.push((get_next_ref_pos()?, Some(prob)));
        }
    }
    // map any remaining read_to_Ref iterator contents onto the result vector:
    result.extend(r2riter.map(|&x| (x, None)));
    Ok(result)
}

/// Intermediate data structures required to produce methylation probability tracks
///
pub struct MethReadsProcessor<'a> {
    /// Buffer to accumulate position specific basecount and methylation data
    meth_prob_pileup_ring_buffer: MethProbPileupRingBuffer,

    /// Position range over which bam input will be processed
    bam_process_range: std::ops::Range<i64>,

    /// Tracks the position of previous bam read input to insure that monotonic increase assumption
    /// is followed, and to trigger buffer processing when moving to a higher start_pos
    last_start: Option<i64>,

    /// Results processed from the data after completing site pileups
    post_pileup: PostPileupData<'a>,

    client_options: MethReadsProcessorOptions,
}

impl<'a> MethReadsProcessor<'a> {
    pub fn new(
        start: i64,
        end: i64,
        mrp_options: &MethReadsProcessorOptions,
        chrom_training_sites: Option<&'a ChromSiteData>,
        chrom_ref: Option<&'a Vec<u8>>,
    ) -> Self {
        assert!(end >= start);
        let is_chrom_start = start == 0;
        Self {
            meth_prob_pileup_ring_buffer: MethProbPileupRingBuffer::new(start),
            bam_process_range: std::ops::Range { start, end },
            last_start: None,
            post_pileup: PostPileupData::new(
                is_chrom_start,
                end,
                mrp_options.min_coverage,
                mrp_options.pileup_mode,
                mrp_options.modsites_mode,
                chrom_training_sites,
                chrom_ref,
            ),
            client_options: mrp_options.clone(),
        }
    }

    fn get_hap_id(record: &bam::Record, haplotype_tag: &[u8]) -> u8 {
        match get_optional_int_aux_tag(record, haplotype_tag) {
            Some(x) => {
                if x == 1 || x == 2 {
                    x as u8
                } else {
                    0u8
                }
            }
            None => 0u8,
        }
    }

    /// Find start and end entries in read_to_ref accounting for trimming:
    fn get_index_range(read_to_ref: &[Option<i64>], trim_edge_matches: i32) -> (usize, usize) {
        // find starting entry in read_to_ref accounting for trimming:
        let mut start_index = {
            let mut index = 0;
            let mut matches = 0;
            loop {
                if matches >= trim_edge_matches || index >= read_to_ref.len() {
                    break;
                }
                if read_to_ref[index].is_some() {
                    matches += 1;
                }
                index += 1;
            }
            index
        };

        // find stop entry in read_to_ref accounting for trimming:
        let mut end_index = {
            let mut index = read_to_ref.len();
            let mut matches = 0;
            loop {
                if matches >= trim_edge_matches || index == 0 {
                    break;
                }
                index -= 1;
                if read_to_ref[index].is_some() {
                    matches += 1;
                }
            }
            index
        };

        if start_index >= end_index {
            start_index = 0;
            end_index = 0;
        }

        (start_index, end_index)
    }

    /// Process a single bam record.
    ///
    /// It is assumed that this method will be called on bam records with monotonically
    /// increasing start positions, which should be fulfilled by scanning reads from
    /// a position sorted bam file.
    ///
    /// Return false when finished processing assigned region
    ///
    pub fn process_bam_record(
        &mut self,
        record: &bam::Record,
        pileup_model: &mut Option<TFLiteModelData>,
    ) -> error_stack::Result<bool, ProcessBamRecordContext> {
        if record.mapq() < self.client_options.min_mapq as u8 {
            return Ok(true);
        }

        let start = record.pos();

        if let Some(val) = self.last_start {
            assert!(start >= val);
        }
        self.last_start = Some(start);

        // Process and remove entries in the pileup structure with key positions less than start-1
        //
        // This value is start - 1 so that we can check the G base after the CpG when determining
        // whether a site is going to be kept.
        let process_pileup_buffer_up_to_pos = start - 1;
        for ref_pos in self.meth_prob_pileup_ring_buffer.head_pos()..process_pileup_buffer_up_to_pos
        {
            if ref_pos > self.meth_prob_pileup_ring_buffer.max_pos() {
                self.meth_prob_pileup_ring_buffer
                    .reset_head_pos(process_pileup_buffer_up_to_pos);
                break;
            }

            let meth_prob_pileup_column = self.meth_prob_pileup_ring_buffer.get(ref_pos);
            let next_meth_prob_pileup_column = self.meth_prob_pileup_ring_buffer.get(ref_pos + 1);
            self.post_pileup.update_site_data(
                ref_pos,
                meth_prob_pileup_column,
                Some(next_meth_prob_pileup_column),
                pileup_model,
            );
            self.meth_prob_pileup_ring_buffer
                .reset_head_pos(ref_pos + 1);
        }

        // Check scan termination conditions:
        if process_pileup_buffer_up_to_pos >= self.bam_process_range.end
            && !self.post_pileup.unprocessed_cpgs()
        {
            return Ok(false);
        }

        // Get CpG methylation info for this read
        let cpg_meth_info = decode_cpg_meth_info(record);

        // For compatibility with python script, skip further processing if MM or ML tag doesn't
        // exist on read:
        if let Err(1) = cpg_meth_info {
            return Ok(true);
        }

        // Get haplotype id for this read
        let hap_id = Self::get_hap_id(record, &self.client_options.hap_tag);

        // Process bam cigar to map read positions back to the reference positions
        let read_to_ref = get_complete_read_to_ref_pos_map(record);

        // This value is used for compatibility with the python script:
        let trim_edge_matches = 20;

        let (start_index, end_index) = Self::get_index_range(&read_to_ref, trim_edge_matches);

        let read_to_ref_with_meth = get_read_to_ref_with_meth(&read_to_ref, &cpg_meth_info)
            .change_context(ProcessBamRecordContext)?;

        // Iterate through read_to_ref_with_meth and add this to the pileup structure counts
        for (read_pos_offset, (ref_pos_option, meth_prob_option)) in read_to_ref_with_meth
            [start_index..end_index]
            .iter()
            .enumerate()
        {
            let ref_pos = match ref_pos_option {
                Some(x) => *x,
                None => {
                    continue;
                }
            };

            let read_pos = read_pos_offset + start_index;

            // Check that the ref position is in the assigned reference segment for this object
            if ref_pos < self.bam_process_range.start {
                continue;
            }

            let hap_pileup =
                &mut self.meth_prob_pileup_ring_buffer.get_mut(ref_pos).haps[hap_id as usize];

            hap_pileup.add_base(record.seq()[read_pos]);
            if let Some(meth_prob) = meth_prob_option {
                hap_pileup.meth_probs.push(*meth_prob);
            }
        }

        Ok(true)
    }

    /// Process any buffered data remaining in the object, then return the summarized results
    ///
    /// Depending on run options, output includes site methylation info and/or sitetraining features
    ///
    pub fn complete_processing(
        &mut self,
        pileup_model: &mut Option<TFLiteModelData>,
    ) -> (ChromSiteMethInfo, ChromSiteMethTrainingData) {
        for ref_pos in self.meth_prob_pileup_ring_buffer.head_pos()
            ..self.meth_prob_pileup_ring_buffer.max_pos() + 1
        {
            let meth_prob_pileup_column = self.meth_prob_pileup_ring_buffer.get(ref_pos);
            let next_meth_prob_pileup_column = self.meth_prob_pileup_ring_buffer.get(ref_pos + 1);
            self.post_pileup.update_site_data(
                ref_pos,
                meth_prob_pileup_column,
                Some(next_meth_prob_pileup_column),
                pileup_model,
            );
            self.meth_prob_pileup_ring_buffer
                .reset_head_pos(ref_pos + 1);
        }

        if let MethReadsProcessorPileupMode::model = self.post_pileup.pileup_mode {
            // If in model pileup mode, finish processing any remaining data in the CpG ring buffer:
            //
            let trailing_window_position_count = MODEL_CPG_WINDOW_COUNT - (CALLED_CPG_INDEX + 1);
            for track_id in 0..HAPLOTYPE_TRACK_COUNT {
                // Don't bother checking for specific conditions, just shove the blanks in and the
                // process logic will automatically determine if its a qualifying condition:
                for _ in 0..trailing_window_position_count {
                    self.post_pileup.model_mode_process_track_cpg_features(
                        track_id,
                        pileup_model,
                        TrackCpGFeatures::blank(),
                    );
                }
            }
        }

        (
            std::mem::take(&mut self.post_pileup.site_meth_info),
            std::mem::take(&mut self.post_pileup.chrom_site_meth_training_data),
        )
    }
}
