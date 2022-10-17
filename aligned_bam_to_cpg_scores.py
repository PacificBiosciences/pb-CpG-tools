#!/usr/bin/env python
# coding: utf-8
import argparse
import concurrent.futures
import logging
import numpy as np
import pandas as pd
import pyBigWig
import pysam
import os
import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
from numpy.lib.stride_tricks import sliding_window_view
from operator import itemgetter
from tqdm import tqdm

os.environ["CUDA_VISIBLE_DEVICES"] = ""


def get_args():
    """
    Get arguments from command line with argparse.
    """
    parser = argparse.ArgumentParser(
        prog='aligned_bam_to_cpg_scores.py',
        description="""Calculate CpG positions and scores from an aligned bam file. Outputs raw and 
        coverage-filtered results in bed and bigwig format, including haplotype-specific results (when available).""")

    parser.add_argument("-b", "--bam",
                        required=True,
                        metavar="input.bam",
                        help="The aligned BAM file. This file must be sorted and indexed.")
    parser.add_argument("-f", "--fasta",
                        required=True,
                        metavar="ref.fasta",
                        help="The reference fasta file.")
    parser.add_argument("-o", "--output_label",
                        required=True,
                        metavar="label",
                        help="Label for output files, which results in [label].bed/bw.")
    parser.add_argument("-p", "--pileup_mode",
                        required=False,
                        choices=["model", "count"],
                        default="model",
                        help="Use a model-based approach to score modifications across sites (model) "
                             "or a simple count-based approach (count). [default = %(default)s]")
    parser.add_argument("-d", "--model_dir",
                        required=False,
                        default=None,
                        metavar="/path/to/model/dir",
                        help="Full path to the directory containing the model (*.pb files) to load. [default = None]")
    parser.add_argument("-m", "--modsites",
                        required=False,
                        choices=["denovo", "reference"],
                        default="denovo",
                        help="Only output CG sites with a modification probability > 0 "
                             "(denovo), or output all CG sites based on the "
                             "supplied reference fasta (reference). [default = %(default)s]")
    parser.add_argument("-c", "--min_coverage",
                        required=False,
                        default=4,
                        type=int,
                        metavar="int",
                        help="Minimum coverage required for filtered outputs. [default: %(default)d]")
    parser.add_argument("-q", "--min_mapq",
                        required=False,
                        default=0,
                        type=int,
                        metavar="int",
                        help="Ignore alignments with MAPQ < N. [default: %(default)d]")
    parser.add_argument("-a", "--hap_tag",
                        required=False,
                        default="HP",
                        metavar="TAG",
                        help="The SAM tag containing haplotype information. [default: %(default)s]")
    parser.add_argument("-s", "--chunksize",
                        required=False,
                        default=500000,
                        type=int,
                        metavar="int",
                        help="Break reference regions into chunks "
                             "of this size for parallel processing. [default = %(default)d]")
    parser.add_argument("-t", "--threads",
                        required=False,
                        default=1,
                        type=int,
                        metavar="int",
                        help="Number of threads for parallel processing. [default = %(default)d]")

    return parser.parse_args()


def validate_args(args):

    def error_exit(msg):
        raise Exception(msg)

    def check_required_file(file, label):
        if not os.path.isfile(file):
            error_exit(f"Can't find {label} file '{file}'")

    check_required_file(args.bam, "input bam")

    def is_bam_index_found(bam_file):
        bam_index_extensions = (".bai", ".csi")
        for ext in bam_index_extensions:
            bam_index_file=bam_file+ext
            if os.path.isfile(bam_index_file):
                return True
        return False

    if not is_bam_index_found(args.bam):
        error_exit(f"Can't find index for bam file '{args.bam}'")

    check_required_file(args.fasta, "reference fasta")

    if args.pileup_mode == "model":
        if args.model_dir is None:
            error_exit("Must supply a model to use when running model-based scoring")
        else:
            if not os.path.isdir(args.model_dir):
                error_exit("{} is not a valid directory path!".format(args.model_dir))
    else:
        if args.model_dir is not None:
            error_exit("Model directory is not used unless model-based scoring is selected")


def setup_logging(output_label):
    """
    Set up logging to file.
    """
    logname = "{}-aligned_bam_to_cpg_scores.log".format(output_label)
    # ensure logging file does not exist, if so remove
    if os.path.exists(logname):
        os.remove(logname)

    # set up logging to file
    logging.basicConfig(filename=logname,
                        format="%(asctime)s: %(levelname)s: %(message)s",
                        datefmt='%d-%b-%y %H:%M:%S',
                        level=logging.DEBUG)


def log_args(args):
    """
    Record argument settings in log file.
    """
    logging.info("Using following argument settings:")
    for arg, val in vars(args).items():
        logging.info("\t--{}: {}".format(arg, val))


def get_regions_to_process(input_bam, input_fasta, chunksize, modsites, pileup_mode, model_dir, min_mapq, hap_tag):
    """
    Breaks reference regions into smaller regions based on chunk
    size specified. Returns a list of lists that can be used for
    multiprocessing. Each sublist contains:
    [bam path (str), fasta path (str), modsites (str),
    reference name (str), start coordinate (int), stop coordinate (int)]

    :param input_bam: Path to input bam file. (str)
    :param input_fasta: Path to reference fasta file. (str)
    :param chunksize: Chunk size (default = 500000). (int)
    :param modsites: Filtering method. (str: "denovo", "reference")
    :param pileup_mode: Site modification calling method. (str: "model", "count")
    :param model_dir: Full path to model directory to load (if supplied), otherwise is None.
    :param min_mapq: Minimum mapping quality score. (int)
    :param hap_tag: The SAM tag label containing haplotype information. (str)
    :return regions_to_process: List of lists containing region sizes. (list)
    """
    logging.info("get_regions_to_process: Starting chunking.")
    # open the input bam file with pysam
    bamIn = pysam.AlignmentFile(input_bam, 'rb')
    # empty list to store sublists with region information
    regions_to_process = []
    # iterate over reference names and their corresponding lengths
    references = zip(bamIn.references, bamIn.lengths)
    for ref, length in references:
        start = 1
        while start < length:
            end = start + chunksize
            if end < length:
                regions_to_process.append(
                    [input_bam, input_fasta, modsites, pileup_mode, model_dir, ref, start, end - 1, min_mapq, hap_tag])
            else:
                regions_to_process.append(
                    [input_bam, input_fasta, modsites, pileup_mode, model_dir, ref, start, length, min_mapq, hap_tag])
            start = start + chunksize
    # close bam
    bamIn.close()
    logging.info("get_regions_to_process: Created {:,} region chunks.\n".format(len(regions_to_process)))

    return regions_to_process


def cg_sites_from_fasta(input_fasta, ref):
    """
    Gets all CG site positions from a given reference region, and
    make positions keys in a dict with empty strings as vals.

    :param input_fasta: A path to reference fasta file. (str)
    :param ref: Reference name. (str)
    :return cg_sites_ref_set: Set with all CG ref positions. (set)
    """
    # open fasta with BioPython and iterated over records
    cg_sites_ref_set = set()

    with open(input_fasta) as fh:
        for record in SeqIO.parse(fh, "fasta"):
            # if record name matches this particular ref,
            if record.id == ref:
                # use regex to find all indices for 'CG' in the reference seq, e.g. the C positions
                cg_sites_ref_set = {i.start() for i in re.finditer('CG', str(record.seq.upper()))}
                # there may be some stretches without any CpGs in a reference region
                # handle these edge cases by adding a dummy value of -1 (an impossible coordinate)
                if not cg_sites_ref_set:
                    cg_sites_ref_set.add(-1)
                # once seq is found, stop iterating
                break
    # make sure the ref region was matched to a ref fasta seq
    if not cg_sites_ref_set:
        logging.error("cg_sites_from_fasta: The sequence '{}' was not found in the reference fasta file.".format(ref))
        raise ValueError('The sequence "{}" was not found in the reference fasta file!'.format(ref))

    return cg_sites_ref_set


def get_mod_sequence(integers):
    """
    A generator that takes an iterable of integers coding mod bases from the SAM Mm tags, and yields an iterable of
    positions of sequential bases.
    Example: [5, 12, 0] -> [6, 19, 20]
    In above example the 6th C, 19th C, and 20th C are modified
    See this example described in: https://samtools.github.io/hts-specs/SAMtags.pdf; Dec 9 2021

    :param integers: Iterable of integers (parsed from SAM Mm tag). (iter)
    :return mod_sequence: Iterator of integers, 1-based counts of position of modified base in set of bases. (iter)
    """
    base_count = 0
    for i in integers:
        base_count += i + 1
        yield base_count


def get_base_indices(query_seq, base, reverse):
    """
    Find all occurrences of base in query sequence and make a list of their
    indices. Return the list of indices.

    :param query_seq: The original read sequence (not aligned read sequence). (str)
    :param base: The nucleotide modifications occur on ('C'). (str)
    :param reverse: True/False whether sequence is reversed. (Boolean)
    :return: List of integers, 0-based indices of all bases in query seq. (list)
    """
    if reverse == False:
        return [i.start() for i in re.finditer(base, query_seq)]
    # if seq stored in reverse, need reverse complement to get correct indices for base
    # use biopython for this (convert to Seq, get RC, convert to string)
    else:
        return [i.start() for i in re.finditer(base, str(Seq(query_seq).reverse_complement()))]


def parse_mmtag(query_seq, mmtag, modcode, base, reverse):
    """
    Get a list of the 0-based indices of the modified bases in the query sequence.

    :param query_seq: The original read sequence (not aligned read sequence). (str)
    :param mmtag: The Mm tag obtained for the read ('C+m,5,12,0;'). (str)
    :param modcode: The modification code to search for in the tag ('C+m'). (str)
    :param base: The nucleotide modifications occur on ('C'). (str)
    :param reverse: True/False whether sequence is reversed. (Boolean)
    :return mod_base_indices: List of integers, 0-based indices of all mod bases in query seq. (list)
    """

    def get_modline():
        """
        tags are written as: C+m,5,12,0;C+h,5,12,0;
        if multiple mod types present in tag, must find relevant one first
        """
        for mmsection in mmtag.split(';'):
            if not mmsection.startswith(modcode):
                continue
            start_index = len(modcode)

            # Note that this method treats either value of the skip-base mode in the same way to
            # stay back-compatible with older bams where this wasn't specified.
            if len(mmsection) > start_index and mmsection[start_index] in "?.":
                start_index += 1
            if len(mmsection) > start_index:
                if mmsection[start_index] != ',':
                    raise Exception(f"Can't parse MM tag segment '{mmsection}' due to unexpected value in position {start_index+1}")
                start_index += 1
            return mmsection[start_index:]
        return None

    modline = get_modline()
    if modline is None or modline == "":
        return []

    # Get the sequence of the mod bases from tag integers
    # this is a 1-based position of each mod base in the complete set of this base from this read
    # e.g., [6, 19, 20] = the 6th, 19th, and 20th C bases are modified in the set of Cs
    mod_sequence = get_mod_sequence((int(x) for x in modline.split(',')))

    # Get all 0-based indices of this base in this read, e.g. every C position
    base_indices = get_base_indices(query_seq, base, reverse)

    # Use the mod sequence to identify indices of the mod bases in the read
    mod_base_indices = []
    base_indices_size = len(base_indices)
    for i in mod_sequence:
        if i <= 0 or i > base_indices_size:
            raise Exception(
                f"Base modification offsets in MM tag for modification type '{modcode}' are inconsistent with read length")
        mod_base_indices.append(base_indices[i - 1])

    return mod_base_indices


def parse_mltag(mltag):
    """
    Convert 255 discrete integer code into mod score 0-1, return as a list.

    This is NOT designed to handle interleaved Ml format for multiple mod types!

    :param mltag: The Ml tag obtained for the read with('Ml:B:C,204,89,26'). (str)
    :return: List of floats, probabilities of all mod bases in query seq. (list)
    """
    return [round(x / 256, 3) if x > 0 else 0 for x in mltag]


def get_mod_dict(query_seq, mmtag, modcode, base, mltag, reverse):
    """
    Make a dictionary from the Mm and Ml tags, in which the
    modified base index (in the query seq) is the key and the
    mod score is the value.

    This is NOT designed to handle interleaved Ml format for multiple mod types!

    :param query_seq: The original read sequence (not aligned read sequence). (str)
    :param mmtag: The Mm tag obtained for the read ('C+m,5,12,0;'). (str)
    :param modcode: The modification code to search for in the tag ('C+m'). (str)
    :param base: The nucleotide modifications occur on ('C'). (str)
    :param mltag: The Ml tag obtained for the read with('Ml:B:C,204,89,26'). (str)
    :param reverse: True/False whether sequence is reversed. (Boolean)
    :return mod_dict: Dictionary with mod positions and scores. (dict)
    """
    mod_base_indices = parse_mmtag(query_seq, mmtag, modcode, base, reverse)
    mod_scores = parse_mltag(mltag)
    if len(mod_base_indices) != len(mod_scores):
        raise Exception("Can't resolve base modifications from MM and ML tags in read")

    mod_dict = dict(zip(mod_base_indices, mod_scores))

    return mod_dict


def get_mm_and_ml_values(read):
    """
    Check for SAM-spec methylation tags in reads

    Return a 2-tuple of MM and ML strings, or None if these tags do not exist
    """

    mmtag, mltag = None, None
    try:
        mmtag = read.get_tag('Mm')
        mltag = read.get_tag('Ml')
    except KeyError:
        pass
    try:
        mmtag = read.get_tag('MM')
        mltag = read.get_tag('ML')
    except KeyError:
        pass

    if mmtag is None or mltag is None:
        # if no SAM-spec methylation tags present, ignore read and log
        logging.warning("pileup_from_reads: read missing MM and/or ML tag(s): {}".format(read.query_name))
        return None
    else:
        return mmtag, mltag

    # Filtering on empty MM or ML tag values makes sense but will change the pileup results, so we need to test impact
    # first.
    #
    # elif mmtag == "" or mltag == "":
    #    logging.warning("pileup_from_reads: MM and/or ML tag(s) have no value: {}".format(read.query_name))
    #    return None


class PileupData:
    def __init__(self):
        self.basemod_data = []

        # These structures are only used for modsites denovo mode
        self.pos_pileup = []
        self.pos_pileup_hap1 = []
        self.pos_pileup_hap2 = []


def process_read(ref, pos_start, pos_stop, hap_tag, is_denovo_modsites, pileup_data, read):
    """
    Process pileup information from a single read

    :param ref: Reference name. (str)
    :param pos_start: Start coordinate for region. (int)
    :param pos_stop: Stop coordinate for region. (int)
    :param hap_tag: Name of SAM tag containing haplotype information. (str)
    :param is_denovo_modsites: True if modsites mode is 'denovo' (bool)
    :param pileup_data:
    :param read:
    """

    # identify the haplotype tag, if any (default tag = HP)
    # values are 1 or 2 (for haplotypes), or 0 (no haplotype)
    # an integer is expected but custom tags can produce strings instead
    try:
        hap_val = read.get_tag(hap_tag)
        try:
            hap = int(hap_val)
        except ValueError:
            logging.error(
                "coordinates {}: {:,}-{:,}: (2) pileup_from_reads: illegal haplotype value {}".format(ref, pos_start,
                                                                                                      pos_stop,
                                                                                                      hap_val))
    except KeyError:
        hap = 0

    meth_tags = get_mm_and_ml_values(read)
    if meth_tags is None:
        return

    mmtag, mltag = meth_tags

    if not pileup_data.basemod_data:
        ref_pos_count = 1 + pos_stop - pos_start

        pileup_data.basemod_data = [[] for _ in range(ref_pos_count)]
        if is_denovo_modsites:
            pileup_data.pos_pileup = [[] for _ in range(ref_pos_count)]
            pileup_data.pos_pileup_hap1 = [[] for _ in range(ref_pos_count)]
            pileup_data.pos_pileup_hap2 = [[] for _ in range(ref_pos_count)]

    is_reverse = bool(read.is_reverse)
    strand = "+"
    if is_reverse:
        strand = "-"
    rev_strand_offset = len(read.query_sequence) - 2

    # note that this could potentially be used for other mod types, but
    # the Mm and Ml parsing functions are not set up for the interleaved format
    # e.g.,  ‘Mm:Z:C+mh,5,12; Ml:B:C,204,26,89,130’ does NOT work
    # to work it must be one mod type, and one score per mod position
    mod_dict = get_mod_dict(read.query_sequence, mmtag, 'C+m', 'C', mltag, is_reverse)

    # iterate over positions
    for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True)[20:-20]:
        # make sure ref position is in range of ref target region
        if ref_pos >= pos_start and ref_pos <= pos_stop:
            ref_offset = ref_pos - pos_start

            # building a consensus is MUCH faster when we iterate over reads (vs. by column then by read)
            # we are building a dictionary with ref position as key and list of bases as val
            if is_denovo_modsites:
                query_base = read.query_sequence[query_pos]
                pileup_data.pos_pileup[ref_offset].append(query_base)
                if hap == 1:
                    pileup_data.pos_pileup_hap1[ref_offset].append(query_base)
                elif hap == 2:
                    pileup_data.pos_pileup_hap2[ref_offset].append(query_base)

            # identify if read is reverse strand or forward to set correct location
            if is_reverse:
                location = (rev_strand_offset - query_pos)
            else:
                location = query_pos

            # check if this position has a mod score in the dictionary,
            # if not assign score of zero
            score = mod_dict.get(location, 0)

            # Add tuple with strand, modification score, and haplotype to the list for this position
            pileup_data.basemod_data[ref_offset].append((strand, score, hap))


def pileup_from_reads(bamIn, ref, pos_start, pos_stop, min_mapq, hap_tag, modsites):
    """
    For a given region, retrieve all reads.
    For each read, iterate over positions aligned to this region.
    Build a list with an entry for each ref position in the region. Each entry has a list of 3-tuples, each of which
    includes information from a read base read aligned to that site. The 3-tuple contains strand information,
    modification score, and haplotype.
    (strand symbol (str), mod score (float), haplotype (int))

    Return the unfiltered list of base modification data.

    :param bamIn: AlignmentFile object of input bam file.
    :param ref: Reference name. (str)
    :param pos_start: Start coordinate for region. (int)
    :param pos_stop: Stop coordinate for region. (int)
    :param min_mapq: Minimum mapping quality score. (int)
    :param hap_tag: Name of SAM tag containing haplotype information. (str)
    :param modsites: Filtering method. (str: "denovo", "reference")
    :return basemod_data: Unfiltered list of base modification data (list)
    :return cg_sites_read_set: Set of positions in read consensus sequence with CG, given as reference position. The
     set is empty unless modsites is 'denovo' (set)
    """
    logging.debug("coordinates {}: {:,}-{:,}: (2) pileup_from_reads".format(ref, pos_start, pos_stop))

    pileup_data = PileupData()

    is_denovo_modsites = modsites == "denovo"

    # iterate over all reads present in this region
    for read in bamIn.fetch(contig=ref, start=pos_start, stop=pos_stop):
        # check if passes minimum mapping quality score
        if read.mapping_quality < min_mapq:
            # logging.warning("pileup_from_reads: read did not pass minimum mapQV: {}".format(read.query_name))
            continue

        try:
            process_read(ref, pos_start, pos_stop, hap_tag, is_denovo_modsites, pileup_data, read)
        except Exception as e:
            raise Exception("Exception thrown while processing read {}: {}\n".format(read.query_name, e))

    cg_sites_read_set = set()
    if is_denovo_modsites:
        for refpos_list in (pileup_data.pos_pileup, pileup_data.pos_pileup_hap1, pileup_data.pos_pileup_hap2):
            last_base = 'N'
            last_index = 0
            for index, v in enumerate(refpos_list):
                # find the most common base, if no reads present use N
                if len(v):
                    base = Counter(v).most_common(1)[0][0]
                else:
                    base = 'N'
                if last_base == 'C' and base == 'G':
                    cg_sites_read_set.add(pos_start + last_index)

                # This restriction recreates the original code behavior:
                # - Advantage: Method can find a CpG aligning across a deletion in the reference
                # - Disadvantage: Method will find 'fake' CpG across gaps in the haplotype phasing
                #
                # The disadvantage is fixable, but first focus on identical output to make verification easy
                if base != 'N':
                    last_base = base
                    last_index = index

    return pileup_data.basemod_data, cg_sites_read_set


def filter_basemod_data(basemod_data, cg_sites_read_set, ref, pos_start, pos_stop, input_fasta, modsites):
    """
    Filter the per-position base modification data, based on the modsites option selected:
    "reference": Keep all sites that match a reference CG site (this includes both
                 modified and unmodified sites). It will exclude all modified sites
                 that are not CG sites, according to the ref sequence.
    "denovo": Keep all sites which have at least one modification score > 0, per strand.
              This can include sites that are CG in the reads, but not in the reference.
              It can exclude CG sites with no modifications on either strand from being
              written to the bed file.

    Return the filtered list.

    :param basemod_data: List of base modification data per position, offset by pos_start. (list)
    :param cg_sites_read_set: Set with reference coordinates for all CG sites in consensus from reads. (set)
    :param ref: A path to reference fasta file. (str)
    :param pos_start: Start coordinate for region. (int)
    :param pos_stop: Stop coordinate for region. (int)
    :param modsites: Filtering method. (str: "denovo", "reference")
    :param ref: Reference name. (str)
    :return filtered_basemod_data: List of 2-tuples for each position retained after filtering. Each 2-tuple is the
    reference position and base mod data list. The list is sorted by reference position (list)
    """
    filtered_basemod_data = []
    if modsites == "reference":
        if basemod_data:
            # Get CG positions in reference
            cg_sites_ref_set = cg_sites_from_fasta(input_fasta, ref)
            # Keep all sites that match a reference CG position and have at least one basemod observation.
            filtered_basemod_data = [(i + pos_start, v) for i, v in enumerate(basemod_data) if
                                     (i + pos_start) in cg_sites_ref_set and v]

        logging.debug(
            "coordinates {}: {:,}-{:,}: (3) filter_basemod_data: sites kept = {:,}".format(ref, pos_start, pos_stop,
                                                                                           len(filtered_basemod_data)))

    elif modsites == "denovo":
        if basemod_data:
            # Keep all sites that match position of a read consensus CG site.
            filtered_basemod_data = [(i + pos_start, v) for i, v in enumerate(basemod_data) if
                                     (i + pos_start) in cg_sites_read_set]

        logging.debug(
            "coordinates {}: {:,}-{:,}: (3) filter_basemod_data: sites kept = {:,}".format(ref, pos_start, pos_stop,
                                                                                           len(filtered_basemod_data)))

    del basemod_data
    del cg_sites_read_set

    return filtered_basemod_data


def calc_stats(df):
    """
    Gets summary stats from a given dataframe p.
    :param df: Pandas dataframe.
    :return: Summary statistics
    """
    total = df.shape[0]
    mod = df[df['prob'] > 0.5].shape[0]
    unMod = df[df['prob'] <= 0.5].shape[0]

    modScore = "." if mod == 0 else str(round(df[df['prob'] > 0.5]['prob'].mean(), 3))
    unModScore = "." if unMod == 0 else str(round(df[df['prob'] <= 0.5]['prob'].mean(), 3))
    percentMod = 0.0 if mod == 0 else round((mod / total) * 100, 1)

    return percentMod, mod, unMod, modScore, unModScore


def collect_bed_results_count(ref, pos_start, pos_stop, filtered_basemod_data):
    """
    Iterates over reference positions and for each position, makes a pandas dataframe from the sublists.
    The dataframe is filtered for strands and haplotypes, and summary statistics are
    calculated with calc_stats().
    For each position and strand/haploytpe combination, a sublist of summary information
    is appended to the bed_results list:
    [(0) ref name, (1) start coord, (2) stop coord, (3) mod probability, (4) haplotype, (5) coverage,
    (6) mod sites, (7) unmod sites, (8) mod score, (9) unmod score]
    This information is used to write the output bed file.

    :param ref: Reference name. (str)
    :param pos_start: Start coordinate for region. (int)
    :param pos_stop: Stop coordinate for region. (int)
    :param filtered_basemod_data: List of 2-tuples for each position remaining after filtration. Each 2-tuple is the
    reference position and base mod dat. The list is sorted by reference position (list)
    :return bed_results: List of sublists with information to write the output bed file. (list)
    """
    logging.debug("coordinates {}: {:,}-{:,}: (4) collect_bed_results_count".format(ref, pos_start, pos_stop))
    # intiate empty list to store bed sublists
    bed_results = []

    # iterate over the ref positions and corresponding vals
    for (refPosition, modinfoList) in filtered_basemod_data:
        # create pandas dataframe from this list of sublists
        df = pd.DataFrame(modinfoList, columns=['strand', 'prob', 'hap'])

        # Filter dataframe based on strand/haplotype combinations, get information,
        # and create sublists and append to bed_results.
        # merged strands / haplotype 1
        percentMod, mod, unMod, modScore, unModScore = calc_stats(df[df['hap'] == 1])
        if mod + unMod >= 1:
            bed_results.append([ref, refPosition, (refPosition + 1), percentMod,
                                "hap1", mod + unMod, mod, unMod, modScore, unModScore])

        # merged strands / haplotype 2
        percentMod, mod, unMod, modScore, unModScore = calc_stats(df[df['hap'] == 2])
        if mod + unMod >= 1:
            bed_results.append([ref, refPosition, (refPosition + 1), percentMod,
                                "hap2", mod + unMod, mod, unMod, modScore, unModScore])

        # merged strands / both haplotypes
        percentMod, mod, unMod, modScore, unModScore = calc_stats(df)
        if mod + unMod >= 1:
            bed_results.append([ref, refPosition, (refPosition + 1), percentMod,
                                "Total", mod + unMod, mod, unMod, modScore, unModScore])

    return bed_results


def get_normalized_histo(probs, adj):
    """
    Create the array data structure needed to apply the model, for a given site.

    :param probs: List of methylation probabilities. (list)
    :param adj: A 0 or 1 indicating whether previous position was a CG. (int)
    :return: List with normalized histogram and coverage (if min coverage met), else returns empty list. (list)
    """

    cov = len(probs)
    if (cov >= 4):
        hist = np.histogram(probs, bins=20, range=[0, 1])[0]
        norm = np.linalg.norm(hist)
        # divide hist by norm and add values to array
        # add either 0 (not adjacent to a prior CG) or 1 (adjacent to a prior CG) to final spot in array
        norm_hist = np.append(hist / norm, adj)
        return [norm_hist, cov]
    else:
        return []


def discretize_score(score, coverage):
    """
    Apply a small correction to the model probability to make it
    compatible with the number of reads at that site. Allows the number
    of modified and unmodified reads to be estimated.

    :param score: Modification probability, from model. (float)
    :param coverage: Number of reads. (int)
    :return mod_reads: Estimated number of modified reads. (int)
    :return unmod_reads: Estimated number of unmodified reads. (int)
    :return adjusted_score: Adjusted probability score, based on percent modified reads. (float)
    """
    # need to round up or round down modified read numbers based on score
    # which allows a push towards 0/50/100 for adjusted score
    if score > 50:
        if score < 65:
            mod_reads = int(np.floor(score / 100 * float(coverage)))
        else:
            mod_reads = int(np.ceil(score / 100 * float(coverage)))
    else:
        if score > 35:
            mod_reads = int(np.ceil(score / 100 * float(coverage)))
        else:
            mod_reads = int(np.floor(score / 100 * float(coverage)))

    unmod_reads = int(coverage) - mod_reads

    if mod_reads == 0:
        adjusted_score = 0.0
    else:
        adjusted_score = round((mod_reads / (mod_reads + unmod_reads)) * 100, 1)

    return mod_reads, unmod_reads, adjusted_score


def apply_model(refpositions, normhistos, coverages, ref, pos_start, pos_stop, model, hap, bed_results):
    """
    Apply model to make modification calls for all sites using a sliding window approach.
    Append to a list of results, ultimately for bed file:
        [(0) ref name, (1) start coord, (2) stop coord, (3) mod probability, (4) haplotype, (5) coverage,
        (6) mod sites, (7) unmod sites, (8) adjusted probability]
    :param refpositions: List with all CG positions. (list)
    :param normhistos: List with all normalized histogram data structures. (list)
    :param coverages: List with all CG coverages. (list)
    :param ref: Reference contig name. (str)
    :param pos_start: Start coordinate for region. (int)
    :param pos_stop: Stop coordinate for region. (int)
    :param model: The tensorflow model object.
    :param hap: Label of haplotype (hap1, hap2, or Total). (str)
    :param bed_results: List of bed results to which these model results will be appended (list)
    """
    if len(normhistos) > 11:
        featPad = np.pad(np.stack(normhistos), pad_width=((6, 4), (0, 0)), mode='constant', constant_values=0)
        featuresWindow = sliding_window_view(featPad, 11, axis=0)

        featuresWindow = np.swapaxes(featuresWindow, 1, 2)
        predict = model.predict(featuresWindow)

        predict = np.clip(predict, 0, 1)

        for i, position in enumerate(refpositions):
            model_score = round(predict[i][0] * 100, 1)
            mod_reads, unmod_reads, adjusted_score = discretize_score(model_score, coverages[i])
            bed_results.append(
                (ref, position, (position + 1), model_score, hap, coverages[i], mod_reads, unmod_reads, adjusted_score))
    else:
        logging.warning(
            "coordinates {}: {:,}-{:,}: apply_model: insufficient data for {}".format(ref, pos_start, pos_stop, hap))


def collect_bed_results_model(ref, pos_start, pos_stop, filtered_basemod_data, model_dir):
    """
    Iterates over reference positions and creates normalized histograms of scores,
    feeds all sites and scores into model function to assign modification probabilities,
    and creates a list of sublists for writing bed files:
        [(0) ref name, (1) start coord, (2) stop coord, (3) mod probability, (4) haplotype, (5) coverage,
        (6) mod sites, (7) unmod sites, (8) adjusted probability]
    This information is returned and ultimately used to write the output bed file.

    :param ref: Reference name. (str)
    :param pos_start: Start coordinate for region. (int)
    :param pos_stop: Stop coordinate for region. (int)
    :param filtered_basemod_data: List of 2-tuples for each position remaining after filtration. Each 2-tuple is the
    reference position and base mod dat. The list is sorted by reference position (list)
    :param model_dir: Full path to directory containing model. (str)
    :return bed_results: List of sublists with information to write the output bed file. (list)
    """
    logging.debug("coordinates {}: {:,}-{:,}: (4) collect_bed_results_model".format(ref, pos_start, pos_stop))

    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    import tensorflow as tf
    logging.getLogger('tensorflow').setLevel(logging.ERROR)

    # this may or may not do anything to help with the greedy thread situation...
    # tf.config.threading.set_intra_op_parallelism_threads(1)
    # tf.config.threading.set_inter_op_parallelism_threads(1)

    model = tf.keras.models.load_model(model_dir, compile=False)

    total_refpositions, total_normhistos, total_coverages = [], [], []
    hap1_refpositions, hap1_normhistos, hap1_coverages = [], [], []
    hap2_refpositions, hap2_normhistos, hap2_coverages = [], [], []

    # set initial C index for CG location to 0
    previousLocation = 0
    # iterate over reference positions and values (list containing [strand, score, hap]) in filtered_basemod_data
    for (refPosition, modinfoList) in filtered_basemod_data:
        # determine if there is an adjacent prior CG, score appropriately
        if (refPosition - previousLocation) == 2:
            adj = 1
        else:
            adj = 0
        # update CG position
        previousLocation = refPosition

        # build lists for combined haplotypes
        # returns [norm_hist, cov] if min coverage met, otherwise returns empty list
        total_result_list = get_normalized_histo([x[1] for x in modinfoList], adj)
        if total_result_list:
            total_normhistos.append(total_result_list[0])
            total_coverages.append(total_result_list[1])
            total_refpositions.append(refPosition)

        # build lists for hap1
        hap1_result_list = get_normalized_histo([x[1] for x in modinfoList if x[2] == 1], adj)
        if hap1_result_list:
            hap1_normhistos.append(hap1_result_list[0])
            hap1_coverages.append(hap1_result_list[1])
            hap1_refpositions.append(refPosition)

        # build lists for hap2
        hap2_result_list = get_normalized_histo([x[1] for x in modinfoList if x[2] == 2], adj)
        if hap2_result_list:
            hap2_normhistos.append(hap2_result_list[0])
            hap2_coverages.append(hap2_result_list[1])
            hap2_refpositions.append(refPosition)

    # initiate empty list to store all bed results
    bed_results = []
    # run model for total, hap1, hap2, and add to bed results if non-empty list was returned
    apply_model(total_refpositions, total_normhistos, total_coverages, ref, pos_start, pos_stop, model, "Total",
                bed_results)
    apply_model(hap1_refpositions, hap1_normhistos, hap1_coverages, ref, pos_start, pos_stop, model, "hap1",
                bed_results)
    apply_model(hap2_refpositions, hap2_normhistos, hap2_coverages, ref, pos_start, pos_stop, model, "hap2",
                bed_results)

    return bed_results


def run_process_region(arguments):
    """
    Process a given reference region to identify modified bases.
    Uses pickled args (input_file, ref, pos_start, pos_stop) to run
    pileup_from_reads() to get all desired sites (based on modsites option),
    then runs collect_bed_results() to summarize information.

    The sublists will differ between model or count method, but they always share the first 7 elements:
    [(0) ref name, (1) start coord, (2) stop coord, (3) mod probability, (4) haplotype, (5) coverage, ...]

    :param arguments: Pickled list. (list)
    :return bed_results: List of sublists with information to write the output bed file. (list)
    """
    # unpack pickled items:
    # [bam path (str), fasta path (str), modsites option (str),
    #  pileup_mode option (str), model directory path (str),
    #  reference contig name (str), start coordinate (int),
    #  stop coordinate (int), minimum mapping QV (int), haplotype tag name (str)]
    input_bam, input_fasta, modsites, pileup_mode, model_dir, ref, pos_start, pos_stop, min_mapq, hap_tag = arguments
    logging.debug("coordinates {}: {:,}-{:,}: (1) run_process_region: start".format(ref, pos_start, pos_stop))
    # open the input bam file with pysam
    bamIn = pysam.AlignmentFile(input_bam, 'rb')
    # get all ref sites with mods and information from corresponding aligned reads
    basemod_data, cg_sites_read_set = pileup_from_reads(bamIn, ref, pos_start, pos_stop, min_mapq, hap_tag, modsites)
    # filter based on denovo or reference sites
    filtered_basemod_data = filter_basemod_data(basemod_data, cg_sites_read_set, ref, pos_start, pos_stop, input_fasta,
                                                modsites)
    # bam object no longer needed, close file
    bamIn.close()

    if filtered_basemod_data:
        # summarize the mod results, depends on pileup_mode option selected
        if pileup_mode == "count":
            bed_results = collect_bed_results_count(ref, pos_start, pos_stop, filtered_basemod_data)
        elif pileup_mode == "model":
            bed_results = collect_bed_results_model(ref, pos_start, pos_stop, filtered_basemod_data, model_dir)
    else:
        bed_results = []

    logging.debug("coordinates {}: {:,}-{:,}: (5) run_process_region: finish".format(ref, pos_start, pos_stop))

    if len(bed_results) > 1:
        return bed_results
    else:
        return


def run_process_region_wrapper(arguments):
    try:
        return run_process_region(arguments)
    except Exception as e:
        sys.stderr.write("Exception thrown in worker process {}: {}\n".format(os.getpid(), e))
        raise


def run_all_pileup_processing(regions_to_process, threads):
    """
    Function to distribute jobs based on reference regions created.
    Collects results and returns list for writing output bed file.

    The bed results will differ based on model or count method, but they always share the first 7 elements:
    [(0) ref name, (1) start coord, (2) stop coord, (3) mod probability, (4) haplotype, (5) coverage, ...]

    :param regions_to_process: List of sublists defining regions (input_file, ref, pos_start, pos_stop). (list)
    :param threads: Number of threads to use for multiprocessing. (int)
    :return filtered_bed_results: List of sublists with information to write the output bed file. (list)
    """
    logging.info("run_all_pileup_processing: Starting parallel processing.\n")
    # run all jobs

    progress_bar = None
    if sys.stderr.isatty():
        progress_bar = tqdm(total=len(regions_to_process), miniters=1, smoothing=0)

    bed_results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(run_process_region_wrapper, r) for r in regions_to_process]

        # Process results in order of completion
        for future in concurrent.futures.as_completed(futures):
            try:
                bed_result = future.result()
            except:
                # cancel all remaining futures
                for f in futures:
                    f.cancel()
                raise

            bed_results.append(bed_result)
            if progress_bar:
                progress_bar.update(1)

    if progress_bar:
        progress_bar.close()

    logging.info("run_all_pileup_processing: Finished parallel processing.\n")
    # results is a list of sublists, may contain None, remove these
    filtered_bed_results = [i for i in bed_results if i]
    # turn list of lists of sublists into list of sublists
    flattened_bed_results = [i for sublist in filtered_bed_results for i in sublist]

    # ensure bed results are sorted by ref contig name, start position
    logging.info("run_all_pileup_processing: Starting sort for bed results.\n")
    if flattened_bed_results:
        flattened_bed_results.sort(key=itemgetter(0, 1))
        logging.info("run_all_pileup_processing: Finished sort for bed results.\n")

    return flattened_bed_results


def write_output_bed(label, modsites, min_coverage, bed_results):
    """
    Writes output bed file(s) based on information in bed_merge_results (default).
    Separates results into total, hap1, and hap2. If haplotypes not available,
    only total is produced.

    The bed_merge_results list will contain slighty different information depending on the pileup_mode option,
    but the first 7 fields will be identical:

    count-based list
        [(0) ref name, (1) start coord, (2) stop coord, (3) mod probability, (4) haplotype, (5) coverage,
        (6) mod sites, (7) unmod sites, (8) mod score, (9) unmod score]

    OR
    model-based list
        [(0) ref name, (1) start coord, (2) stop coord, (3) mod probability, (4) haplotype, (5) coverage,
        (6) mod sites, (7) unmod sites, (8) adjusted probability]

    :param outname: Name of output bed file to write. (str)
    :param modsites: "reference" or "denovo", for the CpG detection mode. (str)
    :param min_coverage: Minimum coverage to retain a site. (int)
    :param bed_results: List of sublists with information to write the output bed file. (list)
    :return output_files: List of output bed file names that were successfully written. (list)
    """
    logging.info("write_output_bed: Writing unfiltered output bed files.\n")
    out_total = "{}.combined.{}.bed".format(label, modsites)
    out_hap1 = "{}.hap1.{}.bed".format(label, modsites)
    out_hap2 = "{}.hap2.{}.bed".format(label, modsites)
    cov_total = "{}.combined.{}.mincov{}.bed".format(label, modsites, min_coverage)
    cov_hap1 = "{}.hap1.{}.mincov{}.bed".format(label, modsites, min_coverage)
    cov_hap2 = "{}.hap2.{}.mincov{}.bed".format(label, modsites, min_coverage)

    # remove any previous version of output files
    for f in [out_total, out_hap1, out_hap2, cov_total, cov_hap1, cov_hap2]:
        if os.path.exists(f):
            os.remove(f)

    with open(out_total, 'a') as fh_total:
        with open(out_hap1, 'a') as fh_hap1:
            with open(out_hap2, 'a') as fh_hap2:
                for i in bed_results:
                    if i[4] == "Total":
                        fh_total.write("{}\n".format("\t".join([str(j) for j in i])))
                    elif i[4] == "hap1":
                        fh_hap1.write("{}\n".format("\t".join([str(j) for j in i])))
                    elif i[4] == "hap2":
                        fh_hap2.write("{}\n".format("\t".join([str(j) for j in i])))

    # write coverage-filtered versions of bed files
    logging.info(
        "write_output_bed: Writing coverage-filtered output bed files, using min coverage = {}.\n".format(min_coverage))
    output_files = []
    for inBed, covBed in [(out_total, cov_total), (out_hap1, cov_hap1), (out_hap2, cov_hap2)]:
        # if haplotypes not present, the bed files are empty, remove and do not write cov-filtered version
        if os.stat(inBed).st_size == 0:
            os.remove(inBed)
        else:
            output_files.append(inBed)
            # write coverage filtered bed file
            with open(inBed, 'r') as fh_in, open(covBed, 'a') as fh_out:
                for line in fh_in:
                    if int(line.split('\t')[5]) >= min_coverage:
                        fh_out.write(line)
            # check to ensure some sites were written, otherwise remove
            if os.stat(covBed).st_size == 0:
                os.remove(covBed)
            else:
                output_files.append(covBed)

    return output_files


def make_bed_df(bed, pileup_mode):
    """
    Construct a pandas dataframe from a bed file.

    count-based list
        [(0) ref name, (1) start coord, (2) stop coord, (3) % mod sites, (4) haplotype, (5) coverage,
        (6) mod sites, (7) unmod sites, (8) mod score, (9) unmod score]

    OR
    model-based list
        [(0) ref name, (1) start coord, (2) stop coord, (3) mod probability, (4) haplotype, (5) coverage,
        (6) mod sites, (7) unmod sites, (8) adjusted probability]

    :param bed: Name of bed file.
    :param pileup_mode: Site modification calling method. (str: "model", "count")
    :return df: Pandas dataframe.
    """
    logging.debug("make_bed_df: Converting '{}' to pandas dataframe.\n".format(bed))
    if pileup_mode == "count":
        df = pd.read_csv(bed, sep='\t', header=None,
                         names=['chromosome', 'start', 'stop', 'mod_probability', 'haplotype', 'coverage',
                                'modified_bases', 'unmodified_bases', 'mod_score', 'unmod_score'])
        df.drop(columns=['modified_bases', 'unmodified_bases', 'mod_score', 'unmod_score', 'haplotype', 'coverage'],
                inplace=True)

    elif pileup_mode == "model":
        df = pd.read_csv(bed, sep='\t', header=None,
                         names=['chromosome', 'start', 'stop', 'mod_probability', 'haplotype', 'coverage',
                                'modified_bases', 'unmodified_bases', 'adj_prob'])
        df.drop(columns=['haplotype', 'coverage', 'modified_bases', 'unmodified_bases', 'adj_prob'], inplace=True)

    # df.sort_values(by=['chromosome', 'start'], inplace=True)

    return df


def get_bigwig_header_info(input_fasta):
    """
    Get chromosome names and lengths from reference fasta.

    :param input_fasta: Name of reference fasta file.
    :return header: List of tuples, containing [ (ref1, length1), (ref2, length2), ...] .
    """
    logging.debug("get_bigwig_header_info: Getting ref:length info from reference fasta.\n")
    header = []
    with open(input_fasta) as fh:
        for record in SeqIO.parse(fh, "fasta"):
            header.append((record.id, len(record.seq)))
    return header


def write_bigwig_from_df(df, header, outname):
    """
    Function to write a bigwig file using a pandas dataframe from a bed file.

    :param df: Pandas dataframe object (created from bed file).
    :param header: List containing (ref name, length) information. (list of tuples)
    :param outname: Name of bigwig output file to write (OUT.bw).
    """
    logging.debug("write_bigwig_from_df: Writing bigwig file for '{}'.\n".format(outname))
    # first filter reference contigs to match those in bed file
    # get all unique ref contig names from bed
    chroms_present = list(df["chromosome"].unique())
    # header is a list of tuples, filter to keep only those present in bed
    # must also sort reference contigs by name
    filtered_header = sorted([x for x in header if x[0] in chroms_present], key=itemgetter(0))
    for i, j in filtered_header:
        logging.debug("\tHeader includes: '{}', '{}'.".format(i, j))
    # raise error if no reference contig names match
    if not filtered_header:
        logging.error("No reference contig names match between bed file and reference fasta!")
        raise ValueError("No reference contig names match between bed file and reference fasta!")

    # open bigwig object, enable writing mode (default is read only)
    bw = pyBigWig.open(outname, "w")
    # must add header to bigwig prior to writing entries
    bw.addHeader(filtered_header)
    # iterate over ref contig names
    for chrom, length in filtered_header:
        logging.debug("\tAdding entries for '{}'.".format(chrom))
        # subset dataframe by chromosome name
        temp_df = df[df["chromosome"] == chrom]
        logging.debug("\tNumber of entries = {:,}.".format(temp_df.shape[0]))
        # add entries in order specified for bigwig objects:
        # list of chr names: ["chr1", "chr1", "chr1"]
        # list of start coords: [1, 100, 125]
        # list of stop coords: ends=[6, 120, 126]
        # list of vals: values=[0.0, 1.0, 200.0]
        bw.addEntries(list(temp_df["chromosome"]),
                      list(temp_df["start"]),
                      ends=list(temp_df["stop"]),
                      values=list(temp_df["mod_probability"]))
        logging.debug("\tFinished entries for '{}'.\n".format(chrom))
    # close bigwig object
    bw.close()


def convert_bed_to_bigwig(bed_files, fasta, pileup_mode):
    """
    Write bigwig files for each output bed file.

    :param bed_files: List of output bed file names. (list)
    :param fasta: A path to reference fasta file. (str)
    :param pileup_mode: Site modification calling method. (str: "model", "count")
    """
    logging.info("convert_bed_to_bigwig: Converting {} bed files to bigwig files.\n".format(len(bed_files)))
    header = get_bigwig_header_info(fasta)
    for bed in bed_files:
        outname = "{}.bw".format(bed.split(".bed")[0])
        df = make_bed_df(bed, pileup_mode)
        write_bigwig_from_df(df, header, outname)


def main():
    args = get_args()
    setup_logging(args.output_label)
    validate_args(args)
    log_args(args)

    print("\nChunking regions for multiprocessing.")
    regions_to_process = get_regions_to_process(args.bam, args.fasta, args.chunksize, args.modsites,
                                                args.pileup_mode, args.model_dir, args.min_mapq, args.hap_tag)

    print("Running multiprocessing on {:,} chunks.".format(len(regions_to_process)))
    bed_results = run_all_pileup_processing(regions_to_process, args.threads)

    print("Finished multiprocessing.\nWriting bed files.")
    bed_files = write_output_bed(args.output_label, args.modsites, args.min_coverage, bed_results)

    print("Writing bigwig files.")
    convert_bed_to_bigwig(bed_files, args.fasta, args.pileup_mode)

    print("Finished.\n")


if __name__ == '__main__':
    main()

###############################################
#
# All unit tests for this script below, these can be run with:
#     python -m unittest aligned_bam_to_cpg_scores.py
#

import unittest


class Testing(unittest.TestCase):
    def test_mm_ml_parse(self):
        # Create an unmapped record with no MM/ML tags:
        a = pysam.AlignedSegment()
        a.query_name = "qname"
        a.query_sequence = "ACGCCGTATCGTCTCGAGGA"
        a.flag = 4
        a.query_qualities = pysam.qualitystring_to_array("DDDDDEEEEEDDDDDEEEEE")

        meth_tags = get_mm_and_ml_values(a)
        self.assertIsNone(meth_tags)

        # # Add empty MM/ML tags to the alignment record:
        # a.set_tag("MM", "", value_type="Z")
        # a.set_tag("ML", "")
        #
        # meth_tags = get_mm_and_ml_values(a)
        # self.assertIsNone(meth_tags)

        # Add MM/ML tags to the alignment record:
        a.set_tag("MM", "C+m,1,0;", value_type="Z")
        a.set_tag("ML", [100, 150])

        meth_tags = get_mm_and_ml_values(a)
        self.assertIsNotNone(meth_tags)

    def test_parse_mmtag(self):
        # Check standard mmtag parse
        query_seq="ACGCCGTATCGTCTCGAGGA"
        mmtag="C+m,1,0;"
        modcode="C+m"
        base="C"
        reverse=False
        result=parse_mmtag(query_seq, mmtag, modcode, base, reverse)
        self.assertEqual(result, [3,4])

        # Check mmtag parse when skip-type '?' is included
        query_seq="ACGCCGTATCGTCTCGAGGA"
        mmtag="C+m?,1,0;"
        modcode="C+m"
        base="C"
        reverse=False
        result=parse_mmtag(query_seq, mmtag, modcode, base, reverse)
        self.assertEqual(result, [3,4])

        # Check mmtag parse when skip-type '.' is included
        query_seq="ACGCCGTATCGTCTCGAGGA"
        mmtag="C+m.,1,0;"
        modcode="C+m"
        base="C"
        reverse=False
        result=parse_mmtag(query_seq, mmtag, modcode, base, reverse)
        self.assertEqual(result, [3,4])

        # Check mmtag parse when skip-type '.' is included with empty tag
        query_seq="ACGCCGTATCGTCTCGAGGA"
        mmtag="C+m.;"
        modcode="C+m"
        base="C"
        reverse=False
        result=parse_mmtag(query_seq, mmtag, modcode, base, reverse)
        self.assertEqual(result, [])

        # Check mmtag parse when skip-type '.' is included with invalid character
        query_seq="ACGCCGTATCGTCTCGAGGA"
        mmtag="C+m.X"
        modcode="C+m"
        base="C"
        reverse=False
        with self.assertRaises(Exception) as context:
            parse_mmtag(query_seq, mmtag, modcode, base, reverse)

        # Check empty mmtag parse
        query_seq="ACGCCGTATCGTCTCGAGGA"
        mmtag=""
        modcode="C+m"
        base="C"
        reverse=False
        result=parse_mmtag(query_seq, mmtag, modcode, base, reverse)
        self.assertEqual(result, [])

        # Check empty C+m mmtag parse
        query_seq="ACGCCGTATCGTCTCGAGGA"
        mmtag="C+m;"
        modcode="C+m"
        base="C"
        reverse=False
        result=parse_mmtag(query_seq, mmtag, modcode, base, reverse)
        self.assertEqual(result, [])
