'''
UMI-based deduplication of mapped BAM files (`preprocess map --dedup`).

A convenience wrapper around `umi_tools dedup`, which is already part of the project's conda
environment (`requirements.yaml`). Protocols this wrapper does not cover can always be handled
by running `umi_tools` directly against the same BAM directory and rejoining the pipeline at
`analyze build`.
'''
import logging
import re
import os
import subprocess

import pysam

logger = logging.getLogger(__name__)


class MissingUMIError(Exception):
    '''
    Raised when deduplication is requested for a BAM whose read names carry no UMI.

    Deliberately fatal rather than a warning: `umi_tools dedup` given UMI-less reads silently
    falls back to collapsing by alignment position, which for short, deeply-covered tRNA
    transcripts discards genuine molecules and leaves no trace in the output.
    '''


# Separators are tried in this order. `_` is unambiguous in practice; `:` is not, because an
# Illumina read name is already colon-delimited, so a colon can only be recognised as a UMI
# separator by what follows the LAST one.
CANDIDATE_SEPARATORS = ("_", ":")

# A UMI appended to a read name is a short nucleotide string. The bounds are deliberately loose
# -- they exist to reject ordinary trailing name fields (tile/x/y coordinates, which are
# numeric), not to validate a particular protocol's UMI length.
UMI_ALPHABET = frozenset("ACGTN")
UMI_MIN_LENGTH = 4
UMI_MAX_LENGTH = 30

# How many reads to inspect before deciding. A UMI is present on every read or none of them, so
# this only needs to be large enough that a consistent trailing field is not a coincidence.
DETECTION_SAMPLE_SIZE = 100


def _looks_like_umi(candidate):
    if not UMI_MIN_LENGTH <= len(candidate) <= UMI_MAX_LENGTH:
        return False
    return set(candidate.upper()) <= UMI_ALPHABET


def detect_umi_separator(bam_path, sample_size=DETECTION_SAMPLE_SIZE):
    '''
    Determines the character separating the read name from its appended UMI, or None when the
    read names carry no UMI at all.

    Decided from the trailing field of the read name rather than from the mere presence of a
    character, since Illumina read names are themselves colon-delimited. A separator is accepted
    only when every sampled read splits into a plausible UMI of the same length, which keeps an
    ordinary trailing name field (tile coordinates and the like) from being read as a UMI.

    This is a heuristic over read names, not a guarantee: a name whose last field happens to be
    a constant-length nucleotide string would be indistinguishable from an appended UMI. It is
    paired with an explicit refusal when no separator is found, so the failure mode is a stopped
    run rather than a silently mis-parsed one.
    '''
    names = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for record in bam:
            names.append(record.query_name)
            if len(names) >= sample_size:
                break

    if not names:
        return None

    for separator in CANDIDATE_SEPARATORS:
        candidates = []
        for name in names:
            head, found, tail = name.rpartition(separator)
            if not found or not _looks_like_umi(tail):
                candidates = []
                break
            candidates.append(tail)
        if candidates and len({len(curr) for curr in candidates}) == 1:
            return separator

    return None


DEFAULT_DEDUP_METHOD = "directional"

# Suffix for the retained pre-deduplication BAM. Also used as the scratch name during the run,
# since umi_tools cannot write in place: the original is moved aside, umi_tools reads it and
# writes the sample's real name, and the moved-aside copy is then removed unless it was asked
# for. Ordering it that way means an interrupted run leaves the original recoverable rather than
# leaving a truncated file under the name downstream commands read.
PREDEDUP_SUFFIX = ".prededup.bam"


def _prededup_path(bam_path):
    return bam_path[: -len(".bam")] + PREDEDUP_SUFFIX if bam_path.endswith(".bam") else bam_path + PREDEDUP_SUFFIX


def dedup_sample(bam_path, method=DEFAULT_DEDUP_METHOD, keep_prededup=False):
    '''
    Deduplicates one mapped BAM using `umi_tools dedup`, leaving the deduplicated reads under
    the sample's ordinary `<sample>.bam` name so downstream commands need no knowledge of it.

    Returns the path of the retained pre-deduplication BAM, or None when it was discarded.
    '''
    separator = detect_umi_separator(bam_path)
    if separator is None:
        raise MissingUMIError(
            f"No UMI found in the read names of {bam_path}. `--dedup` requires reads whose "
            f"UMI was extracted into the read name (see `preprocess trim -u/--umilength`). "
            f"Running umi_tools on these reads would collapse them by alignment position "
            f"instead, which removes genuine tRNA reads rather than PCR duplicates."
        )

    prededup_path = _prededup_path(bam_path)
    os.replace(bam_path, prededup_path)
    # umi_tools requires a coordinate-sorted, indexed input. `map` already sorts, so only the
    # index is missing.
    pysam.index(prededup_path)

    command = [
        "umi_tools", "dedup",
        f"--stdin={prededup_path}",
        f"--stdout={bam_path}",
        f"--log={bam_path[: -len('.bam')]}_dedup.log" if bam_path.endswith(".bam") else f"--log={bam_path}_dedup.log",
        f"--umi-separator={separator}",
        "--method", method,
    ]
    logger.info("Deduplicating %s with umi_tools (method=%s, separator=%r)",
                os.path.basename(bam_path), method, separator)
    result = subprocess.run(command, capture_output=True)

    if result.returncode != 0:
        # Put the mapping back before surfacing the failure. The original was moved aside rather
        # than copied, so it is the only copy that exists at this point.
        if os.path.exists(bam_path):
            os.remove(bam_path)
        os.replace(prededup_path, bam_path)
        _remove_if_present(prededup_path + ".bai")
        stderr = result.stderr.decode("utf-8", errors="replace") if result.stderr else ""
        raise RuntimeError(
            f"umi_tools dedup failed for {bam_path} (exit {result.returncode}). The original "
            f"BAM has been restored. umi_tools reported:\n{stderr}"
        )

    if keep_prededup:
        return prededup_path

    os.remove(prededup_path)
    _remove_if_present(prededup_path + ".bai")
    return None


def _remove_if_present(path):
    if os.path.exists(path):
        os.remove(path)


def umi_tools_version():
    '''
    Reports the installed umi_tools version, or None when it cannot be determined.

    Recorded alongside each run because deduplication method behavior is version-dependent, and
    a deduplicated BAM is otherwise indistinguishable from a non-deduplicated one.
    '''
    try:
        result = subprocess.run(["umi_tools", "--version"], capture_output=True)
    except FileNotFoundError:
        return None
    if result.returncode != 0:
        return None
    return result.stdout.decode("utf-8", errors="replace").strip().split()[-1]


# What umi_tools reports at the end of a dedup run, and the key each lands under. Anything
# missing is left absent rather than defaulted, so an older umi_tools that does not report a
# field is distinguishable from one that reported zero.
_LOG_FIELDS = (
    ('input_reads',            r'Input Reads:\s*(\d+)',                                int),
    ('output_reads',           r'Number of reads out:\s*(\d+)',                        int),
    ('positions',              r'Total number of positions deduplicated:\s*(\d+)',     int),
    ('mean_umis_per_position', r'Mean number of unique UMIs per position:\s*([\d.]+)', float),
    ('max_umis_per_position',  r'Max\. number of unique UMIs per position:\s*(\d+)',   int),
)

# Column order for the emitted table, sample first.
DEDUP_STATS_COLUMNS = (
    'sample', 'input_reads', 'output_reads', 'retained_pct', 'reads_per_molecule',
    'positions', 'reads_per_position', 'mean_umis_per_position', 'max_umis_per_position',
)


def parse_dedup_log(log_path):
    '''
    Extracts umi_tools' end-of-run statistics from one sample's dedup log.

    Returns only the fields actually found. A missing or unreadable log yields an empty dict --
    the statistics are a diagnostic, and their absence must never fail a run that otherwise
    deduplicated correctly.
    '''
    try:
        with open(log_path) as handle:
            text = handle.read()
    except OSError:
        return {}

    stats = {}
    for key, pattern, cast in _LOG_FIELDS:
        found = re.findall(pattern, text)
        if found:
            stats[key] = cast(found[-1])
    return stats


def _ratio(numerator, denominator, scale=1.0):
    if not numerator or not denominator:
        return None
    return scale * numerator / denominator


def dedup_stats_row(samplename, log_path):
    '''
    One sample's row of the deduplication statistics table.

    The derived rates are what make the raw counts readable: `reads_per_molecule` is the
    duplication level, `reads_per_position` distinguishes a deeply-sequenced library from a
    low-complexity one, and `max_umis_per_position` read against the UMI's own 4^n ceiling says
    whether the tag space is being exhausted at the deepest positions.
    '''
    stats = parse_dedup_log(log_path)
    inp = stats.get('input_reads')
    out = stats.get('output_reads')
    pos = stats.get('positions')
    return {
        'sample': samplename,
        'input_reads': inp,
        'output_reads': out,
        'retained_pct': _ratio(out, inp, 100.0),
        'reads_per_molecule': _ratio(inp, out),
        'positions': pos,
        'reads_per_position': _ratio(inp, pos),
        'mean_umis_per_position': stats.get('mean_umis_per_position'),
        'max_umis_per_position': stats.get('max_umis_per_position'),
    }


def dedup_log_path(bam_path):
    '''Where dedup_sample() writes umi_tools' log for a given bam.'''
    return bam_path[:-len('.bam')] + '_dedup.log' if bam_path.endswith('.bam') else bam_path + '_dedup.log'
