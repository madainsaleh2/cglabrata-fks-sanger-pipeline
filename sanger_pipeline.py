#!/usr/bin/env python3
"""
Sanger Sequencing Analysis Pipeline
Nakaseomyces glabratus (Candida glabrata) — FKS1 / FKS2 Mutation Detection

SCOPE: This pipeline DETECTS AND REPORTS mutations only.  It performs no
clinical interpretation — no resistance calls, drug lists, or MIC levels.
Interpretation is done separately against the lab's own curated database, so
the reference data has a single maintained source of truth.

Workflow per sample × hotspot:
  1. Read AB1 chromatogram files (forward + reverse)
  2. Quality-trim both reads (sliding-window Phred)
  3. Align each read to the reference gene (local pairwise)
  4. Build a reference-anchored consensus (F + R merged)
  5. Report nucleotide mutations (substitutions + indels)
  6. Translate the covered codons; report amino-acid mutations
  7. Flag whether each result is trustworthy (coverage + reliability guards)
  8. Export a CSV report of the observed mutations

Reference genome: Nakaseomyces glabratus ATCC 2001
  GCF_010111755.1 (RefSeq) / GCA_010111755.1 (GenBank)
  Submitted by Johns Hopkins Medical Institution, Feb 3 2020
  FKS1 — Chromosome G (NC_088957.1): CDS 107061-112652
  FKS2 — Chromosome K (NC_088961.1): CDS 375102-380795
"""

# ──────────────────────────────────────────────────────────────────────────────
# Imports
# ──────────────────────────────────────────────────────────────────────────────
import argparse
import csv
import sys
import re
from pathlib import Path
from datetime import datetime
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner


# ──────────────────────────────────────────────────────────────────────────────
# Command-line interface
#
# All paths are supplied at run time — nothing is hardcoded, so the pipeline
# runs unchanged on any machine.  Defaults assume the repository layout:
#     ./Raw   ./references   ./Results
# ──────────────────────────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parent

# Reference FASTA filenames expected inside the --refs directory.
REF_FILENAMES = {
    "FKS1": {"nt": "FKS1_CDS.fasta", "protein": "FKS1_protein.fasta"},
    "FKS2": {"nt": "FKS2_CDS.fasta", "protein": "FKS2_protein.fasta"},
}


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="sanger_pipeline.py",
        description=(
            "Detect FKS1/FKS2 hotspot mutations from paired Sanger .ab1 "
            "chromatograms (Nakaseomyces glabratus). Reports observed mutations "
            "only - no resistance interpretation is applied."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Example:\n  python sanger_pipeline.py --raw example_data --out results\n",
    )
    p.add_argument("--raw", type=Path, default=REPO_ROOT / "Raw",
                   help="Directory containing paired .ab1 chromatogram files")
    p.add_argument("--refs", type=Path, default=REPO_ROOT / "references",
                   help="Directory containing the reference FASTA files")
    p.add_argument("--out", type=Path, default=REPO_ROOT / "Results",
                   help="Directory to write the CSV report into")
    p.add_argument("--prefix", default="FKS_Mutation_Report",
                   help="Filename prefix for the output CSV")

    q = p.add_argument_group("quality thresholds")
    q.add_argument("--trim-q", type=int, default=20,
                   help="Minimum mean Phred score for end trimming")
    q.add_argument("--trim-window", type=int, default=10,
                   help="Sliding-window size (bases) used for end trimming")
    q.add_argument("--min-single-strand-q", type=int, default=30,
                   help="Minimum Phred to call a base covered by one strand only")

    g = p.add_argument_group("hotspot reliability gate")
    g.add_argument("--max-local-conflicts", type=int, default=3,
                   help="Max forward/reverse conflicts inside a hotspot before "
                        "the result is reported as inconclusive")
    g.add_argument("--max-local-indels", type=int, default=0,
                   help="Max indels inside a hotspot before the result is "
                        "reported as inconclusive")
    return p.parse_args(argv)


def resolve_references(ref_dir: Path) -> dict:
    """
    Build the {gene: {'nt': Path, 'protein': Path}} map, failing early with a
    clear message if any reference file is missing.
    """
    refs, missing = {}, []
    for gene, names in REF_FILENAMES.items():
        entry = {}
        for kind, fname in names.items():
            path = ref_dir / fname
            if not path.is_file():
                missing.append(str(path))
            entry[kind] = path
        refs[gene] = entry
    if missing:
        sys.exit("ERROR: missing reference file(s):\n  " + "\n  ".join(missing))
    return refs

# ──────────────────────────────────────────────────────────────────────────────
# Hotspot windows (amino-acid positions, 1-based, inclusive)
# Verified against the ATCC 2001 reference protein in this pipeline.
# Literature: Perlin 2007 · Zimbeck 2010 · Shields 2012 · Arendrup 2012
# ──────────────────────────────────────────────────────────────────────────────
HOTSPOT_AA_RANGES = {
    "FKS1": {
        "HS1": (620, 660),   # key: F625, S629, D632
        "HS2": (1335, 1370), # key: W1354
    },
    "FKS2": {
        "HS1": (654, 690),   # key: F659, S663, D666
        "HS2": (1370, 1400), # key: W1388
    },
}

# ──────────────────────────────────────────────────────────────────────────────
# Quality thresholds
#
# TRIM_Q / TRIM_WINDOW  — used to trim read ends (Q20 = 99% per base).
#     Q20 is adequate for END trimming because most positions are covered by
#     BOTH strands, which cross-validate each other.
#
# MIN_SINGLE_STRAND_Q   — the higher bar (Q30 = 99.9%) required to CALL a
#     mutation at a position covered by only ONE strand.  With no opposite
#     strand to confirm the base, the single read must stand on its own.
#     Positions covered by both strands are exempt (they cross-validate).
# ──────────────────────────────────────────────────────────────────────────────
TRIM_Q               = 20
TRIM_WINDOW          = 10
MIN_SINGLE_STRAND_Q  = 30

# ──────────────────────────────────────────────────────────────────────────────
# Hotspot reliability gate
#
# A mutation call is only trustworthy if the hotspot region itself is clean.
# Real echinocandin-resistance mutations are single-base SUBSTITUTIONS.  A run
# of indels or many FW/RV disagreements INSIDE the hotspot is the signature of
# a low-quality read/misalignment (frameshift-like artifact), which produces
# spurious "mutations" and premature stops.  When the hotspot fails this gate,
# the row is reported as INCONCLUSIVE and the raw calls are preserved separately
# so nothing is silently dropped.
#
# Thresholds chosen from this dataset: the one true positive (S49 S663P) has
# 1 local conflict / 0 local indels; every artifact has ≥9 conflicts and ≥1 indel.
# ──────────────────────────────────────────────────────────────────────────────
HOTSPOT_MAX_LOCAL_CONFLICTS = 3
HOTSPOT_MAX_LOCAL_INDELS    = 0


# ──────────────────────────────────────────────────────────────────────────────
# STEP 1 — Read AB1 chromatogram
# ──────────────────────────────────────────────────────────────────────────────
def read_ab1(filepath: Path):
    """
    Read an AB1 (Applied Biosystems) chromatogram file.
    Returns a SeqRecord with .seq and .letter_annotations['phred_quality'].
    """
    record = SeqIO.read(filepath, "abi")
    return record


# ──────────────────────────────────────────────────────────────────────────────
# STEP 2 — Quality trimming (sliding-window Phred)
# ──────────────────────────────────────────────────────────────────────────────
def trim_by_quality(seq: str, quals: list, threshold: int = TRIM_Q, window: int = TRIM_WINDOW):
    """
    Remove low-quality bases from both ends of a Sanger read.

    Algorithm:
      Scan from the left: advance until a window of `window` bases has
      mean Phred ≥ threshold.  Do the same from the right.  Return the
      trimmed sequence and quality list.

    Parameters
    ----------
    seq       : nucleotide string
    quals     : list of integer Phred scores, one per base
    threshold : minimum acceptable mean quality (default Q20 = 99% accuracy)
    window    : number of bases to average (default 10)
    """
    n = len(quals)
    if n < window:
        return "", []

    # Trim LEFT — find first window whose mean is >= threshold.
    # start stays None if NO window anywhere clears the threshold.
    start = None
    for i in range(n - window + 1):
        if sum(quals[i : i + window]) / window >= threshold:
            start = i
            break

    # If nothing passed, the read has no reliable region → reject it entirely.
    # (Previously this fell through returning the FULL raw read — a silent bug
    #  that fed pure noise into mutation calling for uniformly low-quality reads.)
    if start is None:
        return "", []

    # Trim RIGHT — find last window whose mean is >= threshold
    end = None
    for i in range(n - window, start - 1, -1):
        if sum(quals[i : i + window]) / window >= threshold:
            end = i + window
            break

    if end is None or end <= start:
        return "", []

    return seq[start:end], quals[start:end]


# ──────────────────────────────────────────────────────────────────────────────
# STEP 3 — Reverse-complement with matching quality array
# ──────────────────────────────────────────────────────────────────────────────
def reverse_complement(seq: str, quals: list):
    """
    Reverse-complement a nucleotide sequence and reverse its quality array.

    Why:  The reverse primer reads the antisense strand.  To compare it
          to the forward read (and to the sense-strand reference), we must
          flip the sequence direction AND swap each base to its complement
          (A↔T, G↔C).  Quality scores simply reverse — they are per-base
          and follow the same direction flip.
    """
    complement = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    rc_seq   = seq[::-1].translate(complement)
    rc_quals = quals[::-1]
    return rc_seq, rc_quals


# ──────────────────────────────────────────────────────────────────────────────
# STEP 4a — Align one read to the reference (local pairwise)
# ──────────────────────────────────────────────────────────────────────────────
def make_aligner() -> PairwiseAligner:
    """
    Configure a local pairwise aligner suitable for Sanger amplicons.
    Scoring matches standard practice for closely related sequences.
    """
    aln = PairwiseAligner()
    aln.mode            = "local"   # find best sub-region (handles unsequenced ends)
    aln.match_score     =  2        # reward for identical bases
    aln.mismatch_score  = -1        # penalty for SNPs
    aln.open_gap_score  = -3        # penalty to START a gap (indels are rare in Sanger)
    aln.extend_gap_score= -0.5      # per-base extension (prefer one long gap over many)
    return aln


def align_to_reference(query_seq: str, ref_seq: str, aligner: PairwiseAligner):
    """
    Align query to reference.  Returns the best Alignment object.
    Returns None if no alignment reaches a minimum score threshold.
    """
    # An empty query (read rejected during trimming) cannot be aligned.
    if not query_seq:
        return None
    alignments = aligner.align(ref_seq, query_seq)
    best = next(iter(alignments), None)
    if best is None or best.score < 100:   # sanity threshold
        return None
    return best


# ──────────────────────────────────────────────────────────────────────────────
# STEP 4b — Extract a reference-position → query-position map from alignment
# ──────────────────────────────────────────────────────────────────────────────
def alignment_map(alignment) -> tuple[dict, list]:
    """
    Parse a PairwiseAlignment object.

    Returns
    -------
    ref_to_query : dict mapping ref_position (0-based) → query_position (0-based)
                   Only positions that are MATCHED (no gap on either side) are included.
    indels       : list of human-readable indel descriptions (deletions / insertions)

    How alignment.aligned works:
      alignment.aligned is a 2-tuple of arrays.
        [0] = list of (start, end) blocks on the REFERENCE (gaps excluded)
        [1] = list of (start, end) blocks on the QUERY    (gaps excluded)
      Matching blocks are paired element-by-element.
      Gaps are inferred from the spaces *between* consecutive blocks.
    """
    ref_blocks = alignment.aligned[0]
    qry_blocks = alignment.aligned[1]

    ref_to_query = {}
    indels       = []

    for i, ((r_start, r_end), (q_start, q_end)) in enumerate(
        zip(ref_blocks, qry_blocks)
    ):
        # Map every matched position
        for r_pos, q_pos in zip(range(r_start, r_end), range(q_start, q_end)):
            ref_to_query[r_pos] = q_pos

        # Detect gap between this block and the next
        if i + 1 < len(ref_blocks):
            next_r_start = ref_blocks[i + 1][0]
            next_q_start = qry_blocks[i + 1][0]

            ref_gap = next_r_start - r_end   # > 0 → deletion in sample
            qry_gap = next_q_start - q_end   # > 0 → insertion in sample

            if ref_gap > 0 and qry_gap == 0:
                nt_start = r_end + 1       # 1-based
                nt_end   = next_r_start    # 1-based
                indels.append(
                    f"del_nt{nt_start}_{nt_end} (size {ref_gap} bp)"
                )
            elif qry_gap > 0 and ref_gap == 0:
                nt_pos = r_end + 1         # position after last matched ref base
                indels.append(
                    f"ins_nt{nt_pos} (size {qry_gap} bp)"
                )

    return ref_to_query, indels


# ──────────────────────────────────────────────────────────────────────────────
# STEP 4c — Merge forward and reverse maps into one consensus
# ──────────────────────────────────────────────────────────────────────────────
def build_consensus(
    fw_seq: str, fw_quals: list, fw_map: dict,
    rv_seq: str, rv_quals: list, rv_map: dict,
    min_single_strand_q: int = MIN_SINGLE_STRAND_Q,
) -> tuple[dict, list]:
    """
    Merge the forward and reverse read alignments into a single consensus.

    Each entry in the returned dictionary:
        ref_position → {
            "base"          : nucleotide call (str),
            "quality"       : Phred score (int),
            "source"        : "FW" | "RV" | "both_agree" | "both_conflict",
            "single_strand" : True if covered by only one read,
            "callable"      : True if trustworthy enough to call a mutation on,
        }

    Rules:
      - Position covered by one read only → use that base.
      - Covered by both, same base → use higher-quality copy.
      - Covered by both, different bases → flag as CONFLICT; keep higher-quality base.

    Callability (single-strand-aware confidence):
      - Both-strand positions cross-validate → always callable.
      - Single-strand positions must reach `min_single_strand_q` (Q30) to be
        callable; otherwise they are MASKED (a mutation there is not reported,
        because one low-quality read is not enough to trust a base call).
    """
    consensus = {}
    conflicts = []

    # Add forward read positions
    for r_pos, q_pos in fw_map.items():
        consensus[r_pos] = {
            "base":    fw_seq[q_pos],
            "quality": fw_quals[q_pos],
            "source":  "FW",
        }

    # Merge reverse read positions
    for r_pos, q_pos in rv_map.items():
        rv_base = rv_seq[q_pos]
        rv_qual = rv_quals[q_pos]

        if r_pos not in consensus:
            consensus[r_pos] = {
                "base":    rv_base,
                "quality": rv_qual,
                "source":  "RV",
            }
        else:
            existing = consensus[r_pos]
            if existing["base"].upper() == rv_base.upper():
                # Agreement — keep whichever has the higher quality
                if rv_qual > existing["quality"]:
                    consensus[r_pos] = {
                        "base":    rv_base,
                        "quality": rv_qual,
                        "source":  "both_agree",
                    }
                else:
                    consensus[r_pos]["source"] = "both_agree"
            else:
                # Conflict — use the better-quality base and flag
                if rv_qual > existing["quality"]:
                    winner_base, winner_qual = rv_base, rv_qual
                else:
                    winner_base, winner_qual = existing["base"], existing["quality"]
                consensus[r_pos] = {
                    "base":    winner_base,
                    "quality": winner_qual,
                    "source":  "both_conflict",
                }
                conflicts.append(
                    f"nt{r_pos + 1}: FW={existing['base']}(Q{existing['quality']}) "
                    f"RV={rv_base}(Q{rv_qual}) → kept {winner_base}"
                )

    # Tag every position with single-strand status and callability.
    for entry in consensus.values():
        entry["single_strand"] = entry["source"] in ("FW", "RV")
        if entry["single_strand"]:
            entry["callable"] = entry["quality"] >= min_single_strand_q
        else:
            entry["callable"] = True   # both strands agree/resolved → trusted

    return consensus, conflicts


# ──────────────────────────────────────────────────────────────────────────────
# STEP 5 — Find nucleotide mutations
# ──────────────────────────────────────────────────────────────────────────────
def find_nt_mutations(consensus: dict, ref_seq: str) -> tuple[list[str], list[str]]:
    """
    Compare the consensus to the reference at every aligned position.

    Returns
    -------
    mutations       : substitutions in HGVS-like notation, e.g. "T1875C"
                      (position 1875 of the CDS changed from T to C).
    masked_variants : differences that WERE seen but suppressed because the
                      position is single-strand and below Q30 — reported
                      separately so the low-confidence signal is not lost.

    Positions with 'N' are skipped (ambiguous call).  Non-callable single-strand
    positions are masked.  Indels are handled separately by alignment_map().
    """
    mutations       = []
    masked_variants = []
    for ref_pos in sorted(consensus):
        entry       = consensus[ref_pos]
        sample_base = entry["base"].upper()
        ref_base    = ref_seq[ref_pos].upper()
        if sample_base == "N":
            continue
        if sample_base != ref_base:
            nt_pos = ref_pos + 1   # convert to 1-based
            if entry["callable"]:
                mutations.append(f"{ref_base}{nt_pos}{sample_base}")
            else:
                masked_variants.append(
                    f"{ref_base}{nt_pos}{sample_base} "
                    f"(single-strand Q{entry['quality']}, masked)"
                )
    return mutations, masked_variants


# ──────────────────────────────────────────────────────────────────────────────
# STEP 6 — Translate and find amino-acid mutations
# ──────────────────────────────────────────────────────────────────────────────
def find_aa_mutations(
    consensus: dict, ref_seq: str, ref_protein: str
) -> tuple[list[str], list[str], list[str]]:
    """
    Reconstruct each codon covered by the consensus (all 3 positions present),
    translate it, and compare to the reference protein.

    Returns
    -------
    aa_mutations : list of strings e.g. ["F625S", "S629L"]
    synonymous   : list of silent NT changes e.g. ["T1875C (silent)"]
    masked_aa    : AA changes seen but suppressed because ≥1 codon base is a
                   single-strand call below Q30 — reported separately.

    How codon reconstruction works:
      The reference FASTA starts at the ATG of the CDS (position 0).
      Codon N (1-based) occupies reference positions (N-1)*3, (N-1)*3+1, (N-1)*3+2.
      We collect every codon where all 3 nt positions are in the consensus,
      substitute any sample bases, translate, and compare to the reference AA.

    Confidence rule:
      A codon is only CALLED if all three of its nucleotide positions are
      `callable` (both-strand, or single-strand ≥Q30).  If any base is a
      low-confidence single-strand call, the codon is masked — one shaky base
      changes the amino acid, so the whole codon is untrustworthy.
    """
    if not consensus:
        return [], [], []

    positions   = sorted(consensus.keys())
    start_ref   = positions[0]
    end_ref     = positions[-1]

    # First and last complete codon indices (0-based codon index)
    first_codon = (start_ref + 2) // 3           # ceil(start_ref / 3)
    last_codon  = end_ref // 3                    # floor(end_ref  / 3)

    aa_mutations = []
    synonymous   = []
    masked_aa    = []

    for codon_idx in range(first_codon, last_codon + 1):
        nt0 = codon_idx * 3
        nt1 = nt0 + 1
        nt2 = nt0 + 2

        # Need all 3 positions covered by the consensus
        if nt0 not in consensus or nt1 not in consensus or nt2 not in consensus:
            continue

        sample_codon = (
            consensus[nt0]["base"].upper()
            + consensus[nt1]["base"].upper()
            + consensus[nt2]["base"].upper()
        )
        ref_codon = ref_seq[nt0:nt2 + 1].upper()

        if "N" in sample_codon:
            continue

        sample_aa = str(Seq(sample_codon).translate(table=1))
        aa_pos    = codon_idx + 1   # 1-based AA position

        if aa_pos > len(ref_protein):
            continue

        ref_aa = ref_protein[aa_pos - 1]

        # Is every base of this codon trustworthy?
        codon_callable = (
            consensus[nt0]["callable"]
            and consensus[nt1]["callable"]
            and consensus[nt2]["callable"]
        )

        if sample_aa != ref_aa:
            label = (f"{ref_aa}{aa_pos}* (premature stop)"
                     if sample_aa == "*" else f"{ref_aa}{aa_pos}{sample_aa}")
            if codon_callable:
                aa_mutations.append(label)
            else:
                # Lowest quality among the single-strand bases in this codon
                low_q = min(
                    consensus[p]["quality"]
                    for p in (nt0, nt1, nt2)
                    if not consensus[p]["callable"]
                )
                masked_aa.append(f"{label} (single-strand Q{low_q}, masked)")
        elif sample_codon != ref_codon and codon_callable:
            # Same amino acid but different codon → synonymous change
            synonymous.append(f"{ref_codon}→{sample_codon} at nt{nt0+1} (synonymous, {ref_aa}{aa_pos})")

    return aa_mutations, synonymous, masked_aa


# ──────────────────────────────────────────────────────────────────────────────
# STEP 7 — Split AA mutations into hotspot vs. non-hotspot, then classify
# ──────────────────────────────────────────────────────────────────────────────
def split_mutations_by_window(
    aa_mutations: list[str],
    hs_aa_start: int,
    hs_aa_end: int,
) -> tuple[list[str], list[str]]:
    """
    Partition aa_mutations into two lists:
      in_hotspot  — AA position falls within [hs_aa_start, hs_aa_end]
      out_hotspot — everything else

    Only in_hotspot mutations drive the resistance classification.
    Out-of-hotspot changes are natural polymorphisms or amplicon noise
    and should not inflate the clinical call.
    """
    in_hotspot  = []
    out_hotspot = []
    for mut in aa_mutations:
        core = mut.split()[0]
        m = re.match(r"[A-Z*](\d+)[A-Z*\*]", core)
        if m:
            pos = int(m.group(1))
            if hs_aa_start <= pos <= hs_aa_end:
                in_hotspot.append(mut)
            else:
                out_hotspot.append(mut)
        else:
            out_hotspot.append(mut)   # can't parse position → conservatively exclude
    return in_hotspot, out_hotspot


def hotspot_coverage(
    consensus: dict, hs_aa_start: int, hs_aa_end: int
) -> tuple[int, int]:
    """
    Count how many amino-acid positions in the hotspot window have a fully
    CALLABLE codon (all 3 nt positions present and callable).

    Returns (covered_codons, total_codons).

    Why this matters:  a hotspot with no callable coverage will produce zero
    mutations and look mutation-free — but that is a FALSE reassurance if the
    region simply was not sequenced (or was entirely masked for low quality).
    We use this to downgrade such calls to "inconclusive".
    """
    total   = hs_aa_end - hs_aa_start + 1
    covered = 0
    for aa_pos in range(hs_aa_start, hs_aa_end + 1):
        codon_idx = aa_pos - 1          # 0-based
        nt0 = codon_idx * 3
        if all(
            (nt0 + k) in consensus and consensus[nt0 + k]["callable"]
            for k in (0, 1, 2)
        ):
            covered += 1
    return covered, total


def count_in_window(entries: list[str], hs_aa_start: int, hs_aa_end: int) -> int:
    """
    Count how many entries (conflict/indel strings that contain an 'nt<pos>'
    token) fall inside the hotspot's nucleotide window.  The nt window is
    derived from the AA window: nt = (aa-1)*3+1 .. aa*3 (1-based).
    """
    nt_start = (hs_aa_start - 1) * 3 + 1
    nt_end   = hs_aa_end * 3
    n = 0
    for e in entries:
        m = re.search(r"nt(\d+)", e)
        if m and nt_start <= int(m.group(1)) <= nt_end:
            n += 1
    return n


def call_status(hs_aa_mutations: list[str]) -> str:
    """
    Report WHAT WAS FOUND in the hotspot window — no clinical interpretation.

    This pipeline deliberately does NOT judge drug resistance.  It reports the
    observed mutations objectively; interpretation against a curated resistance
    database is done separately (outside this script) so the reference data has
    a single, maintained source of truth.

    Returns
    -------
    "Mutation detected"  — one or more amino-acid changes in the hotspot window
    "No mutation"        — hotspot matches the wild-type reference
    """
    return "Mutation detected" if hs_aa_mutations else "No mutation"


# ──────────────────────────────────────────────────────────────────────────────
# STEP 8 — Compute coverage and range statistics
# ──────────────────────────────────────────────────────────────────────────────
def coverage_stats(consensus: dict, ref_protein: str) -> dict:
    """
    Derive nucleotide range, amino-acid range, and coverage metrics
    from the merged consensus dictionary.
    """
    if not consensus:
        return {
            "nt_start": "N/A", "nt_end": "N/A", "nt_range": "N/A",
            "aa_start": "N/A", "aa_end": "N/A", "aa_range": "N/A",
            "consensus_length": 0,
        }

    positions = sorted(consensus.keys())
    nt_start  = positions[0]  + 1   # 1-based
    nt_end    = positions[-1] + 1   # 1-based

    # AA range: codon containing first nt to codon containing last nt
    aa_start  = positions[0]  // 3 + 1   # 1-based
    aa_end    = positions[-1] // 3 + 1   # 1-based

    return {
        "nt_start":         nt_start,
        "nt_end":           nt_end,
        "nt_range":         f"CDS nt {nt_start}–{nt_end}",
        "aa_start":         aa_start,
        "aa_end":           aa_end,
        "aa_range":         f"AA {aa_start}–{aa_end}",
        "consensus_length": len(positions),
    }


# ──────────────────────────────────────────────────────────────────────────────
# STEP 9 — Per-sample, per-hotspot analysis
# ──────────────────────────────────────────────────────────────────────────────
def analyse_hotspot(
    sample_id: str,
    gene: str,
    hotspot: str,               # canonical hotspot for biology (HS1 / HS2)
    fw_path: Path | None,
    rv_path: Path | None,
    ref_nt_seq: str,
    ref_protein: str,
    aligner: PairwiseAligner,
    hotspot_label: str = None,  # raw label for reporting (e.g. HS1b); defaults to hotspot
    trim_threshold: int = TRIM_Q,
    trim_window: int    = TRIM_WINDOW,
    min_single_strand_q: int = MIN_SINGLE_STRAND_Q,
    max_local_conflicts: int = HOTSPOT_MAX_LOCAL_CONFLICTS,
    max_local_indels: int    = HOTSPOT_MAX_LOCAL_INDELS,
) -> dict:
    """
    Full analysis for one sample × gene × hotspot combination.
    Returns a flat dictionary of all result fields (one CSV row).

    `hotspot` drives the biology (AA window); `hotspot_label` is what gets
    printed/exported, preserving primer-batch suffixes like HS1b.
    Either read path may be None (single-strand run); a missing read is treated
    as an empty sequence and the single-strand-aware calling logic applies.
    Quality thresholds are passed in so they can be set from the command line.
    """
    flags = []

    result = {
        "Sample_ID": sample_id,
        "Gene":      gene,
        "Hotspot":   hotspot_label or hotspot,
    }

    # ── Read AB1 files (either may be absent) ─────────────────────────────────
    def _load(path):
        if path is None:
            return "", []
        rec = read_ab1(path)
        return str(rec.seq), rec.letter_annotations["phred_quality"]

    try:
        fw_seq_raw, fw_quals_raw = _load(fw_path)
        rv_seq_raw, rv_quals_raw = _load(rv_path)
    except Exception as e:
        # Fill every data column with a marker so the CSV row stays well-formed.
        for col in CSV_COLUMNS:
            result.setdefault(col, "FILE_ERROR")
        result["QC_Flags"] = f"FILE_ERROR: {e}"
        return result

    if fw_path is None:
        flags.append("No forward read provided")
    if rv_path is None:
        flags.append("No reverse read provided")

    result["FW_Length_Raw"] = len(fw_seq_raw)
    result["RV_Length_Raw"] = len(rv_seq_raw)

    # ── Trim both reads ───────────────────────────────────────────────────────
    fw_seq, fw_quals = trim_by_quality(fw_seq_raw, fw_quals_raw, trim_threshold, trim_window)
    rv_seq, rv_quals = trim_by_quality(rv_seq_raw, rv_quals_raw, trim_threshold, trim_window)

    result["FW_Length_Trimmed"] = len(fw_seq)
    result["RV_Length_Trimmed"] = len(rv_seq)
    result["FW_Avg_Quality"]    = f"{sum(fw_quals)/len(fw_quals):.1f}" if fw_quals else "N/A"
    result["RV_Avg_Quality"]    = f"{sum(rv_quals)/len(rv_quals):.1f}" if rv_quals else "N/A"

    if fw_path is not None and len(fw_seq) < 100:
        flags.append("FW read too short after trimming")
    if rv_path is not None and len(rv_seq) < 100:
        flags.append("RV read too short after trimming")

    # ── Reverse-complement the reverse read ───────────────────────────────────
    rv_rc_seq, rv_rc_quals = reverse_complement(rv_seq, rv_quals)

    # ── Align both reads to reference ─────────────────────────────────────────
    fw_aln = align_to_reference(fw_seq, ref_nt_seq, aligner)
    rv_aln = align_to_reference(rv_rc_seq, ref_nt_seq, aligner)

    if fw_aln is None:
        flags.append("FW read did not align to reference")
    if rv_aln is None:
        flags.append("RV read did not align to reference")

    # ── Build consensus ───────────────────────────────────────────────────────
    fw_map, fw_indels = (alignment_map(fw_aln) if fw_aln else ({}, []))
    rv_map, rv_indels = (alignment_map(rv_aln) if rv_aln else ({}, []))

    all_indels = fw_indels + rv_indels

    consensus, conflicts = build_consensus(
        fw_seq, fw_quals, fw_map,
        rv_rc_seq, rv_rc_quals, rv_map,
        min_single_strand_q=min_single_strand_q,
    )

    if conflicts:
        flags.append(f"{len(conflicts)} base conflict(s) between FW and RV")

    cov = coverage_stats(consensus, ref_protein)
    result.update({
        "Consensus_Length_bp": cov["consensus_length"],
        "NT_Range":            cov["nt_range"],
        "AA_Range":            cov["aa_range"],
    })

    # ── Strand coverage summary ───────────────────────────────────────────────
    n_both   = sum(1 for e in consensus.values() if not e["single_strand"])
    n_single = sum(1 for e in consensus.values() if e["single_strand"])
    n_masked = sum(1 for e in consensus.values()
                   if e["single_strand"] and not e["callable"])
    total    = len(consensus)

    if total == 0:
        strand_cov = "No coverage"
    elif n_single == 0:
        strand_cov = "Both strands (100%)"
    elif n_both == 0:
        only = "FW" if fw_map else ("RV" if rv_map else "?")
        strand_cov = f"Single strand only ({only})"
    else:
        strand_cov = (f"Mixed ({100*n_both//total}% both, "
                      f"{100*n_single//total}% single)")

    result["Strand_Coverage"]            = strand_cov
    result["Single_Strand_Positions"]    = n_single
    result["LowConf_Masked_Positions"]   = n_masked

    # ── Nucleotide mutations ───────────────────────────────────────────────────
    nt_mutations, nt_masked = find_nt_mutations(consensus, ref_nt_seq)
    result["NT_Mutations"]        = "; ".join(nt_mutations) if nt_mutations else "None"
    result["NT_Masked_LowConf"]   = "; ".join(nt_masked)    if nt_masked    else "None"
    result["Indels"]              = "; ".join(all_indels)   if all_indels   else "None"

    # ── Amino-acid mutations ───────────────────────────────────────────────────
    aa_mutations, synonymous, aa_masked = find_aa_mutations(
        consensus, ref_nt_seq, ref_protein
    )

    # Split into mutations inside vs. outside the hotspot window
    hs_range = HOTSPOT_AA_RANGES[gene][hotspot]
    hs_aa_mut, other_aa_mut = split_mutations_by_window(
        aa_mutations, hs_range[0], hs_range[1]
    )
    # Masked AA changes that fall inside the hotspot window are the dangerous
    # ones: a real resistance mutation could be hiding behind a low-quality base.
    hs_masked, _ = split_mutations_by_window(aa_masked, hs_range[0], hs_range[1])

    result["AA_Mutations_in_Hotspot"]      = "; ".join(hs_aa_mut)    if hs_aa_mut    else "None"
    result["AA_Mutations_outside_Hotspot"] = "; ".join(other_aa_mut) if other_aa_mut else "None"
    result["AA_Masked_LowConf"]            = "; ".join(aa_masked)    if aa_masked    else "None"
    result["Synonymous_NT_Changes"]        = "; ".join(synonymous)   if synonymous   else "None"

    # ── Reliability check ─────────────────────────────────────────────────────
    n_conflicts = len(conflicts)
    if n_conflicts > 50:
        flags.append(f"UNRELIABLE: {n_conflicts} FW/RV conflicts — repeat sequencing advised")
    if not fw_map and not rv_map:
        flags.append("UNRELIABLE: neither read aligned — no data for this hotspot")
    if n_masked:
        flags.append(f"{n_masked} single-strand position(s) masked (Q<{min_single_strand_q})")
    if hs_masked:
        flags.append(
            "REVIEW: possible mutation in hotspot masked by low single-strand "
            f"quality → {'; '.join(hs_masked)} — repeat sequencing advised"
        )

    # ── Hotspot coverage guard ────────────────────────────────────────────────
    # A "No mutation" call is only meaningful if the hotspot was actually read.
    hs_covered, hs_total = hotspot_coverage(consensus, hs_range[0], hs_range[1])
    hs_cov_pct = (100 * hs_covered // hs_total) if hs_total else 0
    result["Hotspot_Coverage"] = f"{hs_covered}/{hs_total} codons ({hs_cov_pct}%)"

    # ── Hotspot reliability gate ──────────────────────────────────────────────
    # Count discordance and indels INSIDE the hotspot window — the artifact
    # signature that produces spurious mutations / premature stops.
    hs_local_conflicts = count_in_window(conflicts,  hs_range[0], hs_range[1])
    hs_local_indels    = count_in_window(all_indels, hs_range[0], hs_range[1])
    result["Hotspot_Local_Conflicts"] = hs_local_conflicts
    result["Hotspot_Local_Indels"]    = hs_local_indels

    hotspot_reliable = (
        hs_local_indels    <= max_local_indels and
        hs_local_conflicts <= max_local_conflicts
    )

    # ── Report what was found (no clinical interpretation) ────────────────────
    # The data-quality guards below still apply: they decide whether the result
    # is TRUSTWORTHY, not whether a mutation confers resistance.
    status = call_status(hs_aa_mut)

    if not hotspot_reliable:
        # The hotspot is too noisy to trust any call from it — substitutions and
        # premature stops here are almost certainly artifacts.  Preserve the raw
        # calls in AA_Mutations_in_Hotspot but do not report them as findings.
        status = "INCONCLUSIVE — unreliable hotspot (indels/discordance)"
        flags.append(
            f"INCONCLUSIVE: hotspot has {hs_local_indels} indel(s) and "
            f"{hs_local_conflicts} FW/RV conflict(s) — calls not trustworthy; "
            "repeat sequencing"
        )
    # Downgrade a "No mutation" verdict when the hotspot is not fully covered —
    # absence of a mutation over unsequenced codons is not evidence of absence.
    elif status == "No mutation" and hs_cov_pct < 100:
        if hs_covered == 0:
            status = "INCONCLUSIVE — hotspot not covered"
            flags.append(
                f"INCONCLUSIVE: 0/{hs_total} hotspot codons callable "
                "— region not sequenced or fully masked; repeat sequencing"
            )
        else:
            status = f"No mutation (partial: {hs_cov_pct}% of hotspot covered)"
            flags.append(
                f"PARTIAL hotspot coverage ({hs_covered}/{hs_total} codons) "
                "— uncovered codons cannot be assessed"
            )

    result["Call_Status"] = status

    # ── Conflicts detail ───────────────────────────────────────────────────────
    result["Conflicts"] = "; ".join(conflicts) if conflicts else "None"

    # ── QC flags ──────────────────────────────────────────────────────────────
    result["QC_Flags"] = "; ".join(flags) if flags else "Pass"

    return result


# ──────────────────────────────────────────────────────────────────────────────
# File discovery — parse sample IDs and hotspot labels from filenames
# ──────────────────────────────────────────────────────────────────────────────
# Two accepted filename orderings (case-insensitive), with '-' or '_' before
# the direction token, an optional letter suffix on the hotspot (e.g. HS1b),
# and a forgiving direction token (FW/RV, plus F/R/W typo variants):
#   1. sample-first : S2-FKS1-HS1b-FW.ab1   /  S49-FKS1-HS1_FW.ab1
#   2. sample-last  : FKS1-HS2b-S2_FW.ab1
_PAT_SAMPLE_FIRST = re.compile(
    r"^(?P<sample>S\d+)-(?P<gene>FKS[12])-(?P<hotspot>HS[12][A-Za-z]?)[-_](?P<dir>FW|RV|F|R|W)\.ab1$",
    re.IGNORECASE,
)
_PAT_SAMPLE_LAST = re.compile(
    r"^(?P<gene>FKS[12])-(?P<hotspot>HS[12][A-Za-z]?)-(?P<sample>S\d+)[-_](?P<dir>FW|RV|F|R|W)\.ab1$",
    re.IGNORECASE,
)


def normalise_direction(token: str) -> tuple[str, bool]:
    """
    Map a raw direction token to canonical 'FW' or 'RV'.
    Returns (direction, was_corrected).  'W' and 'F' → FW; 'R' → RV.
    """
    t = token.upper()
    if t in ("FW", "F", "W"):
        return "FW", (t != "FW")
    if t in ("RV", "R"):
        return "RV", (t != "RV")
    return t, True


def canonical_hotspot(raw_hotspot: str) -> str:
    """
    Strip a trailing primer/batch letter so the biological hotspot is recovered.
    'HS1b' → 'HS1', 'HS2B' → 'HS2'.  Used for AA-range and resistance-DB lookup.
    """
    return raw_hotspot[:3].upper()


def discover_samples(raw_dir: Path) -> tuple[dict, list]:
    """
    Scan the Raw directory for AB1 files under either naming convention.

    Returns
    -------
    samples : nested dict
        { sample_id: { (gene, raw_hotspot): {"FW": Path, "RV": Path} } }
        raw_hotspot preserves the 'b' suffix (e.g. 'HS1b') for provenance;
        the canonical hotspot for biology is derived via canonical_hotspot().
    notes   : list of human-readable parsing notes (typo corrections, skips)
    """
    samples = defaultdict(lambda: defaultdict(dict))
    notes   = []

    for f in sorted(raw_dir.glob("*.ab1")):
        m = _PAT_SAMPLE_FIRST.match(f.name) or _PAT_SAMPLE_LAST.match(f.name)
        if not m:
            notes.append(f"[WARN] Unrecognised filename, skipped: {f.name}")
            continue

        sid  = m.group("sample").upper()
        gene = m.group("gene").upper()
        hs   = m.group("hotspot").upper()          # keep 'b' → e.g. HS1B
        d, corrected = normalise_direction(m.group("dir"))

        if corrected:
            notes.append(
                f"[FIX] {f.name}: direction '{m.group('dir')}' → {d}"
            )

        # Warn if two files map to the same (sample, hotspot, direction) slot.
        if d in samples[sid][(gene, hs)]:
            notes.append(
                f"[WARN] Duplicate {sid} {gene} {hs} {d}: "
                f"{samples[sid][(gene, hs)][d].name} vs {f.name} — keeping latter"
            )
        samples[sid][(gene, hs)][d] = f

    return samples, notes


# ──────────────────────────────────────────────────────────────────────────────
# CSV report generation
# ──────────────────────────────────────────────────────────────────────────────
CSV_COLUMNS = [
    # ── Sample identity ──────────────────────────────────────────────────────
    "Sample_ID",
    "Gene",
    "Hotspot",
    "Hotspot_AA_Range",
    # ── Read quality ─────────────────────────────────────────────────────────
    "FW_Length_Raw",
    "RV_Length_Raw",
    "FW_Length_Trimmed",
    "RV_Length_Trimmed",
    "FW_Avg_Quality",
    "RV_Avg_Quality",
    # ── Coverage ─────────────────────────────────────────────────────────────
    "Consensus_Length_bp",
    "NT_Range",
    "AA_Range",
    "Strand_Coverage",
    "Hotspot_Coverage",
    "Hotspot_Local_Conflicts",
    "Hotspot_Local_Indels",
    "Single_Strand_Positions",
    "LowConf_Masked_Positions",
    # ── Nucleotide-level findings ─────────────────────────────────────────────
    "NT_Mutations",
    "NT_Masked_LowConf",
    "Synonymous_NT_Changes",
    "Indels",
    "Conflicts",
    # ── Amino-acid findings (split by hotspot window) ─────────────────────────
    "AA_Mutations_in_Hotspot",
    "AA_Mutations_outside_Hotspot",
    "AA_Masked_LowConf",
    # ── Result status (what was found — NOT a clinical interpretation) ────────
    "Call_Status",
    # ── QC and provenance ─────────────────────────────────────────────────────
    "QC_Flags",
    "Analysis_Date",
    "Reference_Genome",
]

def write_csv(rows: list[dict], out_path: Path):
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n  Report saved → {out_path}")


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────
def main(argv=None):
    args = parse_args(argv)

    print("=" * 70)
    print("  Sanger Sequencing Pipeline — N. glabratus FKS mutation detection")
    print(f"  Run date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)
    print(f"  Raw reads   : {args.raw}")
    print(f"  References  : {args.refs}")
    print(f"  Output      : {args.out}")
    print(f"  Thresholds  : trim Q>={args.trim_q}/win{args.trim_window}, "
          f"single-strand Q>={args.min_single_strand_q}, "
          f"gate <={args.max_local_conflicts} conflicts / <={args.max_local_indels} indels")

    if not args.raw.is_dir():
        sys.exit(f"ERROR: --raw directory not found: {args.raw}")
    if not args.refs.is_dir():
        sys.exit(f"ERROR: --refs directory not found: {args.refs}")
    args.out.mkdir(parents=True, exist_ok=True)

    refs = resolve_references(args.refs)

    # Load reference sequences
    print("\n[1] Loading reference sequences …")
    ref_data = {}
    for gene, paths in refs.items():
        nt_rec  = SeqIO.read(paths["nt"],      "fasta")
        pr_rec  = SeqIO.read(paths["protein"], "fasta")
        ref_data[gene] = {
            "nt_seq":  str(nt_rec.seq).upper(),
            "protein": str(pr_rec.seq),
            "nt_id":   nt_rec.id,
        }
        print(f"   {gene}  NT  {len(ref_data[gene]['nt_seq'])} bp  |  "
              f"Protein {len(ref_data[gene]['protein'])} aa  ({nt_rec.id})")

    # Configure aligner once (shared across all alignments)
    aligner = make_aligner()

    # Discover sample files
    print("\n[2] Discovering AB1 files …")
    samples, notes = discover_samples(args.raw)
    if not samples:
        print("  ERROR: No AB1 files found in", args.raw)
        sys.exit(1)

    # Sort sample IDs numerically (S2, S3, …, S53) rather than lexically.
    def sample_key(sid):
        m = re.match(r"S(\d+)", sid)
        return int(m.group(1)) if m else 1_000_000
    sorted_samples = sorted(samples, key=sample_key)

    n_combos = sum(len(v) for v in samples.values())
    print(f"   {len(samples)} samples, {n_combos} sample×hotspot combinations")
    if notes:
        print(f"   {len(notes)} parsing note(s):")
        for note in notes:
            print(f"     {note}")

    # Analyse each combination
    print("\n[3] Analysing reads …")
    analysis_date = datetime.now().strftime("%Y-%m-%d")
    rows = []

    for sample_id in sorted_samples:
        for (gene, raw_hotspot) in sorted(samples[sample_id]):
            pair    = samples[sample_id][(gene, raw_hotspot)]
            fw_path = pair.get("FW")
            rv_path = pair.get("RV")
            hotspot = canonical_hotspot(raw_hotspot)   # HS1b → HS1 for biology

            row = analyse_hotspot(
                sample_id     = sample_id,
                gene          = gene,
                hotspot       = hotspot,
                fw_path       = fw_path,
                rv_path       = rv_path,
                ref_nt_seq    = ref_data[gene]["nt_seq"],
                ref_protein   = ref_data[gene]["protein"],
                aligner       = aligner,
                hotspot_label = raw_hotspot,
                trim_threshold      = args.trim_q,
                trim_window         = args.trim_window,
                min_single_strand_q = args.min_single_strand_q,
                max_local_conflicts = args.max_local_conflicts,
                max_local_indels    = args.max_local_indels,
            )

            # Add static metadata
            hs_range = HOTSPOT_AA_RANGES[gene][hotspot]
            row["Hotspot_AA_Range"]  = f"AA {hs_range[0]}–{hs_range[1]}"
            row["Analysis_Date"]     = analysis_date
            row["Reference_Genome"]  = (
                f"GCF_010111755.1 | {ref_data[gene]['nt_id']} | "
                f"ATCC 2001 | Johns Hopkins 2020"
            )

            rows.append(row)

            # Compact one-line console summary (368 files → keep it readable)
            status = row.get("Call_Status", "N/A")
            hs_mut = row.get("AA_Mutations_in_Hotspot", "None")
            flag   = "  ⚠" if status != "No mutation" else ""
            print(f"   {sample_id:<5} {gene} {raw_hotspot:<5} → {status:<48} {hs_mut}{flag}")

    # Write CSV
    print("\n[4] Writing report …")
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path  = args.out / f"{args.prefix}_{timestamp}.csv"
    write_csv(rows, out_path)

    # Summary — counts + the rows that actually need attention
    print("\n" + "=" * 70)
    print("  SUMMARY")
    print("=" * 70)

    from collections import Counter
    tally = Counter(r.get("Call_Status", "N/A") for r in rows)
    print(f"  Total analyses: {len(rows)}")
    for status, n in sorted(tally.items(), key=lambda x: -x[1]):
        print(f"    {n:>3}  {status}")

    # Rows needing attention: mutations found, plus inconclusive/partial results
    actionable = [r for r in rows if r.get("Call_Status", "") != "No mutation"]
    print(f"\n  Rows needing attention: {len(actionable)}")
    if actionable:
        print(f"  {'Sample':<7}{'Gene':<6}{'HS':<6}{'AA in Hotspot':<24}{'Status'}")
        print(f"  {'-'*6:<7}{'-'*5:<6}{'-'*5:<6}{'-'*23:<24}{'-'*40}")
        for r in actionable:
            print(f"  {r['Sample_ID']:<7}{r['Gene']:<6}{r['Hotspot']:<6}"
                  f"{r.get('AA_Mutations_in_Hotspot','N/A')[:22]:<24}"
                  f"{r.get('Call_Status','N/A')}")
    print()


if __name__ == "__main__":
    main()
