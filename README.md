# FKS Sanger Pipeline

Automated detection of *FKS1* and *FKS2* hotspot mutations from paired Sanger
sequencing chromatograms in **_Nakaseomyces glabratus_** (formerly *Candida glabrata*).

The *FKS1* and *FKS2* genes encode subunits of β-1,3-glucan synthase, the target of
the echinocandin antifungals. Mutations cluster in two short, conserved "hotspot"
windows in each gene. This pipeline reads raw `.ab1` traces, builds a
quality-controlled forward/reverse consensus, translates it, and reports
nucleotide and amino-acid differences from the reference.

> **Scope.** This pipeline **detects and reports mutations only.** It performs no
> clinical interpretation — no resistance calls, drug susceptibility, or MIC
> assignment.

> **Reporting standard (haploid).** *N. glabratus* is haploid, so a true mutation is
> homozygous and must be independently confirmed on **both** strands. A variant is
> reported only when the forward and reverse reads **agree on the variant base at the
> mutated position**, with no local indels or forward/reverse discordance. Single-strand
> or strand-discordant calls are kept in the audit output but are **not** reported as
> findings.

---

## What it does

Per sample × amplicon, the pipeline:

1. **Reads** the forward and reverse `.ab1` chromatograms (base calls + Phred scores).
2. **Trims** low-quality read ends with a sliding window (default: mean Q ≥ 20 over 10 bases).
3. **Reverse-complements** the reverse read and aligns both to the reference CDS.
4. **Builds a consensus**, resolving forward/reverse disagreements by Phred score.
5. **Calls variants**, translating each fully-covered codon (standard genetic code).
6. **Classifies** each mutation by gene (FKS1/FKS2) and region (HS1/HS2/Outside).
7. **Scores** each mutation for bidirectional confirmation and local reliability.
8. **Exports** three CSV reports (below).

### Reports

| File | Contents |
|---|---|
| `<prefix>_<ts>.csv` | Amplicon-level QC report — one row per sample × amplicon |
| `<prefix>_by_region_<ts>.csv` | Every called mutation, with FKS/HS classification and confirmation flags (full audit trail) |
| `<prefix>_confirmed_<ts>.csv` | **The reportable set** — both-strand-confirmed mutations from evaluable amplicons only |

Key columns in the mutation reports: `Both_Strand_Confirmed` (Yes/No),
`Strand_Support` (both/single/conflict), `Local_Conflicts`, `Local_Indels`,
`Edge_Dist_bp`, `Isolate_Count` (recurrence), `Confidence` (High/Moderate/Low).

### Quality control

Safeguards judge whether the *sequencing* is trustworthy, not whether a mutation is
clinically significant:

| Guard | Rule | Rationale |
|---|---|---|
| **Bidirectional confirmation** | A mutation is reported only if the forward and reverse reads agree on the variant base at the mutated position | *N. glabratus* is haploid; a real mutation is homozygous and must appear on both strands. This is the governing filter for the reportable set |
| **Coverage guard** | A "no mutation" result is downgraded if the hotspot is not fully covered; 0% coverage → inconclusive | Absence of a mutation across unsequenced codons is not evidence of absence |
| **Reliability gate** | A hotspot with indels or excessive forward/reverse discordance is reported as inconclusive | Indel/discordance clusters are the signature of misalignment in a low-quality read, and generate artefactual substitutions and premature stops |

An internal consensus-quality step masks single-strand positions below Q30
(`--min-single-strand-q`, default 30); this keeps coverage metrics and any deposited
consensus sequences clean, but does **not** determine which mutations are reported —
bidirectional confirmation does. Nothing is silently discarded: every call is retained
in the by-region audit report with its confirmation flags.

---

## Installation

Requires **Python 3.9+**.

```bash
git clone https://github.com/USERNAME/fks-sanger-pipeline.git
cd fks-sanger-pipeline
pip install -r requirements.txt
```

## Usage

```bash
# Run on the bundled example (one isolate carrying FKS2 S663P)
python sanger_pipeline.py --raw example_data --out results

# Run on your own data
python sanger_pipeline.py --raw /path/to/ab1_files --out /path/to/results

# See all options
python sanger_pipeline.py --help
```

### Options

| Option | Default | Description |
|---|---|---|
| `--raw` | `./Raw` | Directory of paired `.ab1` files |
| `--refs` | `./references` | Directory of reference FASTA files |
| `--out` | `./Results` | Output directory for the CSV report |
| `--prefix` | `FKS_Mutation_Report` | Output filename prefix |
| `--trim-q` | `20` | Minimum mean Phred for end trimming |
| `--trim-window` | `10` | Sliding-window size for end trimming |
| `--min-single-strand-q` | `30` | Minimum Phred to call a single-strand position |
| `--max-local-conflicts` | `3` | Max in-hotspot FW/RV conflicts before inconclusive |
| `--max-local-indels` | `0` | Max in-hotspot indels before inconclusive |

---

## Input file naming

Files must be named so the sample, gene, hotspot, and read direction can be parsed.
Both orderings are accepted, with `-` or `_` before the direction:

```
S49-FKS2-HS1_FW.ab1        # sample first
FKS2-HS1b-S49_FW.ab1       # sample last
```

- **Gene**: `FKS1` or `FKS2`
- **Hotspot**: `HS1` or `HS2`, optionally with a primer-batch suffix (`HS1b`),
  which is normalised to the canonical hotspot for analysis but preserved in the report
- **Direction**: `FW` or `RV` (the variants `F`, `R`, `W` are accepted and corrected, with a note)

Unparseable filenames are skipped and listed in the run log.

---

## Reference sequences

Bundled in `references/`, retrieved from NCBI:

| | |
|---|---|
| Organism | *Nakaseomyces glabratus*, strain ATCC 2001 |
| Assembly | `GCF_010111755.1` (RefSeq) / `GCA_010111755.1` (GenBank), ASM1011175v1 |
| Submitter | Johns Hopkins Medical Institution, 3 February 2020 |
| *FKS1* | Chromosome G, `NC_088957.1`:107061–112652 — 5,592 bp / 1,863 aa |
| *FKS2* | Chromosome K, `NC_088961.1`:375102–380795 — 5,694 bp / 1,897 aa |

Each reference FASTA begins at the start codon, so codon numbering is
unambiguous. The in-frame translation of each CDS was confirmed to match the
corresponding reference protein.

### Hotspot windows (amino-acid coordinates, verified against ATCC 2001)

| Gene | HS1 | HS2 |
|---|---|---|
| *FKS1* | 620–660 | 1335–1370 |
| *FKS2* | 654–690 | 1370–1400 |

---

## Output notation

`S663P` = reference Serine at protein position 663 replaced by Proline.
`T1987C` = reference T at CDS nucleotide 1987 replaced by C.

The **amplicon-level report** carries per-amplicon QC (`Call_Status`,
`Hotspot_Coverage`, `Hotspot_Local_Conflicts/Indels`, `NT_Mutations`, `QC_Flags`, …).
The **mutation reports** carry one row per called change with the classification and
confirmation columns listed under [Reports](#reports) above. To obtain the reportable
mutation set, use the `_confirmed_` file (or filter the by-region file to
`Both_Strand_Confirmed == Yes` and `Amplicon_QC == Evaluable`).

---

## Example output

Running on the bundled example data:

```
   S49   FKS2 HS1   → Mutation detected      S663P
```

The `_confirmed_` report contains one row — *FKS2* S663P (`T1987C`), confirmed on
both strands (`Both_Strand_Confirmed = Yes`, `Confidence = High`), at 37/37 hotspot
codons covered.

---

## Dependencies

- [Biopython](https://biopython.org/) — `.ab1` parsing, sequence handling, translation, pairwise alignment

## Citation

If you use this pipeline, please cite the accompanying publication (details to
follow) and Biopython:

> Cock PJA, Antao T, Chang JT, et al. Biopython: freely available Python tools for
> computational molecular biology and bioinformatics. *Bioinformatics.*
> 2009;25(11):1422–1423. doi:10.1093/bioinformatics/btp163

## License

MIT — see [LICENSE](LICENSE).
