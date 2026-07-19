# FKS Sanger Pipeline

Automated detection of *FKS1* and *FKS2* hotspot mutations from paired Sanger
sequencing chromatograms in **_Nakaseomyces glabratus_** (formerly *Candida glabrata*).

The *FKS1* and *FKS2* genes encode subunits of β-1,3-glucan synthase, the target of
the echinocandin antifungals. Mutations cluster in two short, conserved "hotspot"
windows in each gene. This pipeline reads raw `.ab1` traces, builds a
quality-controlled forward/reverse consensus, translates it, and reports every
nucleotide and amino-acid difference from the reference.

> **Scope.** This pipeline **detects and reports mutations only.** It performs no
> clinical interpretation — no resistance calls, drug susceptibility, or MIC
> assignment. Interpretation is done separately against a curated database, so
> that reference data has a single maintained source of truth.

---

## What it does

Per sample × hotspot, the pipeline:

1. **Reads** the forward and reverse `.ab1` chromatograms (base calls + Phred scores).
2. **Trims** low-quality read ends with a sliding window (default: mean Q ≥ 20 over 10 bases).
3. **Reverse-complements** the reverse read and aligns both to the reference CDS.
4. **Builds a consensus**, resolving forward/reverse disagreements by Phred score.
5. **Calls variants**, translating each fully-covered codon (standard genetic code).
6. **Applies quality gates** and reports the result as evaluable or inconclusive.
7. **Exports** a CSV with one row per sample × hotspot.

### Quality control

Three safeguards prevent over-confident calls — they judge whether the *sequencing*
is trustworthy, not whether a mutation is clinically significant:

| Guard | Rule | Rationale |
|---|---|---|
| **Single-strand confidence** | A base covered by only one strand must reach Q30 to be called; otherwise masked | Positions covered by both strands cross-validate each other; a lone read has no such check |
| **Coverage guard** | A "no mutation" result is downgraded if the hotspot is not fully covered; 0% coverage is reported as inconclusive | Absence of a mutation across unsequenced codons is not evidence of absence |
| **Reliability gate** | A hotspot with indels or excessive forward/reverse discordance inside the window is reported as inconclusive | Indel/discordance clusters are the signature of misalignment in a low-quality read, and generate artefactual substitutions and premature stop codons |

Masked and excluded calls are retained in dedicated CSV columns, so nothing is
silently discarded.

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

## Output

A CSV with one row per sample × hotspot. Key columns:

| Column | Meaning |
|---|---|
| `Call_Status` | `Mutation detected` · `No mutation` · `No mutation (partial: X%)` · `INCONCLUSIVE — …` |
| `AA_Mutations_in_Hotspot` | Amino-acid changes inside the hotspot, e.g. `S663P` |
| `AA_Mutations_outside_Hotspot` | Changes in the flanking sequenced region |
| `NT_Mutations` | Nucleotide substitutions, e.g. `T1987C` (CDS numbering) |
| `Synonymous_NT_Changes` | Silent changes with no amino-acid consequence |
| `Hotspot_Coverage` | Callable codons / total codons in the window |
| `Hotspot_Local_Conflicts`, `Hotspot_Local_Indels` | Reliability-gate inputs |
| `NT_Masked_LowConf`, `AA_Masked_LowConf` | Calls suppressed for low single-strand quality |
| `QC_Flags` | Human-readable quality notes |

Notation: `S663P` = reference Serine at protein position 663 replaced by Proline.
`T1987C` = reference T at CDS nucleotide 1987 replaced by C.

---

## Example output

Running on the bundled example data:

```
   S49   FKS2 HS1   → Mutation detected      S663P
```

This isolate carries a serine-to-proline substitution at *FKS2* codon 663
(`T1987C`), with 37/37 hotspot codons covered.

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
