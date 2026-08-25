# Coming from tRAX

What changes when you move a tRAX workflow to tRNAgraph. Skim the tables; read the notes only if a number doesn't look the way you expect.

---

## The three things that matter most

1. **Results live in one `.h5ad` AnnData object**, not a folder of text files. tRAX-format flat files are still written alongside it, so existing scripts keep working.
2. **Most outputs match tRAX numerically.** The exceptions are listed under [Where numbers differ](#where-numbers-differ) below.
3. **Several flags were renamed or removed.** See [Flag changes](#flag-changes).

---

## Bug fixes

Cases where tRAX's behaviour was wrong and tRNAgraph does not reproduce it.

| What | tRAX | tRNAgraph |
| --- | --- | --- |
| **Read classification order** | A read overlapping several features was assigned via Python `set` iteration — hash-order dependent, so results varied between runs | Candidates sorted deterministically; same input always gives the same answer |
| **`avgs.txt` columns** | Labelled per comparison (`A_B`, `A_C`, `B_C`) but every column holds the same value — DESeq2's `baseMean`, which ignores the contrast | `baseMean` emitted once under its own name; the per-group means moved to `groupavgs.txt`, one column per group |
| **`combine.txt` trailing columns** | Per-group medians, labelled with bare group names | Same — previously per-group means under the same labels, so a column-aligned diff compared a mean against a median |
| **`padjs.txt`** | Holds *unadjusted* p-values despite the name — `analyzecounts.R` takes column 5 of each `results()` object, and column 5 of a DESeq2 results table is `pvalue` | Holds genuine Benjamini-Hochberg `padj` |
| **`--mincoverage` scope** | Dropped low-count genes from the coverage file, which silently removed them from *every* downstream output | Renamed `--minfeaturereads`; affects only the VST dispersion-trend fit. Every gene keeps a full coverage row |
| **`aminocounts` / `anticodoncounts`** | Held unique-read counts, duplicating `unique/` — so no all-reads view existed | Main files hold all reads; `unique/` keeps the unique breakdown and still matches tRAX exactly |

---

## Improvements

| Area | tRAX | tRNAgraph | Why |
| --- | --- | --- | --- |
| Differential expression | R `DESeq2` via subprocess | `PyDESeq2`, in-process | Drops the R dependency, and results land straight in the AnnData object instead of being parsed back from text |
| Variance stabilisation | `rlog` | VST | DESeq2's own documentation recommends VST over rlog for anything but small sample counts, because rlog costs far more to compute |
| Trimming | `cutadapt` / `SeqPrep` | `fastp` | One tool covering adapter and quality trimming instead of two |
| Size factors | Computed from all features | Computed from tRNAs only. The all-feature set is still written, as `<exp>-allfeature_SizeFactors.txt` | Non-tRNA abundance can shift independently of tRNAs, which distorts tRNA normalisation when both go into the same reference set |
| Mismatch data | Per-position detail only via an R script | Stored at full per-position granularity in the AnnData object | Keeps the per-position detail queryable for misincorporation work, rather than only as a rendered summary |
| Read basis in plots | Each plot picked its own — some unique reads, some all reads, with nothing on the figure saying which | Every graph type uses unique (transcript-specific) reads by default; `--allreads` switches the whole command at once | Two plots of one dataset can no longer rest on different denominators without saying so |

---

## New features

- **Read-length splits** — `analyze build -c 60` produces `u60`/`o60` variants inside the same `.h5ad`, no separate runs.
- **Clustering** — `analyze cluster` adds UMAP/HDBSCAN coordinates to the object.
- **A plotting layer over the object** — coverage, PCA, volcano, heatmap, radar, seqlogo, correlation, count plots from `trnagraph graph`.
- **Selectable coverage specificity** — `--covtype unique|isodecoder|isotype|notamino|total` picks any of tRAX's four read-assignment categories, each in its own folder, plus a stacked plot of all four at once. tRAX computed the same partition but only ever surfaced it as one figure and a `unique/` folder.
- **A built-in demo pipeline** — `trnagraph tools test` runs end to end on a bundled dataset.

---

## Migration

### Flag changes

| tRAX / old tRNAgraph | Now | Note |
| --- | --- | --- |
| `--uniqueonly` | *(removed)* | Multi-mapped reads are always excluded from coverage, which is what tRAX did in practice |
| `--mincoverage` | `--minfeaturereads` | Also re-scoped — see Bug fixes |
| `--dumpother` | `--filterother` | Rename only |
| `--lazy` | `--force-remap` | Sense inverted: remapping is now the thing you opt into |
| `--nofrag` | *(removed)* | |
| `--nosizefactors` | *(removed)* | Was broken |
| `--maponly` | *(removed)* | From `tools test` |
| `--diffrts total_unique` | `--diffrts total` + `--allreads` | The read basis moved out of the readtype and onto one command-wide flag, so graph types cannot disagree |
| `--pcareadtypes total_unique total` | `--pcareadtypes total` | Same change. The PCA and volcano overview pages still show both bases side by side |

`--skip-env-check` and `--skip-update-check` are *global* options — they go before the subcommand (`trnagraph --skip-env-check analyze build ...`).

### Outputs tRAX produces that tRNAgraph does not

| File | Why |
| --- | --- |
| `trimindex.txt` | fastp uses a different manifest convention |
| `Rlog-<exp>.txt` | Captured R subprocess output; there is no R subprocess |
| `positionmismatches.txt` | Was generated by an R script. The data is in the AnnData object instead |

### Where numbers differ

| Output | Difference | Why |
| --- | --- | --- |
| `typecounts`, `readlengths`, `anticodoncounts`, `aminocounts` | ~0.1–0.5% of reads | Deterministic classification replacing tRAX's hash-order-dependent version |
| `dispersions`, some `combine.txt` values | Small | PyDESeq2 vs. R DESeq2 |
| `SizeFactors.txt` | Different basis | tRNA-controlled rather than all-feature. Compare against `<exp>-allfeature_SizeFactors.txt` for a like-for-like check |
| `coverage`, `pretRNAcoverage` | Extra features present | tRAX dropped genes below `--mincoverage`; tRNAgraph keeps them |
| BAMs and everything downstream | Varies | fastp vs. cutadapt/SeqPrep produce slightly different read sets |

---

_Full flag and option reference: `cli_reference.md`. AnnData schema: `data_structure.md`._
