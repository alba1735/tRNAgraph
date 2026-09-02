# Tutorial

Three worked examples, in increasing order of complexity, covering the whole pipeline from raw
reads to figures. Each one is a complete run you can adapt: the first is a real public dataset
you can download and reproduce exactly, the other two are templates you point at your own FASTQs.

If you only want to confirm your installation works, stop here and run `trnagraph tools test`,
which downloads a small bacterial dataset and runs everything end to end without any setup. This
document is about running the pipeline on _your_ data.

| Example                               | Organism                  | Method   | What it adds                                                      |
| ------------------------------------- | ------------------------- | -------- | ----------------------------------------------------------------- |
| [1](#example-1-the-basic-run)         | _S. cerevisiae_ (sacCer3) | ARM-seq  | The five core commands, start to finish                           |
| [2](#example-2-a-more-considered-run) | _M. musculus_ (mm10)      | ARM-seq  | Reference choice, non-tRNA features, grouped comparisons, styling |
| [3](#example-3-the-full-treatment)    | _H. sapiens_ (hg38)       | OTTR-seq | UMIs and deduplication, read-length splits, config-driven figures |

---

## A note on library preparation

tRNAs are hard to sequence for two reasons: they are heavily modified, and they are stably
folded. Both interfere with reverse transcription, so a standard small-RNA protocol
systematically under-reports exactly the molecules you care about. Different methods attack this
differently, and which one produced your data changes which tRNAgraph outputs are worth reading.

| Method                                                        | Approach                                                                        | What it means for your run                                                                                                                                                                                                                          |
| ------------------------------------------------------------- | ------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Plain small RNA-seq                                           | No special treatment                                                            | Modified positions cause reverse transcription to stop, so tRNA fragments are under-counted and 3′-biased. Nothing in tRNAgraph needs changing; just do not read the fragment landscape as biology.                                                 |
| **ARM-seq** ([DOI](https://doi.org/10.1038/nmeth.3508))       | AlkB demethylase pretreatment removes m¹A, m³C and m¹G before library prep      | Run demethylated and untreated samples as two groups. The comparison _is_ the experiment, so put it in your metadata and give it to `--volgrp`. `-g mismatch` is the readout. Protocol chapter: [DOI](https://doi.org/10.1007/978-1-4939-6807-7_15) |
| **DM-tRNA-seq** ([DOI](https://doi.org/10.1038/nmeth.3478))   | Demethylase _plus_ a thermostable group II intron reverse transcriptase (TGIRT) | Same as ARM-seq from tRNAgraph's point of view, but the processive RT also reads through structure, so expect more full-length reads. Worth a read-length split (Example 3) to separate them.                                                       |
| **OTTR-seq** ([DOI](https://doi.org/10.1073/pnas.2107900118)) | A retroelement RT that jumps templates, fusing both adapters in one reaction    | Low bias and genuinely end-to-end, so 3′ ends are captured faithfully — the CCA end-type breakdown is trustworthy. Often carries UMIs, so see Example 3.                                                                                            |

These are four common choices, not a complete list. Every RNA-seq method and every experiment
differs; treat the table as a starting point for deciding which outputs mean something for your
data, not as a classification your dataset must fit into.

A review covering the demethylase-based methods together: [DOI](https://doi.org/10.1016/j.molmed.2016.10.009).

---

## Before you start: getting your references

Every run needs four reference inputs for `preprocess makedb`:

| Input                   | Flag | Where it comes from                |
| ----------------------- | ---- | ---------------------------------- |
| Genome FASTA            | `-g` | UCSC, Ensembl, ENCODE, NCBI        |
| tRNAscan-SE predictions | `-t` | gtRNAdb, in the per-genome tarball |
| tRNA FASTA              | `-r` | gtRNAdb, same tarball              |
| Name map                | `-m` | gtRNAdb, same tarball              |

Three of the four come from one gtRNAdb download, e.g.
`https://gtrnadb.ucsc.edu/genomes/eukaryota/Scere3/sacCer3-tRNAs.tar.gz`. Unpack it and you have
everything except the genome.

### Gotcha 1: the prediction file is named differently per genome

gtRNAdb does not use one filename across genomes, and the variants are not interchangeable.
Look in the tarball before assuming:

| Genome  | File to pass to `-t`                                                      |
| ------- | ------------------------------------------------------------------------- |
| sacCer3 | `sacCer3-tRNAs.out-noChrM`                                                |
| hg38    | `hg38-tRNAs-detailed.out`                                                 |
| mm10    | `mm10-tRNAs-confidence-set.out` — **not** the detailed one; see Example 2 |

### Gotcha 2: chromosome names must match between genome and GTF

If you pass a GTF to `analyze build --gtf`, its chromosome names have to match the genome you
built the database from. UCSC and ENCODE use `chr1`, `chr2`, `chrM`; Ensembl uses `1`, `2`, `MT`.
Mixing them is easy and the consequences are quiet:

> [!WARNING]
> A chromosome-name mismatch does **not** produce an error. The GTF parses, the features load,
> and then no read ever matches one — so every non-tRNA count comes back zero and the non-tRNA
> volcano and PCA panels silently disappear. The only hint is a message at graph time saying
> "No non-tRNA feature counts found … re-run with `--gtf`", which is actively misleading, because
> you _did_ pass `--gtf`. The tRNA side of the pipeline is driven by the gtRNAdb database rather
> than the GTF, so it looks perfectly healthy throughout.

The fix is to pick one convention and stay in it. Check before you run:

```bash
head -1 genome.fa                          # >chr1  or  >1
awk '!/^#/ {print $1; exit}' features.gtf  # chr1    or   1
```

If they disagree, either fetch the matching pair, or rename one — for an Ensembl GTF against a
UCSC genome:

```bash
awk 'BEGIN{OFS="\t"} /^#/ {print; next} {$1 = ($1=="MT" ? "chrM" : "chr"$1); print}' \
  ensembl.gtf > ensembl.chr.gtf
```

### Adapters

`preprocess trim` autodetects adapters, and for most small-RNA libraries that works. Pin them
with `-a1`/`-a2` when you know what they are, so a run is reproducible rather than dependent on
what fastp infers from a given batch of reads.

For libraries prepared with the standard Illumina small-RNA chemistry — which covers ARM-seq and
the OTTR libraries in Example 3 — the sequences are:

| Read                 | Sequence                | Flag  |
| -------------------- | ----------------------- | ----- |
| R1 / single-end      | `AGATCGGAAGAGCACACGTC`  | `-a1` |
| R2 (paired-end only) | `GATCGTCGGACTGTAGAACTC` | `-a2` |

These are tRAX's own defaults (`trimadapters.py`), which is good provenance: tRAX comes from the
lab that developed ARM-seq and shares authors with that paper, so its defaults encode the
protocol actually in use rather than a generic guess. The ARM-seq paper itself names the library
kit but does not print the sequence.

They also check out against real data. Running fastp in autodetect mode on Example 1's reads
reports the adapter as `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` — the R1 sequence above, followed by
the TruSeq index region. Pinning the shorter form trims marginally more reads (87,815 against
87,254 of 100,000), because a shorter sequence matches more permissively.

You can run that check on your own data before deciding:

```bash
fastp -i yourreads.fastq.gz -o /dev/null -j report.json -h /dev/null
python -c "import json; print(json.load(open('report.json'))['adapter_cutting'])"
```

A wrong pinned adapter is worse than autodetection, so confirm before assuming these apply.

---

## Example 1: the basic run

**Dataset:** _Saccharomyces cerevisiae_ BY4741, ARM-seq, from the paper that introduced the
method (Cozen et al., [DOI](https://doi.org/10.1038/nmeth.3508)). Six runs under SRA accession
**SRP056032**: three AlkB-treated, three untreated.

This example is worth running even if yeast is not your organism: it is small, it finishes in
minutes, and it reproduces a published result, so you can tell whether your installation is
behaving before you trust it on data you care about.

### Get the data

```bash
mkdir -p references raw config && cd references
curl -O https://gtrnadb.ucsc.edu/genomes/eukaryota/Scere3/sacCer3-tRNAs.tar.gz
tar xzf sacCer3-tRNAs.tar.gz
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
gunzip sacCer3.fa.gz
cd ..

# 100,000 reads per run is plenty to see the effect and keeps the run short.
for acc in SRR1874029 SRR1874032 SRR1874034 SRR1874045 SRR1874048 SRR1875056; do
  prefetch $acc --output-file raw/$acc.sra
  fastq-dump -O raw -X 100000 --gzip raw/$acc.sra
done
```

### Write the two input files

`config/manifest.tsv` drives trimming. First column is the sample name _and_ the output
location; second is the FASTQ:

```text
AlkB_1	raw/SRR1874029.fastq.gz
AlkB_2	raw/SRR1874032.fastq.gz
AlkB_3	raw/SRR1874034.fastq.gz
untreated_1	raw/SRR1874045.fastq.gz
untreated_2	raw/SRR1874048.fastq.gz
untreated_3	raw/SRR1875056.fastq.gz
```

`config/metadata.tsv` drives everything after trimming. It needs `fastq`, `sample` and `group`;
any further columns become groupable dimensions later:

```text
fastq	sample	group	treatment
processed/trimmed/AlkB_1_trimmed.fastq.gz	AlkB_1	AlkB	demethylated
processed/trimmed/AlkB_2_trimmed.fastq.gz	AlkB_2	AlkB	demethylated
processed/trimmed/AlkB_3_trimmed.fastq.gz	AlkB_3	AlkB	demethylated
processed/trimmed/untreated_1_trimmed.fastq.gz	untreated_1	untreated	mock
processed/trimmed/untreated_2_trimmed.fastq.gz	untreated_2	untreated	mock
processed/trimmed/untreated_3_trimmed.fastq.gz	untreated_3	untreated	mock
```

> [!IMPORTANT]
> `preprocess trim` writes a `processed/trimmed/trim_metadata.tsv` template for you, but it sets
> `group` equal to `sample` — six groups of one, which makes differential expression meaningless.
> Edit it, or write the file yourself as above. This is the single most common way a first run
> produces nothing useful.

### Run the pipeline

```bash
trnagraph preprocess makedb \
  -g references/sacCer3.fa \
  -t references/sacCer3-tRNAs.out-noChrM \
  -r references/sacCer3-tRNAs.fa \
  -m references/sacCer3-tRNAs_name_map.txt \
  -s euk -o references/sacCer3

trnagraph preprocess trim -i config/manifest.tsv

trnagraph preprocess map -i config/metadata.tsv -d references/sacCer3 -o sacCer3_demo

trnagraph analyze build -i config/metadata.tsv -d references/sacCer3 -o sacCer3_demo

trnagraph graph -i sacCer3_demo/sacCer3_demo.h5ad -o sacCer3_demo/graphs -g all
```

No adapter is given to `trim`: fastp detects it, which is the safe default when you are not
certain what your libraries carry. If you do know, pin it instead of inferring it — for this
dataset that would be `-a1 AGATCGGAAGAGCACACGTC`, the standard Illumina 3′ adapter and the
default tRAX used for exactly this kind of library (see [Adapters](#adapters) below).

`-s euk` matters. It selects the covariance model and Sprinzl position table for the organism
domain; `bact` and `arch` are the other options, and the wrong one mislabels positions.

### What you should see

Around 98–99% of reads align, and `graph -g all` writes roughly 300 files across ten graph types.
The result to check is the one the paper reports — demethylation should recover tRNA reads that
reverse transcription otherwise loses:

```bash
python - <<'PY'
import anndata as ad, pandas as pd
a = ad.read_h5ad("sacCer3_demo/sacCer3_demo.h5ad")
reads = a.obs.groupby("sample", observed=True)["nreads_total_unique_raw"].sum()
group = a.obs.drop_duplicates("sample").set_index("sample")["group"]
means = pd.DataFrame({"reads": reads, "group": group}).groupby("group")["reads"].mean()
print(means, f"\nratio: {means['AlkB'] / means['untreated']:.2f}x")
PY
```

On the 100,000-read subsample above this gives roughly **2.4×** more tRNA reads in the treated
samples. The paper reports the tRNA share of yeast reads rising from 6.9% to 15.1%, i.e. about
2.2×, so a subsample landing near that is a good sign your installation is sound.

You will also see this at graph time:

```text
No non-tRNA feature counts found in AnnData object … Re-run with --gtf to enable these
```

That is correct here — no `--gtf` was passed, so nothing outside the tRNA database was
annotated and every non-tRNA read fell into the `other` category. Example 2 adds it. (This is the
same message that appears in the chromosome-naming trap above, where it is a lie; the difference
is whether you passed `--gtf` at all.)

### Paired-end data

The yeast run above is single-end, but not because ARM-seq is: those libraries were sequenced
paired-end on a MiSeq and the mates were merged before deposit, so what SRA serves is one merged
read per spot. Your own data will usually still be in two files, so it is worth knowing what
tRNAgraph does with them.

Add a third column to the manifest and the sample becomes paired-end:

```text
AlkB_1	raw/AlkB_1_R1.fastq.gz	raw/AlkB_1_R2.fastq.gz
```

Three things then change:

- **fastp merges overlapping mates** into one read, and enables its paired-end adapter detection
  automatically. For tRNA-seq this is what you want: fragments are short, mates overlap heavily,
  and a merged read is the actual molecule. To pin both adapters rather than detect them, pass
  `-a1 AGATCGGAAGAGCACACGTC -a2 GATCGTCGGACTGTAGAACTC` — the same pair tRAX hands to SeqPrep for
  paired-end input.
- **The output is named `_merged.fastq.gz`**, not `_trimmed.fastq.gz`. Your `metadata.tsv` must
  point at that name — this catches people out.
- **Pairs that could not be merged are written to `_unmerged_R1.fastq.gz` / `_unmerged_R2.fastq.gz`
  and go no further.** They are kept so you can inspect them, but they are not recorded in
  `trim_metadata.tsv` and nothing downstream reads them. A large unmerged fraction means your
  inserts are longer than your read length allows for overlap, which is worth knowing.

`--style` is also most informative here: its `colors.trimtype` block controls the trimming QC
plot, and the `Merged`/`Unmerged` categories only exist for paired-end samples. Single-end
samples show `Trimmed` and `Discarded` instead.

---

## Example 2: a more considered run

**Dataset:** your own mouse ARM-seq libraries, treated versus untreated. Everything below is a
template — substitute your FASTQ paths.

Three things change relative to Example 1: the reference set is chosen deliberately, non-tRNA
features are annotated, and the figures are grouped and coloured by the experimental variable
rather than by sample name.

### Choosing the mm10 reference set

This is the part worth reading even if you are not working in mouse.

Rodent genomes contain enormous numbers of SINEs — short interspersed elements — that are
_derived from tRNA genes_ and still look like them to tRNAscan-SE. Passing the full prediction
set builds a database mostly composed of things that are not tRNAs, and every real tRNA read then
multi-maps across hundreds of decoys. Counted in the gtRNAdb mm10 tarball:

| File                                              | Predictions |
| ------------------------------------------------- | ----------- |
| `mm10-tRNAs-detailed.out` (all)                   | **40,912**  |
| `mm10-tRNAs-confidence-set.out` (high confidence) | **408**     |

A hundredfold difference. The composition says where it comes from: **25,655 of the 40,912 —
63% — are annotated Ser**, which is the signature of the mouse B2 SINE family, itself derived
from tRNA-Ser. Another 4,358 are `Undet`. In the high-confidence set Ser is no longer dominant at
all (the largest family is Cys, at 54). For contrast, sacCer3 has 275 predictions in total and no
such problem — yeast simply does not carry these elements.

So for mouse and rat, use the high-confidence set:

```bash
trnagraph preprocess makedb \
  -g references/mm10.fa \
  -t references/mm10-tRNAs-confidence-set.out \
  -r references/mm10-tRNAs.fa \
  -m references/mm10-tRNAs_name_map.txt \
  -s euk -o references/mm10
```

This is a workaround inherited from tRAX rather than a solution — the high-confidence filter
removes real tRNA genes along with the decoys. It is the best available option today; dedicated
SINE handling is on the roadmap.

### Adding non-tRNA features

Passing a GTF to `analyze build` annotates everything outside the tRNA database, which turns the
undifferentiated `other` pile from Example 1 into named categories and unlocks two extra panels
in the volcano and PCA outputs automatically.

```bash
trnagraph preprocess trim -i config/manifest.tsv -a1 AGATCGGAAGAGCACACGTC
trnagraph preprocess map -i config/metadata.tsv -d references/mm10 -o mm10_armseq

trnagraph analyze build \
  -i config/metadata.tsv -d references/mm10 -o mm10_armseq \
  --gtf references/Mus_musculus.GRCm38.102.chr.gtf
```

> [!IMPORTANT]
> Re-read Gotcha 2 before this step. A UCSC mm10 genome (`chr1`) with a stock Ensembl GRCm38 GTF
> (`1`) produces zero non-tRNA counts and no error. Either convert the GTF's names, or take the
> genome and annotation from the same source.

### Metadata that earns its keep

The extra columns are the point. With a `treatment` column, the comparison you actually ran
becomes the axis your figures are drawn on:

```text
fastq	sample	group	treatment	tissue
processed/trimmed/liver_alkb_1_trimmed.fastq.gz	liver_alkb_1	AlkB	demethylated	liver
processed/trimmed/liver_mock_1_trimmed.fastq.gz	liver_mock_1	untreated	mock	liver
```

### Checking what you can group by

Every grouping option below names a column that lives inside the `.h5ad`. Rather than guessing
one, ask the object:

```bash
trnagraph tools info -i mm10_armseq/mm10_armseq.h5ad
```

That lists every `obs` and `var` column with the values it holds, the `uns` keys, and the
`--variant` strings the object can resolve, so you can see that `treatment` really is a column
and that its values are `demethylated` and `mock`. Add `--column treatment` to print one
column's values in full, or `--json` to read it from a script.

This matters because a grouping option that names a column the object does not have **aborts
the run** rather than falling back to something else -- the run reports every bad label at once,
with near matches. A silently substituted column would give you a complete, plausible set of
figures grouped by the wrong thing.

### Graphing with grouping and colour

```bash
trnagraph graph \
  -i mm10_armseq/mm10_armseq.h5ad \
  -o mm10_armseq/graphs \
  -g all -g mismatch \
  --volgrp treatment \
  --pcacolors treatment --pcamarkers sample \
  --covgrp treatment \
  --style config/style.json
```

with `config/style.json`:

```json
{
  "colors": {
    "treatment": { "demethylated": "#FF007F", "mock": "#007FFF" }
  },
  "defaults": { "format": "pdf", "dpi": 300 },
  "volcano": { "rasterize_over": 2000 },
  "pca": { "rasterize_over": 2000 }
}
```

The `colors` block is keyed on the `obs` column name, so `treatment` here has to match what
`--volgrp`/`--pcacolors` name. `rasterize_over` keeps the point layer of a dense plot raster
while leaving text and axes vector — a fully vector PDF carrying thousands of markers is slow to
open and often rejected on submission.

`-g mismatch` is the one to look at for a demethylase experiment. AlkB removes the methyl groups
that cause reverse transcriptase to misincorporate, so the treated samples should show reduced
misincorporation exactly at the modified positions. That difference is the experiment; everything
else is context.

---

## Example 3: the full treatment

**Dataset:** your own human OTTR-seq libraries with UMIs, across an ordered condition such as a
timecourse. Template again.

Everything from Example 2 applies. What is new: UMIs and deduplication, read-length splits, and
figures driven entirely from config files so a whole figure set regenerates from one command.

### Trimming with UMIs

OTTR libraries commonly carry a UMI at the 5′ end of the read. One `trim` call handles the
adapter and the UMI together:

```bash
trnagraph preprocess trim \
  -i config/manifest.tsv \
  -a1 AGATCGGAAGAGCACACGTC \
  -u 7 \
  -l 15
```

`-u 7` moves a 7 nt UMI off the read and into the read name (`@READ_AAAACAC`). Use `--umi3` if
yours is at the 3′ end instead. The adapter is pinned rather than autodetected because it is
known — see [Adapters](#adapters) above for where that sequence comes from and how to confirm it
against your own reads.

### Mapping with deduplication

```bash
trnagraph preprocess map \
  -i config/metadata.tsv -d references/hg38 -o hg38_ottr \
  --dedup \
  --keep-prededup
```

`--dedup` runs `umi_tools dedup` over the mapped BAMs after mapping finishes. The deduplicated
reads take over the ordinary `<sample>.bam` name, so nothing downstream needs to know it
happened; `--keep-prededup` additionally keeps `<sample>.prededup.bam` so you can compare
deduplicated against non-deduplicated results without paying to map twice.

> [!IMPORTANT]
> `--dedup` refuses to run if it cannot find a UMI in the read names. This is deliberate:
> `umi_tools` does not fail on UMI-less input, it collapses reads sharing an alignment position
> instead — and in tRNA-seq, where transcripts are short and deeply covered, those reads are
> overwhelmingly real molecules. It would delete genuine signal and leave no trace.

For a protocol this wrapper does not cover — paired UMIs, a cell barcode, an unusual read-name
layout — `umi_tools` is in the conda environment already. Run it directly on the BAM directory
and continue with `analyze build`; the wrapper is a convenience, not a replacement.

### Building with a read-length split

```bash
trnagraph analyze build \
  -i config/metadata.tsv -d references/hg38 -o hg38_ottr \
  --gtf references/hg38_small_ncRNAs.gtf \
  --vst vst \
  -c 60
```

`-c 60` adds `u60` and `o60` variants — fragments under 60 nt and full-length molecules over it —
into the _same_ `.h5ad` as extra layers, each with independently fitted size factors. No separate
files are produced. This is what lets you ask whether a change is in the mature tRNA pool, in the
fragment pool, or both. Add further cutoffs later with `analyze addsplit -c 50` without
disturbing what is already there.

The GTF here is a curated small-ncRNA subset rather than the full Ensembl annotation. At human
scale that is a deliberate choice: a complete GTF pulls every protein-coding feature into a
comparison that is about non-coding RNA, which buries the signal and slows the build. Subset to
the feature classes your question is about.

### Clustering and graphing every variant

Clustering runs once per variant, chained through the same file so each pass accumulates rather
than overwriting:

```bash
CLUSTERED=hg38_ottr/hg38_ottr.cluster.h5ad
input=hg38_ottr/hg38_ottr.h5ad
for variant in norm:full norm:u60 norm:o60; do
  trnagraph analyze cluster -i $input -o $CLUSTERED --variant $variant -r 42
  input=$CLUSTERED
done

for variant in norm:full norm:u60 norm:o60; do
  trnagraph graph -i $CLUSTERED -o hg38_ottr/graphs \
    --variant $variant \
    --style config/style.json \
    --config config/filter.json
done
```

`-r 42` seeds UMAP. Without it the embedding differs between runs, which matters if the figure is
going into a paper that may be revised. Non-default variants are written into their own
`norm_u60/` and `norm_o60/` subfolders automatically, so one output directory serves all three.

`config/filter.json` restricts what is plotted without rebuilding anything:

```json
{
  "name": "no_undetermined",
  "obs_r": { "amino": ["Und", "Sup"] }
}
```

`obs` includes only the listed values; `obs_r` excludes them. `name` becomes a subfolder, so
several filtered views can coexist under one output directory.

### A note on 3′ ends

OTTR's reverse transcriptase adds non-templated nucleotides to 3′ ends, which raises a reasonable
worry: if reads systematically carry an extra base, does the CCA end-type classification shift by
one? `getendtype()` assigns end types by exact alignment end, with no tolerance, so a
one-base shift would silently relabel every `CCA` read as `CC`.

Measured on a real human OTTR dataset, it does not: **99.66%** of reads end exactly at the
reference 3′ end (`CCA`), against 0.04% `CC` and 0.02% `C`. The extra base does not survive into
the aligned block. So `--ccatail` and the `readends` coverage type mean what they say for OTTR
data, and the end-type breakdown is trustworthy — which is what you would hope from a method
designed to capture ends faithfully.

> [!NOTE]
> This does **not** call for `makedb --forcecca`. For eukaryotes the CCA tail is added to the
> reference by default already — `--forcecca` only overrides the prokaryotic default, where CCA
> is usually genomically encoded and so is not appended. Passing it to a `-s euk` build changes
> nothing. Reach for it only when building a `bact`/`arch` database whose organism does not
> encode CCA genomically. How faithfully a method captures the 3′ end is a property of the
> library, not of how the reference is built.

---

## Reading the output

Whichever example you followed, the run leaves you with the same shape of output:

```text
<output>/
├── <output>.h5ad     # everything, in one AnnData object
├── results/          # the same numbers as text files
└── graphs/           # one subfolder per graph type
```

Where to look first, by question:

| Question                         | Where                                                                                             |
| -------------------------------- | ------------------------------------------------------------------------------------------------- |
| Did mapping work?                | `results/<exp>-mapinfo.txt`; the overall alignment rate printed by `map`                          |
| What are my reads?               | `graphs/count/` — tRNA versus everything else, per sample                                         |
| Do replicates agree?             | `graphs/pca/` and `graphs/correlation/` — check before trusting anything downstream               |
| What changed between groups?     | `graphs/volcano/` and `graphs/heatmap/`                                                           |
| Where on the tRNA are the reads? | `graphs/coverage/` — the stacked overview at the top shows how much signal is transcript-specific |
| Are modifications visible?       | `graphs/mismatch/` — the payoff for demethylase protocols                                         |
| Which tRNAs behave alike?        | `graphs/cluster/` — requires `analyze cluster` first                                              |

By default every plot is drawn from **unique** reads — reads assigned to exactly one tRNA
transcript. `--allreads` switches the whole command to all reads at once, deliberately, so two
plots of the same data can never rest on different denominators.

Everything in `results/` is derived from the `.h5ad`, so for anything the plots do not answer,
load the object directly — see [Advanced Usage](advanced_usage.md) for the Python API and
[Data Structure](data_structure.md) for the full schema.

## Where to go next

- [CLI Reference](cli_reference.md) — every command and flag
- [Data Structure](data_structure.md) — what is in the AnnData object and where
- [Advanced Usage](advanced_usage.md) — config files, style files, custom analysis
- [Improvements from tRAX](improvements_from_trax.md) — if you are coming from tRAX
