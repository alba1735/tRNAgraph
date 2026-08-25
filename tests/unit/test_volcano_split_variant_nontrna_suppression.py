"""End-to-end regression test for roadmap.md's small-RNA/non-tRNA split-variant bug, at the
plotsVolcano.visualizer() seam named in the original report. Even when adata.uns['nontRNA_counts']
is present and non-empty, a split (--variant norm:o60/u60) render must never generate non-tRNA or
combined volcano plots -- only the complete (non-split) variant may, since non-tRNA reads can't be
meaningfully partitioned by a tRNA length cutoff."""
import glob
import os

import anndata as ad
import numpy as np
import pandas as pd

from trnagraph.modules.plotsVolcano import visualizer

# plotsVolcano.visualizer() takes fully-resolved obs column names: adataGraph applies the
# command-wide read basis once, via toolsTG.resolve_readtype(), before calling in.
READTYPES = ['nreads_total_unique_norm']
CUTOFF = 5


def _make_adata():
    # A single differential tRNA alone destabilizes PyDESeq2's dispersion-trend fit (it assumes
    # most features are non-differential) -- a handful of flat tRNAs alongside it, mirroring
    # test_log2fc.py's pattern, keeps the fit well-behaved.
    rng = np.random.default_rng(0)
    groups = ['A', 'B']
    samples = [f'{g}_rep{i}' for g in groups for i in range(3)]
    sample_group = {s: s.rsplit('_rep', 1)[0] for s in samples}
    trna_means = {'trnaHigh': {'A': 800, 'B': 50}}
    trna_means.update({f'trnaFlat{i}': {'A': 400, 'B': 400} for i in range(6)})

    rows = []
    for trna, group_means in trna_means.items():
        for sample in samples:
            group = sample_group[sample]
            mean = group_means[group]
            raw_tu = rng.negative_binomial(n=10, p=10 / (10 + mean))
            raw_t = rng.negative_binomial(n=10, p=10 / (10 + mean))
            rows.append({
                'trna': trna, 'sample': sample, 'group': group,
                'nreads_total_unique_raw': raw_tu, 'nreads_total_unique_norm': raw_tu,
                'nreads_total_raw': raw_t, 'nreads_total_norm': raw_t,
            })
    obs = pd.DataFrame(rows)
    obs.index = [f"{r.trna}_{r['sample']}" for _, r in obs.iterrows()]
    adata = ad.AnnData(X=np.zeros((len(obs), 1)), obs=obs)

    # Non-tRNA counts: present and well above CUTOFF for every sample/group -- if suppression
    # weren't enforced for split variants, this would be more than enough to produce output.
    adata.uns['nontRNA_counts'] = pd.DataFrame(
        {s: [500 + rng.integers(-20, 20)] for s in samples}, index=['ncRNA1']
    )
    # Required for the combined tRNA+non-tRNA ("allRNA") plot specifically -- an unrelated gate
    # from the non-tRNA-only one, needed here only so the full-variant case actually reaches it.
    adata.uns['deseq2_sizefactors_allfeatures'] = {s: 1.0 for s in samples}
    return adata


def _output_files(output_dir):
    return set(os.path.basename(p) for p in glob.glob(os.path.join(output_dir, '**', '*.pdf'), recursive=True))


def test_full_variant_generates_nontrna_and_allrna_volcano_output(tmp_path):
    adata = _make_adata()
    output = str(tmp_path / 'full') + os.sep

    visualizer(adata, 'group', READTYPES, CUTOFF, output, threaded=False, is_full_variant=True)

    files = _output_files(output)
    assert any('nontRNA' in f for f in files), f"expected a nontRNA volcano file, got: {files}"
    assert any('allRNA' in f for f in files), f"expected a combined tRNA+non-tRNA (allRNA) volcano file, got: {files}"


def test_split_variant_never_generates_nontrna_or_allrna_volcano_output_even_with_data_present(tmp_path):
    """The top-level `{grp}_combined_{cutoff}_volcano.pdf` overview always exists regardless (it
    just has 2 tRNA-only panels instead of 4) -- what must never appear for a split variant is the
    non-tRNA-only ('nontRNA') and combined tRNA+non-tRNA ('allRNA') individual plots."""
    adata = _make_adata()
    output = str(tmp_path / 'split') + os.sep

    visualizer(adata, 'group', READTYPES, CUTOFF, output, threaded=False, is_full_variant=False)

    files = _output_files(output)
    assert not any('nontRNA' in f for f in files), f"non-tRNA volcano file must not be generated for a split variant, got: {files}"
    assert not any('allRNA' in f for f in files), f"combined tRNA+non-tRNA volcano file must not be generated for a split variant, got: {files}"
    assert any('tRNA_' in f for f in files), "the ordinary tRNA-only volcano plot should still be generated"
