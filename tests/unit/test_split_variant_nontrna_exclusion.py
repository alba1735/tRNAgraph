"""Regression tests for roadmap.md's "Non-tRNA counts should be excluded from read-length split
variants, not just their plots".

An earlier fix stopped a split variant recomputing its own non-tRNA counts (it reused the full
variant's instead) and suppressed the non-tRNA plots, but the underlying counts were still
carried and the on-disk split outputs were never touched at all. On a real hg38 split build the
`other` biotype came out 12.41M in u60 against 0.99M in o60 -- a 92/8 partition of a feature
class the cutoff has no business partitioning, since a threshold chosen to separate tRNA
fragments from full-length tRNAs is meaningless for features differing in length by orders of
magnitude. Worse, `readcounts.txt` is what the split's DESeq2 run consumes, so its all-feature
size factors were estimated over an arbitrarily truncated non-tRNA pool.

Note `--gtf` is optional: with no GTF no non-tRNA feature is ever counted, so every filter here
must be a well-behaved no-op rather than an error. That is also why the bacterial `vibrChol1`
suite cannot exercise this fix -- its genetypes output is tRNA-only -- and these unit tests
carry the weight.
"""
import pandas as pd

from trnagraph.modules import toolsTG


# Feature names and their genetypes.txt type column, in the vocabulary a real hg38 build
# produces. The type column is the authoritative per-feature annotation; the name-substring
# predicate below must agree with it.
GENETYPES_SAMPLE = [
    ('tRNA-Ala-AGC-1_wholecounts',      'trna_wholecounts'),
    ('tRNA-Ala-AGC-1_fiveprime',        'trna_fiveprime'),
    ('tRNA-Ala-AGC-1-1_wholeprecounts', 'trna_wholeprecounts'),
    ('tRNA-Ala-AGC-1-1',                'tRNA_locus'),
    ('tRNA-Ala-AGC-1',                  'tRNA'),
    ('tRX-Ala-AGC-1_wholecounts',       'trna_wholecounts'),
    ('ENSG00000201098',                 'snRNA'),
    ('ENSG00000206652',                 'Y_RNA'),
    ('ENSG00000283272',                 'miRNA'),
    ('ENSG00000278267',                 'snoRNA_pseudogene'),
    ('ENSG00000275757',                 'rRNA'),
    ('ENSG00000230876',                 'lncRNA'),
]

TRNA_GENETYPES = {'tRNA', 'tRNA_locus'}


def _is_trna_by_genetype(type_value):
    """The authority: genetypes.txt's own type column."""
    return type_value.startswith('trna_') or type_value in TRNA_GENETYPES


def test_trna_feature_predicate_agrees_with_the_genetypes_type_column():
    """The tRNA/non-tRNA split is decided by a name-substring rule that the DESeq2 control-gene
    selection already used. Pin that it agrees with genetypes.txt's independent type annotation,
    so the two definitions cannot drift apart unnoticed. Verified as 0 disagreements across all
    8343 features of a real hg38 build."""
    for name, type_value in GENETYPES_SAMPLE:
        assert toolsTG.is_trna_feature(name) == _is_trna_by_genetype(type_value), name


def test_trna_feature_predicate_covers_both_trna_and_trx():
    assert toolsTG.is_trna_feature('tRNA-Ala-AGC-1_wholecounts')
    assert toolsTG.is_trna_feature('tRX-Ala-AGC-1_wholecounts')
    assert not toolsTG.is_trna_feature('ENSG00000201098')


import os
from types import SimpleNamespace

from trnagraph.modules.adataBuild import _filter_nontrna_rows_from_counts_file


def _write_counts(path, feature_names):
    """A readcounts.txt-shaped file: features as rows, samples as columns."""
    df = pd.DataFrame(
        {'s1': range(len(feature_names)), 's2': range(len(feature_names))},
        index=feature_names,
    )
    df.to_csv(path, sep='\t')
    return df


def test_filter_drops_non_trna_rows_and_keeps_every_trna_row(tmp_path):
    path = str(tmp_path / 'exp-readcounts.txt')
    names = [n for n, _ in GENETYPES_SAMPLE]
    _write_counts(path, names)

    removed = _filter_nontrna_rows_from_counts_file(path)

    kept = pd.read_csv(path, sep='\t', index_col=0)
    assert list(kept.index) == [n for n in names if toolsTG.is_trna_feature(n)]
    assert removed == 6  # the six ENSG rows


def test_filter_is_a_noop_when_there_is_no_gtf_and_so_no_non_trna_rows(tmp_path):
    """--gtf is optional; with no GTF nothing non-tRNA is ever counted, so the filter must leave
    the file untouched rather than erroring or rewriting it needlessly."""
    path = str(tmp_path / 'exp-readcounts.txt')
    names = [n for n, t in GENETYPES_SAMPLE if _is_trna_by_genetype(t)]
    original = _write_counts(path, names)
    before = pd.read_csv(path, sep='\t', index_col=0)

    removed = _filter_nontrna_rows_from_counts_file(path)

    assert removed == 0
    pd.testing.assert_frame_equal(pd.read_csv(path, sep='\t', index_col=0), before)
    assert list(before.index) == list(original.index)


def test_filter_tolerates_a_missing_file(tmp_path):
    """Not every count file exists for every build configuration."""
    assert _filter_nontrna_rows_from_counts_file(str(tmp_path / 'absent.txt')) == 0


from trnagraph.modules.adataBuild import _filter_nontrna_rows_from_type_counts_file

# The full row vocabulary a real hg38 typecounts.txt carries, in file order.
TYPECOUNTS_ROWS = [
    'other', 'rRNA', 'Mt_tRNA', 'Mt_rRNA', 'misc_RNA', 'Y_RNA', 'rRNA_pseudogene',
    'snRNA_pseudogene', 'snoRNA_pseudogene', 'vault_RNA', 'lncRNA', 'misc_RNA_pseudogene',
    'miRNA', 'sRNA', 'scaRNA', 'snRNA', 'snoRNA',
    'pretRNA_antisense', 'pretRNA', 'tRNA_antisense', 'tRNA',
]


def test_type_counts_filter_keeps_only_the_trna_derived_labels(tmp_path):
    """typecounts.txt/typerealcounts.txt are indexed by TYPE label, not feature name. Its
    tRNA-derived rows are the four literals toolsCountReads writes directly, plus Mt_tRNA; every
    other row is a GTF biotype from emblbiotypes."""
    path = str(tmp_path / 'exp-typecounts.txt')
    _write_counts(path, TYPECOUNTS_ROWS)

    _filter_nontrna_rows_from_type_counts_file(path)

    assert list(pd.read_csv(path, sep='\t', index_col=0).index) == [
        'Mt_tRNA', 'pretRNA_antisense', 'pretRNA', 'tRNA_antisense', 'tRNA'
    ]


def test_type_counts_filter_keeps_mt_trna(tmp_path):
    """Mitochondrial tRNAs are tRNAs and are kept in split variants, even though the Mt_tRNA row
    is GTF-driven rather than database-derived (gtRNAdb and tRNAscan-SE exclude mitochondrial
    tRNAs, so a makedb-built database contains none and this row is structurally 0.0 today).
    They are planned to become first-class -- see roadmap.md's mitochondrial tRNA handling item
    -- so the data keeps flowing rather than being dropped now and re-added later. At roughly
    60-75nt a read-length cutoff also partitions them meaningfully, unlike the non-tRNA feature
    classes this filter targets."""
    path = str(tmp_path / 'exp-typecounts.txt')
    _write_counts(path, ['Mt_tRNA', 'Mt_rRNA', 'tRNA', 'rRNA'])

    _filter_nontrna_rows_from_type_counts_file(path)

    assert list(pd.read_csv(path, sep='\t', index_col=0).index) == ['Mt_tRNA', 'tRNA']


from unittest.mock import patch

from trnagraph.modules.adataBuild import AnalysisPipeline


def _pipeline(split_tag):
    """An AnalysisPipeline with only the attributes run_deseq2()/the filters read."""
    pipeline = AnalysisPipeline.__new__(AnalysisPipeline)
    pipeline.logger = __import__('logging').getLogger('trnagraph.modules.adataBuild')
    pipeline.split_tag = split_tag
    pipeline.expinfo = SimpleNamespace(
        genecounts='rc.txt', normalizedcounts='nrc.txt', sizefactors='sf.txt',
        normalizedcounts_allfeatures='af_nrc.txt', allfeaturesizefactors='af_sf.txt',
        trnacounts='tc.txt', normalizedtrnacounts='ntc.txt', trnasizefactors='tsf.txt',
        resultsdir='results', expname='exp', genetypes='gt.txt',
        genetypecounts='typec.txt', genetyperealcounts='typerc.txt',
    )
    return pipeline


def _deseq_prefixes(pipeline):
    calls = []
    with patch.object(AnalysisPipeline, 'run_deseq2_on_file',
                      lambda self, **kw: calls.append(kw['prefix'])):
        pipeline.run_deseq2()
    return calls


def test_full_variant_still_runs_the_all_feature_deseq2_pass():
    assert _deseq_prefixes(_pipeline(split_tag=None)) == ['', 'allfeature_', 'trna_']


def test_split_variant_skips_the_all_feature_deseq2_pass():
    """With non-tRNA features removed from a split's readcounts.txt, "all features" and "tRNA
    only" are the same set, so an allfeature_ pass would duplicate the primary one. Computing it
    from the UNfiltered counts instead -- the previous behaviour -- estimated size factors over a
    non-tRNA pool the cutoff had arbitrarily truncated, and those size factors were reachable by
    users as `--variant allfeatures:<tag>`."""
    assert _deseq_prefixes(_pipeline(split_tag='u60')) == ['', 'trna_']


def test_split_variant_filters_all_four_count_outputs(tmp_path):
    """readcounts/genetypes are feature-indexed; typecounts/typerealcounts are type-indexed.
    All four carry non-tRNA rows for a split today."""
    pipeline = _pipeline(split_tag='u60')
    feature_rows = [n for n, _ in GENETYPES_SAMPLE]
    for attr, rows in [('genecounts', feature_rows), ('genetypes', feature_rows),
                       ('genetypecounts', TYPECOUNTS_ROWS), ('genetyperealcounts', TYPECOUNTS_ROWS)]:
        path = str(tmp_path / f'{attr}.txt')
        _write_counts(path, rows)
        setattr(pipeline.expinfo, attr, path)

    pipeline.apply_split_variant_filters()

    for attr in ('genecounts', 'genetypes'):
        idx = list(pd.read_csv(getattr(pipeline.expinfo, attr), sep='\t', index_col=0).index)
        assert all(toolsTG.is_trna_feature(n) for n in idx), attr
        assert idx, attr
    for attr in ('genetypecounts', 'genetyperealcounts'):
        idx = list(pd.read_csv(getattr(pipeline.expinfo, attr), sep='\t', index_col=0).index)
        assert idx == ['Mt_tRNA', 'pretRNA_antisense', 'pretRNA', 'tRNA_antisense', 'tRNA'], attr


def test_full_variant_filters_nothing(tmp_path):
    """The complete variant keeps every non-tRNA feature -- that is the whole point of --gtf."""
    pipeline = _pipeline(split_tag=None)
    path = str(tmp_path / 'rc.txt')
    _write_counts(path, [n for n, _ in GENETYPES_SAMPLE])
    pipeline.expinfo.genecounts = path
    pipeline.expinfo.genetypes = pipeline.expinfo.genetypecounts = pipeline.expinfo.genetyperealcounts = str(tmp_path / 'absent.txt')

    pipeline.apply_split_variant_filters()

    assert len(pd.read_csv(path, sep='\t', index_col=0)) == len(GENETYPES_SAMPLE)


import numpy as np

from trnagraph.modules.adataBuild import _empty_split_nontrna_counts


def _adata_with_samples(samples, nontrna_columns=None):
    import anndata as ad
    obs = pd.DataFrame({'sample': samples}, index=[f'g_{s}' for s in samples])
    adata = ad.AnnData(X=np.zeros((len(samples), 2)), obs=obs)
    if nontrna_columns is not None:
        adata.uns['nontRNA_counts'] = pd.DataFrame(columns=nontrna_columns)
    return adata


def test_split_nontrna_counts_is_empty_but_keeps_the_sample_columns():
    """A split's non-tRNA counts must be empty, not a copy of the full variant's -- but a
    well-formed empty frame carrying the sample axis, so anything inspecting the object's shape
    sees an empty result rather than something it has to special-case."""
    adata = _adata_with_samples(['s1', 's2'], nontrna_columns=['s1', 's2'])

    empty = _empty_split_nontrna_counts(adata)

    assert list(empty.columns) == ['s1', 's2']
    assert len(empty) == 0


def test_split_nontrna_counts_falls_back_to_the_objects_samples_with_no_gtf():
    """--gtf is optional. With no GTF the full variant's own non-tRNA frame is already empty and
    may carry no columns at all, so the sample axis has to come from the object itself."""
    adata = _adata_with_samples(['s1', 's2'], nontrna_columns=[])

    empty = _empty_split_nontrna_counts(adata)

    assert list(empty.columns) == ['s1', 's2']
    assert len(empty) == 0


def test_split_nontrna_counts_handles_a_missing_full_variant_frame():
    adata = _adata_with_samples(['s1'])

    assert list(_empty_split_nontrna_counts(adata).columns) == ['s1']


import pytest


def test_allfeatures_on_a_split_fails_saying_it_is_deliberate():
    """The layer is absent by design, not by accident. The generic "was this normalization
    computed for this split?" message invites the reader to go looking for a build flag that
    would produce it -- there isn't one."""
    import anndata as ad
    adata = _adata_with_samples(['s1', 's2'])
    adata.layers['norm'] = np.zeros((2, 2))
    adata.layers['norm_u60'] = np.zeros((2, 2))
    adata.uns['size_splits'] = {'u60': {}}

    spec = toolsTG.parse_variant(adata, 'allfeatures:u60')
    with pytest.raises(ValueError) as excinfo:
        toolsTG.build_variant_view(adata, spec)

    message = str(excinfo.value)
    assert 'complete' in message and 'split' in message
    assert 'was this normalization computed' not in message


def test_allfeatures_on_the_full_variant_is_unaffected():
    import anndata as ad
    adata = _adata_with_samples(['s1', 's2'])
    adata.layers['norm_allfeatures'] = np.zeros((2, 2))

    spec = toolsTG.parse_variant(adata, 'allfeatures:full')
    view = toolsTG.build_variant_view(adata, spec)

    assert view.X.shape == (2, 2)


from trnagraph.modules.toolsSchemas import VariantContribution


def test_variant_contribution_accepts_no_all_feature_data():
    """A split no longer writes an allfeature SizeFactors file, so the loader that reads one back
    must tolerate its absence and the contribution must be constructible without it -- otherwise
    dropping the allfeature pass would break the split build outright."""
    frame = pd.DataFrame({'a': [1.0]})
    contribution = VariantContribution(
        x_raw=frame, x_norm=frame, obsm_counts=frame,
        sizefactors_trna={'s1': 1.0},
        type_counts=frame, type_real_counts=frame,
        amino_counts=frame, anticodon_counts=frame, nontrna_counts=frame,
    )

    assert contribution.x_norm_allfeatures is None
    assert contribution.sizefactors_allfeatures is None


import ast
import inspect

from trnagraph.modules import adataBuild


def test_every_split_pipeline_construction_passes_split_tag():
    """There are two entry points that build a split variant -- `analyze build
    --readlengthsplit` and `analyze addsplit` -- and they construct AnalysisPipeline
    independently. A split built without split_tag silently keeps its non-tRNA rows and its
    all-feature pass, which is the whole bug. Assert structurally that no construction inside a
    split loop can omit it, rather than trusting both call sites to be kept in step by hand."""
    source = inspect.getsource(adataBuild)
    tree = ast.parse(source)

    constructions = [
        node for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == 'AnalysisPipeline'
    ]
    assert len(constructions) >= 3, "expected the full-variant plus both split call sites"

    # A construction is split-building iff bamdir/results_dir_name were derived from a tag; the
    # observable proxy is that it passes split_tag. Exactly one construction (the full variant's)
    # may omit it.
    without_tag = [
        c for c in constructions
        if not any(kw.arg == 'split_tag' for kw in c.keywords)
    ]
    assert len(without_tag) == 1, (
        f"{len(without_tag)} AnalysisPipeline constructions omit split_tag; only the "
        f"full/complete variant's may."
    )


from trnagraph.modules.adataBuild import _load_split_nontrna_read_counts


def test_nontrna_read_counts_load_is_empty_when_the_allfeature_file_is_absent(tmp_path):
    """Regression for a crash on the first real split build after all-feature normalization
    became complete-variant only. AnnDataBuilder.__init__ derived uns['nontRNA_counts'] by
    reading the allfeature normalizedreadcounts file unconditionally -- a second all-feature read
    beyond the size-factors one -- so the split loader died with FileNotFoundError on a file the
    split no longer writes. For a split the answer is an empty frame anyway: non-tRNA features
    are excluded from splits entirely."""
    absent = str(tmp_path / 'allfeature' / 'exp-allfeature_normalizedreadcounts.txt')
    fallback_columns = ['s1', 's2']

    result = _load_split_nontrna_read_counts(absent, fallback_columns, lambda df: df)

    assert list(result.columns) == fallback_columns
    assert len(result) == 0


def test_nontrna_read_counts_load_drops_trna_and_trx_rows_when_the_file_exists(tmp_path):
    """The complete variant's path is unchanged: read the all-feature-normalized counts and keep
    only the non-tRNA rows."""
    path = str(tmp_path / 'exp-allfeature_normalizedreadcounts.txt')
    pd.DataFrame(
        {'s1': [1.0, 2.0, 3.0, 4.0], 's2': [1.0, 2.0, 3.0, 4.0]},
        index=['tRNA-Ala-AGC-1_wholecounts', 'tRX-Ala-AGC-1_wholecounts',
               'ENSG00000201098', 'ENSG00000206652'],
    ).to_csv(path, sep='\t')

    result = _load_split_nontrna_read_counts(
        path, ['s1', 's2'], lambda df: df.set_index(df.columns[0]) if df.index.dtype != 'object' else df
    )

    assert list(result.index) == ['ENSG00000201098', 'ENSG00000206652']


def test_no_all_feature_file_is_read_without_an_existence_guard():
    """Split variants write no allfeature/ outputs, so every read of one has to tolerate its
    absence. Two such reads existed and only one was guarded on the first attempt -- the
    unguarded second crashed the first real split build after the change. Assert the invariant
    structurally rather than relying on having grepped for every attribute name."""
    source = inspect.getsource(adataBuild)
    tree = ast.parse(source)

    offenders = []
    for func in [n for n in ast.walk(tree) if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]:
        reads = [
            n for n in ast.walk(func)
            if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
            and n.func.attr in ('read_csv', 'read_table') and n.args
            and 'allfeature' in ast.unparse(n.args[0]).lower()
        ]
        if not reads:
            continue
        has_guard = any(
            isinstance(n, ast.Attribute) and n.attr == 'exists'
            for n in ast.walk(func)
        )
        if not has_guard:
            offenders.extend(f"{func.name}:{r.lineno}" for r in reads)

    assert not offenders, (
        f"all-feature file read with no os.path.exists guard in: {offenders}. "
        f"A read-length split variant writes no allfeature/ outputs."
    )
