"""End to end for one agreement run, with the DESeq2 boundary stubbed.

The fits themselves are toolsTG's and are tested there; what this pins is the orchestration --
one figure per reference-anchored contrast, a combined overview beside them, a table per
contrast, and the membership stored on the object so the figure and the file cannot disagree.
"""
import matplotlib
matplotlib.use('Agg')
import anndata as ad
import numpy as np
import pandas as pd
import pytest

from trnagraph.modules import plotsAgreement
from trnagraph.modules.toolsSchemas import MultivariateConfig


REF = 'Day 0'
LEVELS = ['Day 0', 'Day 35', 'Day 70']


def _adata():
    obs = pd.DataFrame({
        'timepoint': pd.Categorical(LEVELS * 2, categories=LEVELS, ordered=True),
        'trna': [f'tRNA-{i}' for i in range(6)],
        'nreads_fiveprime_unique_norm': np.linspace(10, 100, 6),
    })
    adata = ad.AnnData(X=np.zeros((6, 2)), obs=obs)
    adata.uns['trnagraphruninfo'] = {'trnagraph_directory': None}
    return adata


def _frames(features=('tRNA-A', 'tRNA-B', 'tRNA-C')):
    data = {}
    for other, (lfc, p) in {
        'Day 35': ([4.0, -4.0, 0.1], [1e-8, 1e-8, 0.9]),
        'Day 70': ([5.0, -0.2, 0.2], [1e-9, 0.7, 0.8]),
    }.items():
        data[f'log2_{REF}-{other}'] = pd.Series(lfc, index=list(features))
        data[f'pval_{REF}-{other}'] = pd.Series(p, index=list(features))
    return pd.DataFrame(data, index=list(features))


@pytest.fixture
def stub_fits(monkeypatch):
    frame = _frames()
    monkeypatch.setattr(plotsAgreement, '_log2fc_frames',
                        lambda *a, **k: {'fiveprime': frame})
    return frame


def _run(tmp_path, adata=None, **kw):
    kw.setdefault('readtypes', ['nreads_fiveprime_unique_norm'])
    kw.setdefault('cutoff', 20)
    # adataGraph derives this once and hands it down; the module no longer builds it itself.
    kw.setdefault('results_dir', str(tmp_path).replace('graphs', 'results'))
    plotsAgreement.visualizer(adata if adata is not None else _adata(),
                              MultivariateConfig(grouping='timepoint'),
                              f'{tmp_path}/', **kw)
    return tmp_path


def test_one_figure_per_reference_anchored_contrast(stub_fits, tmp_path):
    _run(tmp_path)
    drawn = sorted(p.name for p in (tmp_path / 'individual').glob('*.pdf'))

    assert drawn == ['Day 0-Day 35_agreement.pdf', 'Day 0-Day 70_agreement.pdf']


def test_a_combined_overview_sits_beside_the_individual_folder(stub_fits, tmp_path):
    _run(tmp_path)

    assert (tmp_path / 'timepoint_combined_agreement.pdf').exists()


def test_the_membership_is_stored_on_the_object(stub_fits, tmp_path):
    adata = _adata()
    _run(tmp_path, adata=adata, config_name='basic')
    stored = adata.uns['multivariate']['basic']['agreement']

    assert sorted(stored) == ['Day 0-Day 35', 'Day 0-Day 70']


def test_the_stored_provenance_names_the_thresholds(stub_fits, tmp_path):
    adata = _adata()
    _run(tmp_path, adata=adata, config_name='basic')
    provenance = adata.uns['multivariate']['basic']['agreement']['Day 0-Day 70']['provenance']

    assert float(provenance['log2fc']) == 1.5
    assert float(provenance['padj']) == 0.001
    assert provenance['reference'] == REF


def test_a_table_is_written_per_contrast_beside_the_figures(stub_fits, tmp_path):
    """The tables used to be filed under the directory recorded in the object's build
    provenance and skipped outright when it was gone -- which it usually is, so they went
    missing exactly when they were most wanted. They now sit in the results/ twin of the
    figures' own path, which every run has by construction."""
    _run(tmp_path / 'graphs')

    written = sorted(p.name for p in (tmp_path / 'results').glob('*.tsv'))
    assert written == ['Day 0-Day 35_agreement.tsv', 'Day 0-Day 70_agreement.tsv']


def test_the_table_no_longer_depends_on_build_provenance(stub_fits, tmp_path):
    """An object whose build directory is gone -- the demo object records a scratch path, the
    hg38 object a server path -- still gets its tables."""
    adata = _adata()
    adata.uns['trnagraphruninfo'] = {'trnagraph_directory': str(tmp_path / 'gone')}
    _run(tmp_path / 'graphs', adata=adata)

    assert sorted(p.name for p in (tmp_path / 'results').glob('*.tsv')) == [
        'Day 0-Day 35_agreement.tsv', 'Day 0-Day 70_agreement.tsv']


def test_a_single_level_column_is_refused_by_name(stub_fits, tmp_path):
    adata = _adata()
    adata.obs['timepoint'] = pd.Categorical(['Day 0'] * 6, categories=['Day 0'])

    with pytest.raises(plotsAgreement.InvalidAgreementError):
        _run(tmp_path, adata=adata)


def test_every_file_written_is_reported_when_running_unthreaded(stub_fits, tmp_path, caplog):
    """The messages list was returned only in the threaded path and dropped otherwise, so a
    sequential run wrote figures and tables in silence -- caught in the demo, where the two
    multivariate types run outside the pool and agreement logged nothing between 'Generating'
    and 'generated!' while venn listed every file."""
    import logging

    adata = _adata()
    build = tmp_path / 'build'
    build.mkdir()
    adata.uns['trnagraphruninfo'] = {'trnagraph_directory': str(build)}
    with caplog.at_level(logging.INFO, logger='trnagraph.modules.plotsAgreement'):
        _run(tmp_path, adata=adata)

    logged = ' '.join(caplog.messages)
    assert 'Day 0-Day 70_agreement.pdf' in logged
    assert 'Day 0-Day 70_agreement.tsv' in logged
    assert 'combined' in logged


def test_the_threaded_path_returns_the_same_messages_instead(stub_fits, tmp_path):
    adata = _adata()
    returned = plotsAgreement.visualizer(
        adata, MultivariateConfig(grouping='timepoint'), f'{tmp_path}/',
        readtypes=['nreads_fiveprime_unique_norm'], cutoff=20, threaded=True)

    assert 'Day 0-Day 70_agreement.pdf' in returned
