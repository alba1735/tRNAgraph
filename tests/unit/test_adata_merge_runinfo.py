"""Regression tests for anndataMerger.merge()'s conflicting-run-info enforcement (roadmap.md
Phase 2: "Merge Logic" -- mirrors `analyze addsplit`'s --force-gated database/gtf provenance
check, previously only implemented for addsplit and missing entirely from `tools merge`)."""
import logging
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from trnagraph.modules.adataMerge import anndataMerger


def _make_adata(database='db1', gtf='genes1.gtf', n_obs=2, index_prefix='obs'):
    obs = pd.DataFrame({'sample': [f's{i}' for i in range(n_obs)]}, index=[f'{index_prefix}{i}' for i in range(n_obs)])
    adata = ad.AnnData(X=np.zeros((n_obs, 1)), obs=obs)
    adata.uns['trnagraphruninfo'] = {'flags': {'database': database, 'gtf': gtf}}
    return adata


def _make_merger(adata1, adata2, force=False, output='/tmp/out.h5ad'):
    merger = anndataMerger.__new__(anndataMerger)
    merger.logger = logging.getLogger("trnagraph.modules.adataMerge")
    merger.adata1 = adata1
    merger.adata2 = adata2
    merger.args = SimpleNamespace(dropno=False, droprna=False, output=output, force=force)
    return merger


def test_merge_raises_on_conflicting_database_without_force():
    adata1 = _make_adata(database='db1', index_prefix='obj1_')
    adata2 = _make_adata(database='db2', index_prefix='obj2_')
    merger = _make_merger(adata1, adata2, force=False)

    with pytest.raises(ValueError, match='database'):
        merger.merge()


def test_merge_raises_on_conflicting_gtf_without_force():
    adata1 = _make_adata(gtf='genesA.gtf', index_prefix='obj1_')
    adata2 = _make_adata(gtf='genesB.gtf', index_prefix='obj2_')
    merger = _make_merger(adata1, adata2, force=False)

    with pytest.raises(ValueError, match='gtf'):
        merger.merge()


def test_merge_proceeds_past_conflict_check_with_force(caplog):
    adata1 = _make_adata(database='db1', index_prefix='obj1_')
    adata2 = _make_adata(database='db2', index_prefix='obj2_')
    merger = _make_merger(adata1, adata2, force=True)

    # --force should bypass the conflict check (logged as a warning) and let execution reach the
    # next stage of merge() -- which then fails on missing uns['amino_counts'] in this minimal
    # fixture, proving the conflict check itself did not block it.
    with pytest.raises(KeyError, match='amino_counts'), \
         caplog.at_level(logging.WARNING, logger="trnagraph.modules.adataMerge"):
        merger.merge()
    warning_messages = "\n".join(r.message for r in caplog.records if r.levelno == logging.WARNING)
    assert 'database' in warning_messages


def test_merge_does_not_raise_when_provenance_matches():
    adata1 = _make_adata(database='db1', gtf='genes1.gtf', index_prefix='obj1_')
    adata2 = _make_adata(database='db1', gtf='genes1.gtf', index_prefix='obj2_')
    merger = _make_merger(adata1, adata2, force=False)

    # Matching provenance should not raise from the conflict check -- next failure (missing
    # uns['amino_counts']) proves it passed through cleanly.
    with pytest.raises(KeyError, match='amino_counts'):
        merger.merge()
