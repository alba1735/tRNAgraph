"""Every graph type plots unique (transcript-specific) counts unless --allreads says otherwise.

The defect this guards against is *silent* inconsistency: before this, `plotsCluster` and
`plotsSeqlogo` were hardcoded to unique counts while `plotsCompare` was hardcoded to all
reads, `plotsRadar` and `plotsCorrelation` emitted both every run, and `--pcareadtypes`
defaulted to a mixture -- so two plots of one dataset could rest on different denominators
with nothing on either figure saying which.

"Unique" means TRANSCRIPT-SPECIFIC (bowtie2's YR tag), not genome-level uniqueness --
`toolsCountReads.getbamcounts()` gates `adduniquecount()` on `isuniquetrnamapping()`, and
`toolsGetCoverage` fills 'uniquecoverage' from the same predicate. The genome MAPQ >= 2
filter is separate and always on (see test_filtermultimapped_default.py).

The source scan below is the part that matters for the future: a runtime test can only
check the graph types that exist today, whereas the scan fails when a *newly added* module
hardcodes a basis. Anything that genuinely cannot follow the flag must be added to
EXEMPT_SITES with a reason, so an exemption is always a deliberate, reviewed act.
"""
import ast
import inspect
import re
from pathlib import Path

import pytest

from trnagraph.modules import toolsTG

MODULES_DIR = Path(inspect.getfile(toolsTG)).parent

# Modules whose basis is fixed by the data rather than by choice, or which compare bases
# deliberately. Each entry must say why.
EXEMPT_SITES = {
    # adata.uns's amino/anticodon/type/nontRNA count tables are read straight from the
    # tRAX-style count files and have no unique-only counterpart at all, so count plots
    # cannot follow the flag without new computation in adataBuild.
    "plotsCount.py": "uns count tables have no unique-only counterpart",
    # These two pin a both-bases comparison on purpose; see OVERVIEW_TRNA_READTYPES.
    "plotsPca.py": "OVERVIEW_TRNA_READTYPES is a deliberate labelled comparison",
    "plotsVolcano.py": "OVERVIEW_TRNA_READTYPES is a deliberate labelled comparison",
    # Reads resolved column names in and only filters the two totals out of the stacked
    # combine; it never selects a basis.
    "plotsHeatmap.py": "filters the totals out of the combined heatmap, does not select a basis",
}

# A literal obs column naming a read basis, e.g. 'nreads_total_unique_norm'. Matched in
# full, not searched for, so that prose mentioning a column (in a warning message) and the
# neighbouring all-feature columns are not caught:
#
#   'nreads_total_norm_allfeatures' -- the all-feature-controlled normalization used by the
#       non-tRNA and combined panels, which have no unique basis at all because uns's
#       non-tRNA counts are not broken down by transcript specificity.
#   'nreads_total_raw' -- raw counts, the input those panels renormalize. Raw-versus-
#       normalized is --variant's axis, not --allreads'.
BASIS_COLUMN = re.compile(r"nreads_[a-z]+(_unique)?_norm")


def _plot_modules():
    return sorted(
        p for p in MODULES_DIR.glob("plots*.py")
        # plotsLegacy*.py are orphaned tRAX ports wired into nothing (see roadmap.md's
        # "Legacy script cleanup"); they are dead code, not a graph type anyone can run.
        if not p.name.startswith("plotsLegacy")
    )


def _string_literals(path):
    """Every string literal in the file, excluding comments and docstrings.

    Comments are excluded because they legitimately discuss the column names, and
    docstrings because they document them -- neither can change what gets plotted.
    """
    tree = ast.parse(path.read_text())
    docstrings = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            doc = ast.get_docstring(node, clean=False)
            if doc is not None:
                docstrings.add(doc)
    return [
        node.value for node in ast.walk(tree)
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
        and node.value not in docstrings
    ]


@pytest.mark.parametrize("path", _plot_modules(), ids=lambda p: p.name)
def test_plot_module_does_not_hardcode_a_read_basis(path):
    if path.name in EXEMPT_SITES:
        pytest.skip(f"{path.name}: {EXEMPT_SITES[path.name]}")
    offenders = sorted({s for s in _string_literals(path) if BASIS_COLUMN.fullmatch(s)})
    assert not offenders, (
        f"{path.name} names read-basis obs columns directly: {offenders}. Resolve them "
        f"through toolsTG.resolve_readtype() so --allreads switches every graph type at "
        f"once, or add {path.name} to EXEMPT_SITES with a reason."
    )


def test_every_exempt_site_still_exists():
    """An exemption for a module that no longer exists is stale permission."""
    names = {p.name for p in _plot_modules()}
    assert not (set(EXEMPT_SITES) - names), f"stale exemptions: {set(EXEMPT_SITES) - names}"


class TestResolveReadtype:
    def test_unique_is_the_default_basis(self):
        assert toolsTG.read_basis(False) == toolsTG.READ_BASIS_UNIQUE
        assert toolsTG.read_basis(True) == toolsTG.READ_BASIS_ALL

    @pytest.mark.parametrize("readtype", toolsTG.DUAL_BASIS_READTYPES)
    def test_both_bases_resolve_for_every_dual_basis_readtype(self, readtype):
        assert toolsTG.resolve_readtype(readtype, toolsTG.READ_BASIS_UNIQUE) == f"nreads_{readtype}_unique_norm"
        assert toolsTG.resolve_readtype(readtype, toolsTG.READ_BASIS_ALL) == f"nreads_{readtype}_norm"

    @pytest.mark.parametrize("readtype", ["total_unique", "wholecounts_unique"])
    def test_a_readtype_naming_a_basis_is_rejected(self, readtype):
        """--diffrts/--pcareadtypes no longer carry the basis; accepting it would let one
        graph type sit on a different denominator than the rest."""
        with pytest.raises(ValueError, match="names a read basis"):
            toolsTG.resolve_readtype(readtype, toolsTG.READ_BASIS_UNIQUE)

    def test_a_full_column_name_is_rejected(self):
        with pytest.raises(ValueError, match="already a full obs column name"):
            toolsTG.resolve_readtype("nreads_total_norm", toolsTG.READ_BASIS_UNIQUE)

    @pytest.mark.parametrize("readtype", toolsTG.ALL_READS_ONLY_READTYPES)
    def test_all_reads_only_readtypes_fall_back_rather_than_fail(self, readtype, caplog):
        """adataBuild never writes a unique column for the pre-tRNA/antisense categories, so
        asking for one under the unique basis must degrade to all reads, loudly."""
        assert toolsTG.resolve_readtype(readtype, toolsTG.READ_BASIS_UNIQUE) == f"nreads_{readtype}_norm"
        assert "no transcript-specific" in caplog.text

    def test_an_unknown_basis_is_rejected(self):
        with pytest.raises(ValueError, match="Unknown read basis"):
            toolsTG.resolve_readtype("total", "sometimes")


class TestResolveCovtype:
    def test_default_follows_the_basis(self):
        assert toolsTG.resolve_covtype(None, toolsTG.READ_BASIS_UNIQUE) == "uniquecoverage"
        assert toolsTG.resolve_covtype(None, toolsTG.READ_BASIS_ALL) == "coverage"

    @pytest.mark.parametrize("basis", [toolsTG.READ_BASIS_UNIQUE, toolsTG.READ_BASIS_ALL])
    def test_an_explicit_covtype_is_always_honoured(self, basis):
        """tRAX's four coverage categories are a partition a user may legitimately inspect
        in either mode, so the basis only supplies the default."""
        assert toolsTG.resolve_covtype("mismatchedbases", basis) == "mismatchedbases"


def test_cli_graph_exposes_allreads_and_bare_readtypes():
    from trnagraph import cli
    params = inspect.signature(cli.graph).parameters
    assert "allreads" in params, "--allreads is the single command-wide basis switch"
    for name in ("diffrts", "pcareadtypes"):
        for value in params[name].default.default:
            assert "unique" not in value, (
                f"--{name} default {value!r} names a read basis; the basis belongs to "
                f"--allreads alone so the two cannot disagree."
            )
    assert params["covtype"].default.default is None, (
        "--covtype must default to None so resolve_covtype() can pick per basis while still "
        "telling an explicit value apart from an unset one."
    )
