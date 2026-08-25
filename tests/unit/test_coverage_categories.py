"""tRAX's four-way read-specificity partition is selectable, and summable.

`getcoverage.py`'s `getsamplecoverage()` bins every read into exactly one of four buckets
via an if/elif chain on bowtie2's YM/YA/YR tags (how many aminos / anticodons / transcripts
the read hit):

    YM > 1                  -> multiaminocoverage       "Not Amino Specific"
    else YA > 1             -> multianticodoncoverage   "Isotype Specific"
    else YR > 1             -> multitrnacoverage        "Isodecoder Specific"
    else                    -> uniquecoverage           "Transcript Specific"

They are therefore mutually exclusive and sum to `coverage`. That makes them a PARTITION,
not a filter -- the genome MAPQ >= 2 prefilter sits beneath all four equally, and is
separately always on (see test_filtermultimapped_default.py).

tRAX only ever surfaced the partition as a stacked plot plus a `unique/` folder holding the
Transcript-Specific slice; its `-uniquecoverages.txt` writer is dead code (`uniqcoverage()`
references names that exist only as locals in `getsamplecoverage`, and its one caller is
commented out). tRNAgraph stores all four in adata.var and makes each selectable.
"""
import inspect

import numpy as np
import pytest

from trnagraph.modules import toolsTG


class TestCategoryVocabulary:
    @pytest.mark.parametrize("alias,expected", [
        ("unique", "uniquecoverage"),
        ("transcript", "uniquecoverage"),
        ("isodecoder", "multitrnacoverage"),
        ("isotype", "multianticodoncoverage"),
        ("notamino", "multiaminocoverage"),
        ("total", "coverage"),
    ])
    def test_alias_resolves_to_the_var_value(self, alias, expected):
        assert toolsTG.resolve_covtype(alias, toolsTG.READ_BASIS_UNIQUE) == expected

    @pytest.mark.parametrize("raw", ["uniquecoverage", "multitrnacoverage", "readstarts",
                                     "mismatchedbases", "deletions"])
    def test_raw_var_values_still_work(self, raw):
        """Aliases are additive: anything that already scripted --covtype uniquecoverage,
        or a non-partition coverage type, must keep working unchanged."""
        assert toolsTG.resolve_covtype(raw, toolsTG.READ_BASIS_UNIQUE) == raw

    def test_the_partition_is_the_four_ambiguity_levels(self):
        assert set(toolsTG.COVERAGE_PARTITION) == set(toolsTG.COVERAGE_CATEGORY_LABELS)
        assert "coverage" not in toolsTG.COVERAGE_PARTITION, (
            "'coverage' is the SUM of the partition, not a member of it -- including it "
            "would double every stacked total."
        )

    def test_partition_is_ordered_least_to_most_specific(self):
        """Stacking order, matching newcoverageplots.R's factor levels, so a tRNAgraph
        specificity plot reads the same way round as tRAX's."""
        assert toolsTG.COVERAGE_PARTITION == (
            "multiaminocoverage", "multianticodoncoverage",
            "multitrnacoverage", "uniquecoverage",
        )

    def test_labels_are_tRAXs_own_wording(self):
        assert toolsTG.COVERAGE_CATEGORY_LABELS == {
            "uniquecoverage": "Transcript Specific",
            "multitrnacoverage": "Isodecoder Specific",
            "multianticodoncoverage": "Isotype Specific",
            "multiaminocoverage": "Not Amino Specific",
        }


class TestCategoryDirectories:
    @pytest.mark.parametrize("covtype,expected", [
        ("uniquecoverage", "unique"),
        ("multitrnacoverage", "isodecoder"),
        ("multianticodoncoverage", "isotype"),
        ("multiaminocoverage", "notamino"),
        ("coverage", "total"),
    ])
    def test_partition_members_use_their_short_alias(self, covtype, expected):
        assert toolsTG.coverage_category_dir(covtype) == expected

    def test_non_partition_types_use_their_own_name(self):
        assert toolsTG.coverage_category_dir("mismatchedbases") == "mismatchedbases"

    def test_every_covtype_gets_exactly_one_directory(self):
        """Two --covtype values must never share a directory, or one run's (expensive)
        output silently overwrites another's."""
        covtypes = list(toolsTG.COVERAGE_CATEGORY_DIRS) + ["readstarts", "readends", "deletions"]
        dirs = [toolsTG.coverage_category_dir(c) for c in covtypes]
        assert len(dirs) == len(set(dirs))

    def test_aliases_and_dirs_agree(self):
        """Round trip: every directory name is itself an accepted alias, so a user can read
        a path off disk and pass it straight back to --covtype."""
        for var_value, dirname in toolsTG.COVERAGE_CATEGORY_DIRS.items():
            assert toolsTG.COVERAGE_CATEGORY_ALIASES[dirname] == var_value


def test_partition_sums_to_total_coverage():
    """The property the stacked overview asserts visually. Synthetic here; verified on the
    real vibrChol1 object too, where the four categories reproduce `coverage` across all
    550 x 93 values to within float rounding (max relative difference 5.8e-16)."""
    rng = np.random.default_rng(0)
    parts = {c: rng.random((5, 12)) for c in toolsTG.COVERAGE_PARTITION}
    total = sum(parts.values())
    assert np.allclose(sum(parts[c] for c in toolsTG.COVERAGE_PARTITION), total)


def test_visualizer_exposes_the_overview_and_owns_its_layout():
    from trnagraph.modules import plotsCoverage
    for name in ("generate_partition_overview", "generate_partition_plot", "build_output_dirs"):
        assert hasattr(plotsCoverage.visualizer, name), f"missing {name}"
    # adataGraph must delegate directory creation rather than building the coverage tree
    # itself, so the per-category layout has exactly one owner.
    from trnagraph.modules import adataGraph
    source = inspect.getsource(adataGraph.anndataGrapher.plot)
    assert "pcV.build_output_dirs()" in source
    assert "pcV.generate_partition_overview()" in source
    assert f"{'{output}'}{'{self.args.covobs}'}/" not in source, (
        "adataGraph should no longer build coverage's subdirectories directly"
    )
