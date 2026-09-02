"""`-g compare` structural defects.

Three bugs in plotsCompare.visualizer, all structural rather than data-dependent, so none of
them needed unusual input to fire:

  * the figure was created once per countgrp, OUTSIDE the loop over comparison pairs, so each
    pair's bars were drawn onto axes still carrying the previous pair's;
  * the `return` sat INSIDE that loop, so a pooled run saved one figure and stopped, while the
    serial path (which takes the other branch) saved all of them -- two execution modes
    producing different output from the same command;
  * --comparegrp1 and --comparegrp2 both default to 'group', which pivots on a duplicated
    column and raised `ValueError: Grouper for 'group' not 1-dimensional` from inside pandas.

The fold change is taken BETWEEN --comparegrp2 values within each --comparegrp1 value
(log2fc_compare_df: log2(mdf[cgrp1][cgrp2[1]]) - log2(mdf[cgrp1][cgrp2[0]])), so --comparegrp1
is the series/colour axis. The CLI help and cli_reference described it the other way round.
"""
import anndata as ad
import numpy as np
import pandas as pd

from trnagraph.modules import plotsCompare


# Three aminos that move fourfold per timepoint and three that stay flat. The flat half is
# load-bearing now that the fold changes come from PyDESeq2: it re-derives size factors from
# the raw counts, so a dataset in which EVERY feature moves together is read as a library-size
# difference and normalized back to ~0. A stable majority is what anchors the normalization,
# which is the same reason test_log2fc.py's fixture carries flat tRNAs beside its differential
# one. All three pairs put the same three aminos over the |log2FC| >= 1 threshold, which is
# what keeps the per-figure bar count equal across pairs.
MOVING_AMINOS = ["Ala", "Arg", "Asn"]
FLAT_AMINOS = ["Asp", "Cys", "Gln"]
AMINOS = MOVING_AMINOS + FLAT_AMINOS
ISOS = {"Ala": "AGC", "Arg": "ACG", "Asn": "GTT", "Asp": "GTC", "Cys": "GCA", "Gln": "TTG"}
GROUPS = ["ctrl", "drug"]
TIMEPOINTS = ["t0", "t1", "t2"]
REPLICATES = 2


def _counts(amino, t_index, rng):
    """One cell's counts. Drawn from a negative binomial rather than set to a constant: with
    zero within-group variance the dispersion estimates collapse and apeGLM shrinkage returns
    nonsense (a flat feature came back at log2FC 2.0 instead of 0.0), which is an artifact of
    degenerate synthetic data rather than anything the pipeline does to real counts."""
    mean = 100.0 * (4 ** t_index) if amino in MOVING_AMINOS else 400.0
    value = float(rng.negative_binomial(n=20, p=20 / (20 + mean)))
    return {"nreads_total_unique_norm": value, "nreads_total_unique_raw": value,
            "nreads_total_norm": value, "nreads_total_raw": value}


def _make_adata():
    """Two comparegrp1 values x three comparegrp2 values x two replicates."""
    rng = np.random.default_rng(0)
    rows = []
    for amino in AMINOS:
        for group in GROUPS:
            for t_index, timepoint in enumerate(TIMEPOINTS):
                for replicate in range(REPLICATES):
                    rows.append({
                        "amino": amino,
                        "iso": ISOS[amino],
                        "trna": f"tRNA-{amino}-{ISOS[amino]}-1",
                        "sample": f"{group}_{timepoint}_{replicate}",
                        "group": group,
                        "timepoint": timepoint,
                        **_counts(amino, t_index, rng),
                    })
    obs = pd.DataFrame(rows)
    obs.index = [f"obs{i}" for i in range(len(obs))]
    var = pd.DataFrame({"coverage": ["uniquecoverage"], "gap": [False]}, index=["v0"])
    return ad.AnnData(X=np.zeros((len(obs), 1)), obs=obs, var=var)


def _run(monkeypatch, threaded=False):
    """Run the visualizer, recording what each save saw on the axes.

    Patching save_current is how the drawing surface is observed: 'the bars accumulated' is
    only visible at the moment a figure is written, and is invisible in the returned value or
    in the file list.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    saved = []

    def fake_save(path, settings=None):
        # Bars only, not the grey significance backdrop behind each significant row: that is
        # drawn with linewidth=0, and how many rows get one legitimately varies from pair to
        # pair once the fold changes come from a real dispersion fit. The bars themselves are
        # one per feature per --comparegrp1 value on every figure, which is what accumulation
        # would have multiplied.
        figure = plt.gcf()
        saved.append((path, sum(1 for p in figure.axes[0].patches if p.get_linewidth() > 0)))

    monkeypatch.setattr(plotsCompare.toolsTG, "save_current", fake_save)
    plotsCompare.visualizer(_make_adata(), "group", "timepoint", None, "out/",
                            threaded=threaded)
    return saved


def test_each_comparison_is_drawn_on_its_own_axes(monkeypatch):
    """Three timepoints give three pairs per countgrp, and every pair plots one bar per feature
    per --comparegrp1 value. Under the bug the second figure carried the first's bars too and
    the third carried both, so the counts grew 1x, 2x, 3x down the file list."""
    saved = _run(monkeypatch)

    counts = [count for _, count in saved]
    expected = len(AMINOS) * len(GROUPS)
    assert counts == [expected] * len(counts), (
        f"each figure must carry only its own {expected} bars; got {counts}")


def test_a_threaded_run_saves_every_comparison_and_reports_all_of_them(monkeypatch):
    """Guard, not a red->green cycle: the early `return` was removed by the same edit that
    moved the figure creation. It is pinned because the two execution modes silently produced
    different output sets -- the serial path took the `else` branch and saved everything, the
    pooled path returned after the first file -- and nothing outside this test compares them."""
    serial = _run(monkeypatch, threaded=False)
    threaded = _run(monkeypatch, threaded="Generating compare plots...\n")

    expected_pairs = len(TIMEPOINTS) * (len(TIMEPOINTS) - 1) // 2
    assert len(serial) == expected_pairs * 2, "one figure per pair, for amino and for iso"
    assert [path for path, _ in threaded] == [path for path, _ in serial]


def test_a_configured_line_width_reaches_the_bar_edges(monkeypatch):
    """compare's bar edges default to the bar spacing itself (`bardiff`), which is how the
    plot was tuned. A --style line_width replaces that with an absolute width; the grid
    hlines/vlines are furniture and keep theirs."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    widths = []

    def fake_save(path, settings=None):
        widths.extend(patch.get_linewidth() for patch in plt.gcf().axes[0].patches)

    monkeypatch.setattr(plotsCompare.toolsTG, "save_current", fake_save)
    plotsCompare.visualizer(_make_adata(), "group", "timepoint", None, "out/",
                            threaded=False, settings={"line_width": 0.25})

    assert widths, "no bars were drawn"
    # The zero-width entries are the grey significance backdrop drawn behind each row
    # (linewidth=0 by design). It is furniture, so line_width deliberately does not reach it;
    # every bar that actually has an edge carries the configured width.
    assert set(w for w in widths if w) == {0.25}
    assert 0.0 in widths, 'the strokeless backdrop must stay strokeless'


# ---------------------------------------------------------------------------------------------
# Re-basing compare on the same engine as the heatmap and volcano (roadmap.md). It used to call
# toolsTG.log2fc_compare_df -- a hand-rolled scipy.stats.ttest_ind_from_stats over per-group
# means and SDs of already-normalized counts -- so two plots of one dataset could report
# different fold changes for the same contrast, and its `sdf.dropna()` discarded every feature
# as soon as any group had a single replicate, because that group's SD is NaN.


def _make_single_replicate_adata():
    """One observation per comparegrp1/comparegrp2 cell. Real timecourse designs are routinely
    unreplicated at a timepoint, and a per-cell standard deviation does not exist there."""
    rows = []
    for amino in AMINOS:
        for group in GROUPS:
            for t_index, timepoint in enumerate(TIMEPOINTS):
                raw = 100.0 * (4 ** t_index)
                rows.append({
                    "amino": amino, "iso": ISOS[amino],
                    "trna": f"tRNA-{amino}-{ISOS[amino]}-1",
                    "sample": f"{group}_{timepoint}", "group": group, "timepoint": timepoint,
                    "nreads_total_unique_norm": raw, "nreads_total_unique_raw": raw,
                    "nreads_total_norm": raw, "nreads_total_raw": raw,
                })
    obs = pd.DataFrame(rows)
    obs.index = [f"obs{i}" for i in range(len(obs))]
    var = pd.DataFrame({"coverage": ["uniquecoverage"], "gap": [False]}, index=["v0"])
    return ad.AnnData(X=np.zeros((len(obs), 1)), obs=obs, var=var)


def test_a_single_replicate_per_cell_is_refused_by_name(monkeypatch, caplog):
    """An unreplicated design cannot be fitted by a negative-binomial model -- there are no
    residual degrees of freedom to estimate dispersion from -- so the honest outcome is a
    named refusal, not a plot. What changes is that it is now SAID: the t-test path's per-cell
    standard deviation was NaN, dropna() emptied the frame, and compare wrote nothing while
    reporting nothing. The check is made here rather than left to PyDESeq2's own ValueError,
    which names neither the column nor the stratum and aborts the entire graph run."""
    import logging

    import matplotlib
    matplotlib.use("Agg")

    saved = []
    monkeypatch.setattr(plotsCompare.toolsTG, "save_current",
                        lambda path, settings=None: saved.append(path))
    with caplog.at_level(logging.WARNING, logger="trnagraph.modules.plotsCompare"):
        plotsCompare.visualizer(_make_single_replicate_adata(), "group", "timepoint", None,
                                "out/", threaded=False)

    assert saved == [], "an unfittable design must not produce figures"
    messages = "\n".join(r.message for r in caplog.records)
    assert "no replicates to estimate dispersion" in messages
    assert "timepoint" in messages and "group" in messages, (
        f"the warning must name both columns; got: {messages}")
