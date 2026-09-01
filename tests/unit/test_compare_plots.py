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
import pytest

from trnagraph.modules import plotsCompare, toolsTG


AMINOS = ["Ala", "Arg", "Asn"]
ISOS = {"Ala": "AGC", "Arg": "ACG", "Asn": "GTT"}
GROUPS = ["ctrl", "drug"]
TIMEPOINTS = ["t0", "t1", "t2"]
REPLICATES = 2


def _make_adata():
    """Two comparegrp1 values x three comparegrp2 values x two replicates.

    Replicates matter: log2fc_compare_df takes a per-cell standard deviation and drops any row
    whose SD is NaN, so a single observation per cell silently empties the frame.
    """
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
                        # Fourfold per timepoint, so EVERY pair clears the |log2FC| >= 1
                        # threshold that decides whether the significance backdrop is drawn.
                        # Pairs straddling that threshold would draw different numbers of bars
                        # for legitimate reasons, which would blur the accumulation this test
                        # is looking for. Distinct per replicate so the SD is defined.
                        "nreads_total_unique_norm": 100.0 * (4 ** t_index) + replicate,
                        "nreads_total_norm": 100.0 * (4 ** t_index) + replicate,
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
        figure = plt.gcf()
        saved.append((path, len(figure.axes[0].patches)))

    monkeypatch.setattr(plotsCompare.toolsTG, "save_current", fake_save)
    plotsCompare.visualizer(_make_adata(), "group", "timepoint", None, "out/",
                            threaded=threaded)
    return saved


def test_each_comparison_is_drawn_on_its_own_axes(monkeypatch):
    """Three timepoints give three pairs per countgrp, and every pair plots the same number of
    bars. Under the bug the second figure carried the first's bars too and the third carried
    both, so the counts grew 1x, 2x, 3x down the file list."""
    saved = _run(monkeypatch)

    counts = [count for _, count in saved]
    assert len(set(counts)) == 1, (
        f"each figure must carry only its own bars; got patch counts {counts}")


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
