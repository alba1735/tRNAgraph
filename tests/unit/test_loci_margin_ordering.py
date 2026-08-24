"""Pre-tRNA position numbers must be computed before the loci margin is added.

`gettnanums()` adds the edge margin itself, emitting `head{N}..head1` before the
consensus and `tail1..tail{N}` after it. Passing it an alignment that has already
been margined therefore applies the margin twice: `add_margin()`'s padding
characters are not in the `+=*` structure set, so the consensus loop consumes
them as ordinary gap positions, the returned list overruns the alignment width,
and the trailing `tail` labels drop off the end when it is zipped against the
coverage array. The `head` labels survive because they sit at the front.

tRAX gets the order right -- `getcoverage.py:700` computes the numbers, then
`:707` margins the alignment. tRNAgraph had margined first, which cost every
pre-tRNA its entire 3' flanking region: on the human hg38 dataset
`pretRNAcoverage.txt` carried all 30 `head` positions but zero of the 30 `tail`
positions tRAX emits (222 distinct positions against tRAX's 252).
"""
from trnagraph.modules import toolsTG
from trnagraph.modules.toolsGetCoverage import gettnanums

MARGIN = 3


def _alignment(consensus):
    return toolsTG.RnaAlignment(
        {"trna1": "A" * len(consensus)}, ":" * len(consensus), consensus=consensus
    )


def test_unmargined_alignment_yields_both_head_and_tail_labels():
    nums = gettnanums(_alignment("=" * 20), margin=MARGIN, orgtype="euk")
    assert nums[:MARGIN] == ["head3", "head2", "head1"]
    assert nums[-MARGIN:] == ["tail1", "tail2", "tail3"]


def test_pre_margined_alignment_overruns_and_loses_the_tail():
    align = _alignment("=" * 20)
    correct = gettnanums(align, margin=MARGIN, orgtype="euk")
    doubled = gettnanums(align.add_margin(MARGIN), margin=MARGIN, orgtype="euk")
    # add_margin pads both sides, so the consensus loop emits 2*MARGIN extra entries
    assert len(doubled) == len(correct) + 2 * MARGIN
    # ...which is what pushes the tail labels past the coverage array's width
    assert doubled[len(correct) - MARGIN:len(correct)] != ["tail1", "tail2", "tail3"]


def test_main_margins_the_loci_alignment_after_numbering_it():
    """Guards the call-site ordering that actually produced the bug."""
    import inspect
    from trnagraph.modules import toolsGetCoverage
    src = inspect.getsource(toolsGetCoverage.main)
    add_margin_at = src.index("locistk_obj.add_margin(lociedgemargin)")
    numbering_at = src.index("locipositionnums = gettnanums(")
    assert numbering_at < add_margin_at, (
        "locistk_obj must be margined only after gettnanums() has numbered it"
    )
