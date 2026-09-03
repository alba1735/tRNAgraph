"""How many sets a diagram has decides how it is drawn, and whether it can be drawn at all.

Two and three sets are laid out AREA-PROPORTIONALLY by matplotlib-venn: the overlap is legible
without reading a number, which is most of the value of a Venn at that size. Four to six cannot
be drawn proportionally with circles at all, so they use fixed ellipses where the areas mean
nothing and every region is labelled with its count instead.

Six is where this stops. A 4- and 5-set Venn have known-good ellipse arrangements; SIX does
not -- it is why the reference implementations switch to triangles at that size -- and a
six-ellipse figure cannot represent all 63 regions, so counts would be silently absent from the
picture. That is the exact failure this module is built to avoid.

Rather than refuse the analysis, a diagram too large to draw honestly falls back to the UpSet
plot that is written alongside every complex diagram anyway, which stays exact at any number of
sets. The user still gets their answer; they get it in the representation that can carry it.

One set is refused outright: a Venn of one circle is a count, not a comparison.
"""
import pytest

from trnagraph.modules import plotsVenn


@pytest.mark.parametrize('n', [2, 3])
def test_small_diagrams_are_area_proportional(n):
    assert plotsVenn.venn_layout(n) == 'proportional'


@pytest.mark.parametrize('n', [4, 5])
def test_larger_diagrams_use_fixed_ellipses(n):
    assert plotsVenn.venn_layout(n) == 'ellipse'


@pytest.mark.parametrize('n', [6, 7, 12])
def test_beyond_five_falls_back_to_upset_rather_than_drawing_a_lie(n):
    """Six ellipses cannot show all 63 regions, so a count would be missing from the figure
    with nothing to say so. UpSet carries every region exactly, at any size."""
    assert plotsVenn.venn_layout(n) == 'upset_only'


def test_one_set_is_refused():
    with pytest.raises(Exception) as excinfo:
        plotsVenn.venn_layout(1)

    assert '1' in str(excinfo.value)
