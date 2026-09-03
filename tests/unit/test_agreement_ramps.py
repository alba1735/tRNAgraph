"""Colour for the agreement tiers, one ramp per direction.

Two new `gradients` roles carry these. Left unset -- which is the normal case -- each direction's
ramp is DERIVED from the `colors` entry of the level it favours, so the agreement volcano uses
the same blue for Day 0 and the same yellow for Day 70 as every other figure in the set, with
nothing named in the style file.

More agreement is darker, following the `ordered` role's documented convention (light = least,
dark = most). The reference figure runs the other way, palest for the strongest tier; a style
file can have that back by naming the roles explicitly.
"""
import matplotlib
matplotlib.use('Agg')
import pytest
from matplotlib.colors import rgb_to_hsv, to_hex, to_rgba

from trnagraph.modules import plotsAgreement, plotsPalette


COLORS = {'Day 0': '#007FFF', 'Day 70': '#FFDC49'}


def test_a_ramp_gives_one_colour_per_tier():
    ramp = plotsAgreement.direction_ramp('Day 0', 3, colormap=COLORS)

    assert len(ramp) == 3
    assert len(set(ramp)) == 3


def test_the_strongest_tier_is_the_level_s_own_colour():
    """Derivation's whole point: the top tier is the colour that level wears elsewhere."""
    ramp = plotsAgreement.direction_ramp('Day 0', 3, colormap=COLORS)

    assert to_hex(ramp[-1]) == to_hex(COLORS['Day 0'])


def test_weaker_tiers_are_darker_than_stronger_ones():
    """The published figure's direction: the top tier wears the level's own light colour and
    weaker tiers darken. Blending toward WHITE was tried first and washed out entirely off a
    light base like a yellow, which has almost no headroom in lightness alone."""
    ramp = plotsAgreement.direction_ramp('Day 70', 4, colormap=COLORS)
    luminance = [0.299 * r + 0.587 * g + 0.114 * b for r, g, b, _ in map(to_rgba, ramp)]

    assert luminance == sorted(luminance), 'least agreement should be darkest'


def test_a_yellow_level_darkens_through_orange_toward_red():
    """Hue rotates as well as darkening, which is where the headroom comes from."""
    ramp = plotsAgreement.direction_ramp('Day 70', 3, colormap=COLORS)
    hues = [rgb_to_hsv(to_rgba(c)[:3])[0] for c in ramp]

    assert hues[-1] == pytest.approx(rgb_to_hsv(to_rgba(COLORS['Day 70'])[:3])[0])
    # Hue wraps past zero on the way to red, so compare the colours, not the raw numbers: the
    # middle tier must be an orange (red channel highest, blue lowest), never the cyan that a
    # naive wrapped interpolation produced.
    r, g, b, _ = to_rgba(ramp[1])
    assert r > g > b, 'the middle tier must be an orange'
    assert hues[0] > 0.85 or hues[0] < 0.05, 'the weakest tier is a red'


def test_a_blue_level_darkens_through_purple():
    ramp = plotsAgreement.direction_ramp('Day 0', 3, colormap=COLORS)
    hues = [rgb_to_hsv(to_rgba(c)[:3])[0] for c in ramp]

    assert hues[0] > hues[-1], 'blue rotates up toward purple'


def test_both_families_rotate_toward_the_same_end_of_the_wheel():
    """One rule, not a per-colour table: every ramp walks the short way round to magenta, which
    is what turns a yellow into a red and a blue into a purple without either being special."""
    yellow = rgb_to_hsv(to_rgba(plotsAgreement.direction_ramp(
        'Day 70', 3, colormap=COLORS)[0])[:3])[0]
    blue = rgb_to_hsv(to_rgba(plotsAgreement.direction_ramp(
        'Day 0', 3, colormap=COLORS)[0])[:3])[0]

    assert yellow < 0.10 or yellow > 0.85, 'yellow ends in the red/magenta arc'
    assert 0.70 < blue < 1.0, 'blue ends in the purple/magenta arc'


def test_a_single_tier_is_just_the_level_colour():
    """N=1 happens with two levels and one contrast; a ramp of one must not divide by zero."""
    ramp = plotsAgreement.direction_ramp('Day 0', 1, colormap=COLORS)

    assert len(ramp) == 1
    assert to_hex(ramp[0]) == to_hex(COLORS['Day 0'])


def test_an_explicit_gradient_role_overrides_the_derivation():
    from trnagraph.modules import toolsSchemas, toolsTG

    style = toolsSchemas.StyleFile(gradients={'agreement_down': 'Greens'})
    settings = toolsTG.resolve_plot_style(style, 'volcano')
    ramp = plotsAgreement.direction_ramp('Day 0', 3, colormap=COLORS,
                                         settings=settings, role='agreement_down')

    assert to_hex(ramp[-1]) != to_hex(COLORS['Day 0'])
    assert len(ramp) == 3


def test_a_level_with_no_colour_entry_falls_back_without_failing():
    """A style file need not name every level, and an unnamed one must still be drawable."""
    ramp = plotsAgreement.direction_ramp('Day 35', 3, colormap=COLORS, role='agreement_up')

    assert len(ramp) == 3
    assert len(set(ramp)) == 3


def test_no_style_file_at_all_still_produces_a_ramp():
    ramp = plotsAgreement.direction_ramp('Day 0', 2)

    assert len(ramp) == 2


def test_the_two_roles_are_registered_as_gradients():
    """resolve_gradients() hands every role to every plot, so an unregistered role would fail
    inside a worker rather than at load."""
    for role in ('agreement_up', 'agreement_down'):
        assert role in plotsPalette.GRADIENT_ROLES
    assert plotsPalette.gradient({}, 'agreement_up') is not None


def test_the_two_role_lists_cannot_drift_apart():
    """plotsPalette derives its roles from _GRADIENT_DEFAULTS and toolsSchemas declares them for
    validation. They were separate hand-written literals, and adding a role to one but not the
    other made gradient() raise KeyError for a role the style file happily accepted."""
    from trnagraph.modules import toolsSchemas

    assert set(plotsPalette.GRADIENT_ROLES) == set(toolsSchemas.GRADIENT_ROLES)


def test_the_fallback_ramp_runs_the_same_direction_as_the_derived_one():
    """Both must put the STRONGEST tier at the light end. They disagreed: the derived ramp
    darkened as agreement fell while the built-in fallback darkened as it rose, so the same
    legend row meant a different shade depending on whether a style file happened to name that
    level. Caught on real data, where the fallback was silently in use."""
    derived = plotsAgreement.direction_ramp('Day 70', 3, colormap=COLORS)
    fallback = plotsAgreement.direction_ramp('Unnamed', 3, role='agreement_up')

    def luminance(ramp):
        return [0.299 * r + 0.587 * g + 0.114 * b for r, g, b, *_ in map(to_rgba, ramp)]

    assert luminance(derived) == sorted(luminance(derived))
    assert luminance(fallback) == sorted(luminance(fallback))


def test_the_fallback_still_spans_a_visible_range():
    ramp = plotsAgreement.direction_ramp('Unnamed', 3, role='agreement_down')

    assert len(set(map(to_hex, ramp))) == 3


def test_settings_record_which_roles_the_style_file_actually_named():
    """resolve_gradients fills EVERY role with its default, so a resolved value is not evidence
    the user set anything. The derivation needs to tell those apart -- without this it saw a
    non-None colormap for agreement_up on every run and never derived from `colors` at all."""
    from trnagraph.modules import toolsSchemas, toolsTG

    style = toolsSchemas.StyleFile(gradients={'agreement_up': 'Greens'})
    settings = toolsTG.resolve_plot_style(style, 'volcano')

    assert settings['gradients_set'] == {'agreement_up'}


def test_no_style_file_names_no_roles():
    from trnagraph.modules import toolsTG

    assert toolsTG.resolve_plot_style(None, 'volcano')['gradients_set'] == set()


def test_the_derivation_runs_when_the_role_was_only_defaulted():
    from trnagraph.modules import toolsTG

    settings = toolsTG.resolve_plot_style(None, 'volcano')
    ramp = plotsAgreement.direction_ramp('Day 70', 3, colormap=COLORS, settings=settings,
                                         role='agreement_up')

    assert to_hex(ramp[-1]) == to_hex(COLORS['Day 70']), 'should derive, not fall back'


def test_an_explicitly_named_role_still_wins_over_the_derivation():
    from trnagraph.modules import toolsSchemas, toolsTG

    style = toolsSchemas.StyleFile(gradients={'agreement_up': 'Greens'})
    settings = toolsTG.resolve_plot_style(style, 'volcano')
    ramp = plotsAgreement.direction_ramp('Day 70', 3, colormap=COLORS, settings=settings,
                                         role='agreement_up')

    assert to_hex(ramp[-1]) != to_hex(COLORS['Day 70'])
