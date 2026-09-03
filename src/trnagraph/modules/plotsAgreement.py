#!/usr/bin/env python3
'''
The agreement volcano: one contrast on the axis, cross-contrast agreement in colour.

A per-readtype volcano answers "did this tRNA change between these two conditions". It cannot
answer "did it change the same way at every timepoint", because each figure sees one contrast.
This plot puts one contrast on the axis and colours each feature by how many of the
reference-anchored contrasts it is significant in, so a consistent responder and a one-off
separate visibly.

Marker shape carries read type, so the three read types that partition the reads appear on one
figure rather than three.

Gated behind the `multivariate` --config block, like `-g venn`: the thresholds and the reference
level are choices about a particular experimental design.
'''

import logging
import os

import numpy as np
import pandas as pd

from . import plotsPalette

logger = logging.getLogger(__name__)


class InvalidAgreementError(ValueError):
    '''An agreement figure that cannot be built as specified.'''


def contrast_universe(levels, reference):
    '''
    The reference-anchored contrasts for `levels`, in the order the levels are given.

    Reference-anchored rather than all-pairwise: four timepoints give three contrasts against
    Day 0, not six combinations. Every contrast then shares a denominator, so "significant in 2
    of 3" describes one comparison repeated, not a mixture of unrelated ones.

    Level order is preserved rather than sorted, since an ordered category exists precisely so
    that figures read in experimental order.
    '''
    levels = list(levels)
    if reference not in levels:
        raise InvalidAgreementError(
            f'Reference level {reference!r} is not present in {levels}. Set '
            f'`multivariate.reference` to one of them, or leave it unset to take the first.')
    others = [level for level in levels if level != reference]
    if not others:
        raise InvalidAgreementError(
            f'{reference!r} is the only level present, so there is nothing to contrast it '
            f'against. An agreement figure needs at least two levels.')
    return [(reference, level) for level in others]


def _contrast_columns(contrast):
    '''The (log2FC, p-value) column names one contrast occupies in a toolsTG log2FC frame.'''
    name = f'{contrast[0]}-{contrast[1]}'
    return f'log2_{name}', f'pval_{name}'


def agreement_table(frames, contrasts, drawn, log2fc, padj, call_padj=None):
    '''
    One row per feature per read type: where it sits on the drawn contrast, and how many of the
    whole contrast universe agree with it.

    `frames` maps a read type to its log2FC frame; `drawn` is the contrast on the x-axis.

    The two thresholds do different jobs. `call_padj` (0.05 by default) decides whether the
    feature is CALLED on the contrast being drawn, which is what gives it a direction and a
    colour -- so no coloured point ever sits below that line. `padj` (0.001) is stricter and
    decides whether a contrast COUNTS toward agreement, which is the stronger claim the colour
    tier makes. A feature can therefore be called on this axis and agree with nothing, scoring
    0 of N: a real result on this contrast that was never strong anywhere.

    Agreement also requires the same DIRECTION as the drawn contrast -- a tRNA up at one
    timepoint and down at another has not been consistently anything. A feature not called on
    the drawn contrast has no direction and no count; it is still returned, and drawn in the
    neutral colour, because it is a real measurement rather than a missing one.

    Features DESeq2 could not call at all (NaN from independent filtering or Cook's-outlier
    detection) are dropped, following the volcano: they were never given a significance call,
    which is not the same as being insignificant.
    '''
    from . import plotsThresholds

    call_padj = plotsThresholds.PVAL_SIG if call_padj is None else call_padj
    reference, other = drawn
    rows = []
    for readtype, frame in frames.items():
        log2_col, pval_col = _contrast_columns(drawn)
        if log2_col not in frame.columns:
            raise InvalidAgreementError(
                f'The log2FC frame for {readtype!r} has no contrast {reference}-{other}. '
                f'It carries: {sorted(c[5:] for c in frame.columns if c.startswith("log2_"))}.')
        drawn_lfc = frame[log2_col]
        drawn_pval = frame[pval_col]
        callable_mask = drawn_lfc.notna() & drawn_pval.notna()

        # Significance and sign per contrast, aligned on the frame's index, so agreement is a
        # sum over booleans rather than a loop over features.
        significant, signs = {}, {}
        for contrast in contrasts:
            lfc_col, p_col = _contrast_columns(contrast)
            if lfc_col not in frame.columns:
                raise InvalidAgreementError(
                    f'The log2FC frame for {readtype!r} has no contrast '
                    f'{contrast[0]}-{contrast[1]}, which the agreement universe needs.')
            lfc, pval = frame[lfc_col], frame[p_col]
            significant[contrast] = (lfc.abs() >= log2fc) & (pval <= padj)
            signs[contrast] = np.sign(lfc)

        drawn_called = ((drawn_lfc.abs() >= log2fc) & (drawn_pval <= call_padj)
                        & callable_mask)
        drawn_sign = signs[drawn]
        agreeing = {contrast: (significant[contrast] & (signs[contrast] == drawn_sign)
                               & drawn_called)
                    for contrast in contrasts}
        n_agree = sum(agreeing.values()).astype(int)

        for feature in frame.index[callable_mask]:
            called = bool(drawn_called.loc[feature])
            if not called:
                direction = None
            else:
                direction = other if drawn_lfc.loc[feature] > 0 else reference
            rows.append({
                'feature': feature,
                'readtype': readtype,
                'log2fc': float(drawn_lfc.loc[feature]),
                'pval': float(drawn_pval.loc[feature]),
                'direction': direction,
                'n_agree': int(n_agree.loc[feature]),
                'n_contrasts': len(contrasts),
                'contrasts': ','.join(f'{a}-{b}' for (a, b) in contrasts
                                      if bool(agreeing[(a, b)].loc[feature])),
            })
    return pd.DataFrame(rows, columns=['feature', 'readtype', 'log2fc', 'pval', 'direction',
                                       'n_agree', 'n_contrasts', 'contrasts'])


def direction_ramp(level, n_tiers, colormap=None, settings=None, role='agreement_up'):
    '''
    The `n_tiers` colours for one direction, palest tier first.

    Derived from `colormap[level]` -- the same entry every other figure uses for that level --
    so a style file that names Day 0 blue gets a blue agreement ramp with nothing else declared.
    An explicit `gradients.agreement_up`/`agreement_down` overrides the derivation; a level with
    no colour entry falls back to that role's built-in ramp.

    The STRONGEST tier is the level's own colour, and weaker tiers darken while rotating toward
    magenta: yellow runs back through orange to red, blue forward through violet to purple. That
    is the published figure's direction, and it is also the only one with room to work -- a light
    base colour has almost no lightness left to give.
    '''
    # `gradients_set`, not `gradients`: every role resolves to a Colormap whether or not the
    # user named it, so reading the resolved value made this branch always win and the
    # derivation never ran.
    named = (settings or {}).get('gradients_set') or set()
    base = (colormap or {}).get(level)
    if role in named or base is None:
        # Spelled out rather than passed through as a variable: test_style_palette.py scans for
        # literal role names at gradient() call sites, so that a role a style file accepts but
        # no plot actually reads cannot slip through as a knob that does nothing.
        cmap = (plotsPalette.gradient(settings, 'agreement_up') if role == 'agreement_up'
                else plotsPalette.gradient(settings, 'agreement_down'))
        # REVERSED, and away from the extremes. Sequential colormaps run light to dark, but the
        # strongest tier is the LIGHT one here, so sampling forwards made the fallback run the
        # opposite way to the derived ramp -- the same legend row meaning a different shade
        # depending on whether the style file happened to name that level. The ends are avoided
        # because they are near-white and near-black, neither of which reads as a category.
        return [cmap(value) for value in np.linspace(0.95, 0.35, n_tiers)]

    # The derivation itself lives in plotsPalette, shared with the Venn's level x variant
    # circles: the same relationship between two colours has to look the same on both figures,
    # or a reader learns the encoding twice.
    return plotsPalette.related_ramp(base, n_tiers)


#: Colour for a feature with no significance call on the drawn contrast. Shared with the
#: volcano, so "grey means not called" reads the same way across the figure set.
UNCALLED_COLOR = plotsPalette.DE_NONSIGNIFICANT

#: Opacity for called and uncalled points, matching plotsVolcano so the two read alike.
_SIG_ALPHA = 0.9
_UNCALLED_ALPHA = 0.25


def _tier_colors(table, drawn, colormap, settings):
    '''
    A colour per row: the direction's ramp indexed by how many contrasts agree.

    Both ramps are built at the full tier depth so that "2 of 3" is the same shade whichever
    direction it belongs to -- a tier's colour has to mean the same thing on both sides, or the
    two halves of the legend are not comparable.

    Depth is N+1, not N: a feature called on the drawn axis at the loose threshold but strong in
    no contrast scores 0 of N, which is a real result and gets the palest stop rather than being
    dropped into the uncalled grey.
    '''
    from matplotlib.colors import to_rgba

    reference, other = drawn
    depth = max(int(table['n_contrasts'].max()), 1) + 1
    ramps = {
        other: direction_ramp(other, depth, colormap=colormap, settings=settings,
                              role='agreement_up'),
        reference: direction_ramp(reference, depth, colormap=colormap, settings=settings,
                                  role='agreement_down'),
    }
    colors = []
    for _, row in table.iterrows():
        if row['direction'] is None:
            colors.append(to_rgba(UNCALLED_COLOR, alpha=_UNCALLED_ALPHA))
        else:
            colors.append(to_rgba(ramps[row['direction']][int(row['n_agree'])],
                                  alpha=_SIG_ALPHA))
    return colors


def _label_rows(table, toplabels):
    '''
    Which rows get a text label: the most significant instance of each called feature.

    Priority is computed per POINT, since a tRNA can be a strong hit on one read type and
    marginal on another, but capped at one label per FEATURE -- otherwise a single tRNA drawn
    across three read types spends three of the label budget saying its own name.
    '''
    called = table[table['direction'].notna()]
    if called.empty or toplabels == 0:
        return called.iloc[0:0]
    score = called['log2fc'].abs() * -np.log10(
        called['pval'].astype(float).clip(lower=np.nextafter(0, 1)))
    ranked = called.assign(_score=score).sort_values('_score', ascending=False)
    best = ranked.drop_duplicates(subset='feature', keep='first')
    return best if toplabels is None else best.head(toplabels)


def draw_agreement(ax, table, drawn, log2fc, padj, colormap=None, settings=None,
                   toplabels=100, xlim=None, show_legend=True, call_padj=None):
    '''
    Draw one agreement volcano onto `ax`.

    `table` is agreement_table()'s output for the contrast being drawn. Marker shape carries
    read type and colour carries the agreement count, so the two encodings stay independent:
    a reader recovers "which reads" from the shape and "how consistent" from the colour.
    '''
    from . import plotsThresholds, plotsVolcano

    call_padj = plotsThresholds.PVAL_SIG if call_padj is None else call_padj
    if table.empty:
        raise InvalidAgreementError(
            'No features survived to plot for this contrast. Every feature was dropped as '
            'uncallable, or the read-count cutoff removed them all.')

    reference, other = drawn
    x = table['log2fc'].astype(float)
    # np.log10 rather than plotsThresholds.neglog10, which is math.log10 and scalar-only: it
    # exists to place a threshold LINE, and the thresholds themselves still go through it below.
    y = -np.log10(table['pval'].astype(float).clip(lower=np.nextafter(0, 1)))
    colors = _tier_colors(table, drawn, colormap, settings)

    cap, capped = plotsVolcano.resolve_x_limit(x, explicit=xlim)
    x_plot = x.clip(-cap, cap)
    offscale = x.abs() > cap
    settings = settings or {}
    marker_size = settings.get('marker_size') or 40
    threshold = settings.get('rasterize_over')
    rasterized = bool(threshold) and len(table) > threshold

    # One scatter call per read type, which is what makes the shape channel work.
    for readtype in dict.fromkeys(table['readtype']):
        marker = plotsPalette.readtype_marker(readtype)['marker']
        rows = (table['readtype'] == readtype).to_numpy() & ~offscale.to_numpy()
        if not rows.any():
            continue
        ax.scatter(x_plot[rows], y[rows], c=[c for c, keep in zip(colors, rows) if keep],
                   marker=marker, s=marker_size, linewidths=0, rasterized=rasterized)

    # Off-scale points keep the volcano's boundary-triangle treatment. Read type is recoverable
    # from position, and inventing a second set of glyphs for the edge would be worse (D15).
    if capped:
        for side, marker in plotsVolcano.OFFSCALE_MARKERS.items():
            rows = (offscale & ((x > 0) if side == 'right' else (x < 0))).to_numpy()
            if not rows.any():
                continue
            ax.scatter(x_plot[rows], y[rows], c=[c for c, keep in zip(colors, rows) if keep],
                       marker=marker, s=marker_size * 1.4, linewidths=0, rasterized=rasterized)

    for boundary in (-log2fc, log2fc):
        ax.axvline(x=boundary, color=plotsPalette.REFERENCE_LINE, linestyle='--')
    for pvalue in (plotsThresholds.PVAL_SIG, plotsThresholds.PVAL_STRONG):
        ax.axhline(y=plotsThresholds.neglog10(pvalue), color=plotsPalette.REFERENCE_LINE,
                   linestyle='--')

    ax.set_xlim(-cap * plotsVolcano.VOLCANO_AXIS_MARGIN, cap * plotsVolcano.VOLCANO_AXIS_MARGIN)
    labelled = _label_rows(table, toplabels)
    if not labelled.empty:
        names = pd.Index(labelled['feature'])
        positions_x = pd.Series(x_plot.to_numpy()[labelled.index], index=names)
        positions_y = pd.Series(y.to_numpy()[labelled.index], index=names)
        plotsVolcano._place_labels(ax, names, positions_x, positions_y)

    ax.set_xlabel('log2(fold change)')
    ax.set_ylabel('-log10(p-value)')

    if show_legend:
        ax.legend(handles=_legend_handles(table, drawn, colormap, settings, log2fc, padj,
                                          call_padj),
                  loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=7)
    return ax


def _legend_handles(table, drawn, colormap, settings, log2fc, padj, call_padj):
    '''
    Tiers as literal counts, then the read-type shapes, then the thresholds in force.

    Counts rather than invented tier names ("strong"/"moderate"): a name has to be redefined
    every time the number of contrasts changes, and "2 of 3" never stops being true.
    '''
    from matplotlib.lines import Line2D

    reference, other = drawn
    n_contrasts = max(int(table['n_contrasts'].max()), 1)
    depth = n_contrasts + 1
    handles = []
    for level, role in ((other, 'agreement_up'), (reference, 'agreement_down')):
        ramp = direction_ramp(level, depth, colormap=colormap, settings=settings, role=role)
        for tier in range(n_contrasts, -1, -1):
            handles.append(Line2D([0], [0], marker='o', linestyle='None', markersize=7,
                                  markerfacecolor=ramp[tier], markeredgecolor='none',
                                  label=f'{level} favored ({tier} of {n_contrasts})'))
    handles.append(Line2D([0], [0], marker='o', linestyle='None', markersize=7,
                          markerfacecolor=UNCALLED_COLOR, markeredgecolor='none',
                          alpha=_UNCALLED_ALPHA, label='Not significant'))

    seen = set()
    for readtype in dict.fromkeys(table['readtype']):
        spec = plotsPalette.readtype_marker(readtype)
        if spec['label'] in seen:
            continue
        seen.add(spec['label'])
        handles.append(Line2D([0], [0], marker=spec['marker'], linestyle='None', markersize=7,
                              markerfacecolor=plotsPalette.LEGEND_SHAPE_SWATCH,
                              markeredgecolor='none', label=spec['label']))

    # Both thresholds, each labelled with the job it does. Printing only the strict one while
    # colouring points down to the loose one is a caption that contradicts its own figure.
    handles.append(Line2D([0], [0], linestyle='None', marker='None',
                          label=f'colored: p <= {call_padj:g}, |log2(FC)| >= {log2fc:g}'))
    handles.append(Line2D([0], [0], linestyle='None', marker='None',
                          label=f'agreement counted: p <= {padj:g}'))
    return handles


#: The read type excluded from every multi-readtype overlay. It is the SUM of the others, so
#: drawing it beside them plots the same reads twice and a tRNA appears to agree with itself.
EXCLUDED_READTYPE = 'total'


def agreement_readtypes(readtypes):
    '''
    The read types this figure may carry: everything given except `total`.

    Excluded by construction rather than by documentation. Callers hold resolved obs columns
    (`nreads_total_unique_norm`), so the bare name is matched rather than the whole string.
    '''
    kept = [rt for rt in readtypes
            if plotsPalette.bare_readtype(rt) != EXCLUDED_READTYPE]
    if not kept:
        raise InvalidAgreementError(
            f'Every read type requested was {EXCLUDED_READTYPE!r}, which an agreement figure '
            f'cannot draw: it is the sum of the others, so plotting it beside them would count '
            f'the same reads twice. Name at least one of the read types that partition them.')
    return kept


def resolve_reference(levels, block):
    '''
    The level every contrast is measured against.

    An explicit `multivariate.reference` wins. Otherwise the FIRST category, which for an
    ordered categorical is the one `order` declared -- the same level DESeq2 treats as the
    reference, so the figure and the fit cannot disagree about which way a fold change points.
    '''
    import pandas as pd

    if isinstance(levels, pd.Categorical):
        ordered = list(levels.categories)
    else:
        ordered = list(dict.fromkeys(levels))
    if not ordered:
        raise InvalidAgreementError(
            'The grouping column has no levels, so there is nothing to contrast.')
    reference = getattr(block, 'reference', None)
    if reference is None:
        return ordered[0]
    if reference not in ordered:
        raise InvalidAgreementError(
            f'`multivariate.reference` names {reference!r}, which is not a level of the '
            f'grouping column. Available: {ordered}.')
    return reference


#: Columns the membership table carries, in order.
_TABLE_COLUMNS = ('feature', 'readtype', 'direction', 'n_agree', 'n_contrasts',
                  'contrasts', 'log2fc', 'pval')


def write_agreement_table(path, table, provenance):
    '''
    Write one contrast's called features to a TSV, most consistent first.

    The figure shows the shape; this is what a reader can act on -- which tRNAs agreed, which
    way, and across which contrasts, none of it recoverable off a scatter plot. Only CALLED
    features are listed: a table of every uncalled tRNA in the object is not a deliverable, it
    is the object.

    A commented provenance header names every parameter behind the numbers, so a file left in
    results/ after a rebuild identifies itself rather than quietly disagreeing with the object.
    Nothing reads this back -- the AnnData object stays the source of truth.
    '''
    called = table[table['direction'].notna()].copy()
    called = called.sort_values(['n_agree', 'pval'], ascending=[False, True])

    lines = ['# tRNAgraph agreement membership']
    lines += [f'# {key}: {value}' for key, value in provenance.items()]
    lines.append('\t'.join(_TABLE_COLUMNS))
    for _, row in called.iterrows():
        lines.append('\t'.join([
            str(row['feature']), str(row['readtype']), str(row['direction']),
            str(int(row['n_agree'])), str(int(row['n_contrasts'])),
            str(row['contrasts']), f"{row['log2fc']:.4f}", f"{row['pval']:.4g}",
        ]))
    path = str(path)
    with open(path, 'w') as handle:
        handle.write('\n'.join(lines) + '\n')
    return path


def _log2fc_frames(adata, grouping, readtypes, cutoff, config_name, overwrite=False,
                   shrink='apeGLM'):
    '''
    The cached log2FC frame for each read type, keyed by its bare name.

    The boundary to the DE engine, and the reason it is a function rather than inline: the fits
    are toolsTG's and are tested there, so an orchestration test stubs this instead of standing
    up a PyDESeq2-fittable fixture.

    n_cpus stays at 1. adataGraph precomputes these before its worker pool for exactly this
    reason -- fitting inside a pool deadlocks -- so by the time this runs the entries are cached
    and it is a lookup rather than a fit.
    '''
    from . import toolsTG

    frames = {}
    for readtype in readtypes:
        frame, _ = toolsTG.adataLog2FC(adata, grouping, readtype, readcount_cutoff=cutoff,
                                       config_name=config_name, overwrite=overwrite,
                                       n_cpus=1, shrink=shrink).main()
        frames[plotsPalette.bare_readtype(readtype)] = frame
    return frames


def visualizer(adata, block, output, config_name='default', settings=None, readtypes=None,
               cutoff=20, colormap=None, threaded=False, toplabels=100, xlim=None,
               overwrite=False, shrink='apeGLM', results_dir=None):
    '''
    Draw an agreement volcano for every reference-anchored contrast, and write their tables.

    One figure per contrast, each coloured by agreement across the whole contrast universe, plus
    a combined overview. Membership is computed once per contrast and BOTH the figure and the
    TSV render from it, so the two cannot disagree.
    '''
    import matplotlib.pyplot as plt

    from . import plotsThresholds, toolsTG

    grouping = block.grouping
    if grouping not in adata.obs.columns:
        raise InvalidAgreementError(
            f'`multivariate.grouping` names {grouping!r}, which is not a column of obs.')
    reference = resolve_reference(adata.obs[grouping], block)
    levels = (list(adata.obs[grouping].cat.categories)
              if hasattr(adata.obs[grouping], 'cat')
              else list(dict.fromkeys(adata.obs[grouping])))
    contrasts = contrast_universe(levels, reference)

    readtypes = agreement_readtypes(readtypes or [])
    frames = _log2fc_frames(adata, grouping, readtypes, cutoff, config_name,
                            overwrite=overwrite, shrink=shrink)

    individual_output = f'{output}individual/'
    messages = [toolsTG.builder(individual_output)]
    if results_dir:
        messages.append(toolsTG.builder(results_dir))

    provenance_base = {'config': config_name, 'cutoff': cutoff, 'grouping': grouping,
                       'reference': reference, 'log2fc': block.log2fc, 'padj': block.padj,
                       'called_padj': plotsThresholds.PVAL_SIG,
                       'readtypes': sorted(plotsPalette.bare_readtype(rt)
                                           for rt in readtypes)}

    tables = {}
    for contrast in contrasts:
        name = f'{contrast[0]}-{contrast[1]}'
        table = agreement_table(frames, contrasts, contrast,
                                log2fc=block.log2fc, padj=block.padj)
        tables[contrast] = table
        called = table[table['direction'].notna()]
        # Stored per CONTRAST, which is what the entry key means here -- an agreement figure
        # spans every contrast, so keying on a variant tag would make each overwrite the last.
        toolsTG.write_multivariate(
            adata, config_name, 'agreement', name,
            {str(direction): list(rows['feature'])
             for direction, rows in called.groupby('direction', observed=True)},
            dict(provenance_base, contrast=name))

        fig, ax = plt.subplots(figsize=toolsTG.figsize_for(settings, (7, 5)))
        draw_agreement(ax, table, contrast, log2fc=block.log2fc, padj=block.padj,
                       colormap=colormap, settings=settings, toplabels=toplabels, xlim=xlim)
        ax.set_title(f'{grouping}: {contrast[0]} vs {contrast[1]}', fontsize=10)
        path = f'{individual_output}{name}_agreement.pdf'
        fig.savefig(path, bbox_inches='tight')
        plt.close(fig)
        messages.append(f'Saving agreement volcano: {path}')

        if results_dir:
            table_path = os.path.join(results_dir, f'{name}_agreement.tsv')
            write_agreement_table(table_path, table,
                                  dict(provenance_base, contrast=name,
                                       object=adata.uns.get('trnagraphruninfo', {})))
            messages.append(f'Saving agreement table: {table_path}')

    messages.append(_draw_combined(tables, output, grouping, block, colormap, settings,
                                   toplabels, xlim))
    # Logged individually when running outside the pool, returned as one block when inside it --
    # the same split plotsVenn uses. Returning them only in the threaded path meant a sequential
    # run wrote every figure and table in silence.
    for message in messages:
        if not threaded:
            logger.info(message)
    return '\n'.join(messages) + '\n' if threaded else None


def _draw_combined(tables, output, grouping, block, colormap, settings, toplabels, xlim):
    '''
    One page carrying every contrast side by side.

    Deliberately NOT wired to the per-plot style knobs: combined pages compute their own
    geometry from how many panels they lay out, which is the convention the coverage and
    volcano overviews already follow.
    '''
    import matplotlib.pyplot as plt

    contrasts = list(tables)
    fig, axes = plt.subplots(1, len(contrasts), figsize=(7 * len(contrasts), 5), squeeze=False)
    for ax, contrast in zip(axes[0], contrasts):
        draw_agreement(ax, tables[contrast], contrast, log2fc=block.log2fc, padj=block.padj,
                       colormap=colormap, settings=settings, toplabels=toplabels, xlim=xlim,
                       show_legend=contrast is contrasts[-1])
        ax.set_title(f'{contrast[0]} vs {contrast[1]}', fontsize=10)
    path = f'{output}{grouping}_combined_agreement.pdf'
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    return f'Saving combined agreement overview: {path}'
