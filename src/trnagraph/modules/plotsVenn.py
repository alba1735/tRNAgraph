#!/usr/bin/env python3
'''
Set-overlap diagrams over feature populations.

Two Venns are drawn automatically whenever `graph -g venn` runs and the data supports them:
fragment vs full-length (the u<N>/o<N> read-length split) and 5' vs 3' (the two end-specific
read types). Both answer the same shape of question -- is a given tRNA present as one species,
the other, or both -- which is why they are worth drawing for any dataset rather than being
declared per experiment.

The plot FAMILY is still gated behind a --config block: these are deliberate analyses, and a
user should not have them produced automatically and read conclusions off them.
'''

import logging
import os

import matplotlib.pyplot as plt

from . import toolsTG

from .toolsSchemas import VennPlan, VennSet

plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

logger = logging.getLogger(__name__)


class InvalidVennError(ValueError):
    '''A Venn that cannot be drawn as specified -- wrong number of sets.'''

#: The BARE read type the fragment-vs-full-length Venn is measured on. Total rather than an
#: end-specific count: the split already partitions by length, so the question is whether the
#: tRNA is seen at all within each length class.
#:
#: Bare, and resolved through toolsTG.resolve_readtype() against the run's read basis, so these
#: diagrams follow --allreads like every other graph type. Hardcoding a basis here would let a
#: Venn rest on a different denominator than the volcano beside it with nothing saying so.
SPLIT_READTYPE = 'total'

#: The two bare end-specific read types the second automatic Venn contrasts.
END_READTYPES = ('fiveprime', 'threeprime')


def draw_venn(ax, sets, colors=None):
    '''
    Draw a 2- or 3-set Venn onto `ax` and return matplotlib-venn's diagram object.

    Region counts come from exclusive_regions(), the same computation the TSV is written from,
    so the picture and the table beside it cannot report different numbers. matplotlib-venn is
    used here rather than the fixed-ellipse drawing needed at 4+ sets because at two and three
    it lays the circles out AREA-PROPORTIONALLY: the overlap is legible without reading a count.
    '''
    from matplotlib_venn import venn2, venn3

    labels = list(sets)
    if not 2 <= len(labels) <= 3:
        raise InvalidVennError(
            f'draw_venn handles 2 or 3 sets area-proportionally, got {len(labels)}. '
            f'Larger diagrams use the fixed-ellipse layout.')
    members = [set(sets[label]) for label in labels]
    draw = venn2 if len(labels) == 2 else venn3
    kwargs = {'set_colors': colors} if colors else {}
    return draw(members, set_labels=labels, ax=ax, **kwargs)


def resolve_results_dir(adata):
    '''
    Where this object's multivariate tables belong: (directory, skip message).

    `results/multivariate/` is a SIBLING of the per-variant result directories rather than a
    child of one, because a Venn spanning u60 and o60 belongs to neither and filing it under one
    would say something false about what produced it.

    The location comes from the object's own build provenance. When that directory is gone --
    the build output was moved, or built somewhere temporary, which is the case for the demo
    object -- this returns (None, message) and the caller skips the convenience copy. There is
    deliberately NO fallback directory: inventing one scatters output to wherever a graph run
    happened to be pointed, and the membership is still in the object, which is the source of
    truth regardless.
    '''
    runinfo = adata.uns.get('trnagraphruninfo') or {}
    build_dir = runinfo.get('trnagraph_directory') if hasattr(runinfo, 'get') else None
    if not build_dir:
        return None, ("Skipping the multivariate result tables: this object has no "
                      "uns['trnagraphruninfo']['trnagraph_directory'], so there is no results/ "
                      "tree to write them beside. The membership is stored in the object itself.")
    if not os.path.isdir(str(build_dir)):
        return None, (f"Skipping the multivariate result tables: the build directory recorded in "
                      f"this object, {build_dir}, no longer exists. The membership is stored in "
                      f"the object itself.")
    return os.path.join(str(build_dir), 'results', 'multivariate'), None


def exclusive_regions(sets):
    '''
    {region label: [feature, ...]} where every feature appears in exactly ONE region.

    A region is named by the full combination of sets containing the feature, joined with " & ",
    so a 2-set diagram yields "A", "B" and "A & B". This is the same membership an UpSet plot
    encodes -- which is what will let a complex Venn and its UpSet companion report identical
    numbers rather than two computations that merely ought to agree.

    Empty regions are omitted: a row reading `A & B  0` says nothing a reader needs, and with N
    sets there are 2**N - 1 possible regions, most of them empty on real data.
    '''
    labels = list(sets)
    members = {label: set(features) for label, features in sets.items()}
    regions = {}
    for feature in sorted(set().union(*members.values()) if members else set()):
        containing = tuple(label for label in labels if feature in members[label])
        if containing:
            regions.setdefault(' & '.join(containing), []).append(feature)
    return regions


def write_membership_table(path, sets, provenance):
    '''
    Write one Venn's membership to a TSV: which features fall in each exclusive region.

    The figure shows the counts; this is what a reader can act on. A commented provenance header
    names the object and every parameter behind the numbers, so a file left behind in results/
    after a rebuild identifies itself rather than quietly disagreeing with the object. results/
    is a convenience for people and papers -- the AnnData object remains the source of truth and
    nothing reads this back.
    '''
    regions = exclusive_regions(sets)
    lines = ['# tRNAgraph Venn membership']
    lines += [f'# {key}: {value}' for key, value in provenance.items()]
    lines.append('\t'.join(['region', 'n', 'features']))
    for region, features in regions.items():
        lines.append('\t'.join([region, str(len(features)), ','.join(features)]))
    path = str(path)
    with open(path, 'w') as handle:
        handle.write('\n'.join(lines) + '\n')
    return path


def require_multivariate_config(config):
    '''
    Return the `multivariate` config block, or refuse the run naming what is missing.

    These plots are opt-in by design (see adataGraph.GRAPH_TYPES_ALL). Refusing here, by name,
    is the difference between a user learning they need to declare the analysis and a user
    getting an empty output directory and no explanation.
    '''
    block = getattr(config, 'multivariate', None) if config is not None else None
    if block is None:
        raise InvalidVennError(
            "`-g venn` needs a `multivariate` block in your --config file, which declares the "
            "grouping column and thresholds the sets are built from. These analyses are "
            "deliberately opt-in rather than part of `-g all`: the sets and cutoffs are choices "
            "about your experiment, and figures produced from unchosen ones invite wrong "
            "conclusions. Minimal example: "
            '{"name": "my_analysis", "multivariate": {"grouping": "group"}}')
    return block


def _split_pairs(adata):
    '''Matched (under, over) split tags, keyed by cutoff -- u60 with o60, not with o50.'''
    tags = list(adata.uns.get('size_splits', {}))
    cutoffs = {}
    for tag in tags:
        if tag[:1] in ('u', 'o') and tag[1:].isdigit():
            cutoffs.setdefault(tag[1:], {})[tag[:1]] = tag
    return {cutoff: (pair['u'], pair['o'])
            for cutoff, pair in sorted(cutoffs.items()) if {'u', 'o'} <= set(pair)}


def simple_venn_plans(adata, read_basis=None, readtype=SPLIT_READTYPE):
    '''
    The automatic Venns this object can support, and a message for each it cannot.

    Returns (plans, skipped). A missing prerequisite is never an error: a plain build has no
    read-length split and a build without end-specific counts has no 5'/3' contrast, and in both
    cases the other Venn is still worth drawing. But the skip is REPORTED by name, because a
    figure that silently fails to appear is easily mistaken for a biological absence.
    '''
    plans, skipped = [], []
    basis = read_basis or toolsTG.READ_BASIS_UNIQUE
    readtype = toolsTG.resolve_readtype(readtype, basis, adata)
    end_readtypes = [toolsTG.resolve_readtype(rt, basis, adata) for rt in END_READTYPES]

    pairs = _split_pairs(adata)
    if not pairs:
        skipped.append(
            "Skipping the fragment_vs_full_length Venn: this object has no read-length split "
            "variant. Add one with `trnagraph analyze addsplit -c <cutoff>`, or build with "
            "--readlengthsplit.")
    elif readtype not in adata.obs.columns:
        skipped.append(f"Skipping the fragment_vs_full_length Venn: obs has no '{readtype}'.")
    else:
        for cutoff, (under, over) in pairs.items():
            name = ('fragment_vs_full_length' if len(pairs) == 1
                    else f'fragment_vs_full_length_{cutoff}')
            plans.append(VennPlan(
                name=name,
                title=f'Fragment (<{cutoff}nt) vs full-length (>={cutoff}nt) tRNAs',
                sets=[VennSet(label=f'Fragment ({under})', readtype=readtype, tag=under),
                      VennSet(label=f'Full length ({over})', readtype=readtype, tag=over)]))

    missing = [rt for rt in end_readtypes if rt not in adata.obs.columns]
    if missing:
        skipped.append(
            f"Skipping the fiveprime_vs_threeprime Venn: obs has no {missing}. End-specific "
            f"counts are written by `analyze build`; an object built before they existed, or "
            f"one filtered down to other read types, will not carry them.")
    else:
        plans.append(VennPlan(
            name='fiveprime_vs_threeprime',
            title="5' vs 3' tRNA fragments",
            sets=[VennSet(label="5' counts", readtype=end_readtypes[0]),
                  VennSet(label="3' counts", readtype=end_readtypes[1])]))

    return plans, skipped


if __name__ == '__main__':
    pass


def _set_members(adata, spec, cutoff):
    '''
    The features present in one circle.

    Presence is pooled across every sample, not split by group: these two automatic Venns ask
    whether a tRNA is seen at all as a fragment / as a 5' end, which is a property of the
    dataset rather than of one condition. A synthetic single-level grouping reuses
    presence_sets() rather than duplicating its cutoff rule.
    '''
    obs = toolsTG.variant_obs(adata, spec.tag)
    if spec.level is not None:
        obs = obs[obs[spec.level_column] == spec.level]
    obs = obs.assign(_pooled='all')
    return toolsTG.presence_sets(obs, '_pooled', spec.readtype, cutoff)['all']


def visualizer(adata, block, output, config_name='default', settings=None,
               read_basis=None, threaded=False):
    '''
    Draw every Venn this object supports, store the membership, and write the tables.

    Order matters: membership is computed once, stored on the object, and BOTH the figure and
    the TSV are rendered from that one computation, so the three can never disagree.
    '''
    individual_output = f'{output}individual/'
    messages = []
    plans, skipped = simple_venn_plans(adata, read_basis=read_basis)
    for message in skipped:
        messages.append(message)
        logger.warning(message)
    if not plans:
        return '\n'.join(messages) + '\n' if threaded else None

    results_dir, results_message = resolve_results_dir(adata)
    if results_message:
        messages.append(results_message)
        logger.warning(results_message)
    else:
        os.makedirs(results_dir, exist_ok=True)

    messages.append(toolsTG.builder(individual_output))
    provenance_base = {'config': config_name, 'presence_cutoff': block.presence_cutoff,
                       'membership': 'expression-presence'}

    for plan in plans:
        sets = {spec.label: _set_members(adata, spec, block.presence_cutoff)
                for spec in plan.sets}
        provenance = dict(provenance_base, plan=plan.name,
                          readtypes=sorted({s.readtype for s in plan.sets}),
                          tags=sorted({s.tag for s in plan.sets}))
        # Keyed by the DIAGRAM, not by a variant tag: a Venn spans variants, so every one
        # in a run would otherwise share the same key and overwrite its predecessor.
        toolsTG.write_multivariate(adata, config_name, 'venn', plan.name, sets, provenance)

        fig, ax = plt.subplots(figsize=toolsTG.figsize_for(settings, (6, 6)))
        draw_venn(ax, sets)
        ax.set_title(f'{plan.title}\n(present at mean >= {block.presence_cutoff:g} normalized reads)',
                     fontsize=10)
        path = f'{individual_output}{plan.name}_venn.pdf'
        fig.savefig(path, bbox_inches='tight')
        plt.close(fig)
        messages.append(f'Saving Venn: {path}')

        if results_dir:
            table = os.path.join(results_dir, f'{plan.name}_venn.tsv')
            write_membership_table(table, sets, dict(provenance, object=adata.uns.get(
                'trnagraphruninfo', {}).get('expname', 'unknown')))
            messages.append(f'Saving Venn membership table: {table}')

    for message in messages:
        if not threaded:
            logger.info(message)
    return '\n'.join(messages) + '\n' if threaded else None
