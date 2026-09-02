"""Which features count as present in a group, for the expression-presence Venn.

A Venn circle is only as meaningful as the rule deciding what falls inside it, so this is the
most load-bearing choice in the plot. It reuses the vocabulary the DE path already applies --
adataLog2FC keeps features whose mean normalized count clears a cutoff -- so that a feature in a
Venn circle and a feature on a volcano passed the SAME test and the two figures can be read
against each other.

The difference is the axis the mean is taken over. adataLog2FC takes one mean across all groups
as a single global filter, because a fold change needs the feature present somewhere. Presence
is a per-group question: a tRNA detected at Day 0 and absent at Day 70 is exactly what the Venn
exists to show, and a global filter would place it in both circles or neither.

Thresholds deliberately are NOT inherited from the earlier notebook analyses (MIN_MEAN = 10,
MAX_DEPTH_CUTOFF = 20): those were fitted to older sequencing data with different size factors.
"""
import pandas as pd

from trnagraph.modules import toolsTG


def _obs():
    """Two tRNAs x two groups x two replicates.

    Ala is strongly present in both groups. Gly is present in ctrl only -- its treat replicates
    average 5, below any sensible cutoff.
    """
    rows = []
    for group, ala, gly in [('ctrl', 100.0, 80.0), ('treat', 120.0, 5.0)]:
        for replicate in range(2):
            rows.append({'trna': 'tRNA-Ala-AGC-1', 'sample': f'{group}_{replicate}',
                         'group': group, 'nreads_total_norm': ala})
            rows.append({'trna': 'tRNA-Gly-GCC-1', 'sample': f'{group}_{replicate}',
                         'group': group, 'nreads_total_norm': gly})
    return pd.DataFrame(rows)


def test_a_feature_present_in_one_group_only_lands_in_one_set():
    sets = toolsTG.presence_sets(_obs(), 'group', 'nreads_total_norm', cutoff=20)

    assert sets['ctrl'] == ['tRNA-Ala-AGC-1', 'tRNA-Gly-GCC-1']
    assert sets['treat'] == ['tRNA-Ala-AGC-1'], 'Gly averages 5 in treat, below the cutoff'


def test_every_group_gets_a_set_even_when_empty():
    """An empty circle is information -- nothing cleared the cutoff in that group -- and a
    missing key would instead read as 'this group was not analysed'."""
    sets = toolsTG.presence_sets(_obs(), 'group', 'nreads_total_norm', cutoff=1000)

    assert sets == {'ctrl': [], 'treat': []}


def test_membership_is_sorted_for_a_stable_stored_form():
    sets = toolsTG.presence_sets(_obs(), 'group', 'nreads_total_norm', cutoff=0)

    for members in sets.values():
        assert members == sorted(members)


def test_group_order_follows_a_declared_categorical():
    """The Venn's circles are labelled in this order, so it has to honour the config's declared
    order rather than pandas' discovery order."""
    obs = _obs()
    toolsTG.apply_category_order(obs, {'group': ['treat', 'ctrl']})

    sets = toolsTG.presence_sets(obs, 'group', 'nreads_total_norm', cutoff=20)

    assert list(sets) == ['treat', 'ctrl']
