"""The per-feature table the agreement figure writes beside itself.

D18: the figure shows the shape, this is the deliverable. Which tRNAs agreed, which way they
moved, and which contrasts they were significant in -- none of it readable off a scatter plot.
A commented header names every parameter behind the numbers, so a file left in results/ after a
rebuild identifies itself rather than quietly disagreeing with the object.
"""
import pandas as pd

from trnagraph.modules import plotsAgreement


REF, OTHER = 'Day 0', 'Day 70'
CONTRASTS = [(REF, 'Day 35'), (REF, OTHER)]


def _table():
    features = ['up-both', 'down-one', 'quiet']
    data = {}
    for (a, b), (lfc, p) in {
        (REF, 'Day 35'): ([4.0, -0.2, 0.1], [1e-8, 0.7, 0.9]),
        (REF, OTHER): ([5.0, -4.0, 0.2], [1e-9, 1e-9, 0.8]),
    }.items():
        data[f'log2_{a}-{b}'] = pd.Series(lfc, index=features)
        data[f'pval_{a}-{b}'] = pd.Series(p, index=features)
    frame = pd.DataFrame(data, index=features)
    return plotsAgreement.agreement_table({'fiveprime': frame}, CONTRASTS, (REF, OTHER),
                                          log2fc=1.5, padj=0.001)


def _write(tmp_path, **provenance):
    path = tmp_path / 'agreement.tsv'
    provenance = provenance or {'config': 'test', 'cutoff': 20}
    plotsAgreement.write_agreement_table(path, _table(), provenance)
    return path.read_text()


def test_only_called_features_are_listed(tmp_path):
    """A table of every uncalled tRNA in the object is not a deliverable, it is the object."""
    body = [line for line in _write(tmp_path).splitlines() if not line.startswith('#')]

    assert body[0].split('\t')[0] == 'feature'
    assert [row.split('\t')[0] for row in body[1:]] == ['up-both', 'down-one']


def test_each_row_carries_direction_and_the_contrasts_it_agreed_in(tmp_path):
    rows = {line.split('\t')[0]: line.split('\t')
            for line in _write(tmp_path).splitlines() if not line.startswith('#')}
    header = [line.split('\t') for line in _write(tmp_path).splitlines()
              if line.startswith('feature')][0]

    up = dict(zip(header, rows['up-both']))
    assert up['direction'] == OTHER
    assert up['n_agree'] == '2'
    assert f'{REF}-{OTHER}' in up['contrasts']


def test_a_feature_that_agreed_with_nothing_else_says_so(tmp_path):
    rows = {line.split('\t')[0]: line.split('\t')
            for line in _write(tmp_path).splitlines() if not line.startswith('#')}
    header = [line.split('\t') for line in _write(tmp_path).splitlines()
              if line.startswith('feature')][0]

    down = dict(zip(header, rows['down-one']))
    assert down['n_agree'] == '1'
    assert down['contrasts'] == f'{REF}-{OTHER}'


def test_the_header_names_every_parameter_behind_the_numbers(tmp_path):
    text = _write(tmp_path, config='basic', cutoff=20, log2fc=1.5, padj=0.001,
                  reference=REF, contrast=f'{REF}-{OTHER}')
    header = [line for line in text.splitlines() if line.startswith('#')]

    joined = ' '.join(header)
    for expected in ('basic', '20', '1.5', '0.001', REF):
        assert expected in joined


def test_rows_are_ordered_most_consistent_first(tmp_path):
    """A reader opens this looking for the strongest agreement, not for alphabetical order."""
    body = [line.split('\t') for line in _write(tmp_path).splitlines()
            if not line.startswith('#')][1:]

    assert [row[0] for row in body] == ['up-both', 'down-one']
