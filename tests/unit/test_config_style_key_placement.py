"""A key in the wrong file should say which file it belongs in.

`--config` and `--style` deliberately share no keys: a style file is meant to be reused
across differently-parameterized runs, which it could not be if it also fixed the cutoffs, so
selection/analysis options live in one file and presentation in the other. The cost of that
split is that a user reaching for the wrong file gets pydantic's bare "Extra inputs are not
permitted", which says nothing about where the key does belong.

Real report: `heatorient` put in style.json's `heatmap` block, which failed with
`extra_forbidden` and no indication that it is a --config flag.
"""
import json

import pytest

from trnagraph.modules import adataGraph, toolsTG


def _write(tmp_path, name, payload):
    path = tmp_path / name
    path.write_text(json.dumps(payload))
    return str(path)


def _load_config(path):
    from types import SimpleNamespace

    grapher = object.__new__(adataGraph.anndataGrapher)
    grapher.args = SimpleNamespace(config=path)
    grapher.logger = __import__('logging').getLogger('test')
    return grapher._load_config()


def test_a_graph_flag_in_a_style_file_says_it_belongs_in_the_config_file(tmp_path):
    path = _write(tmp_path, 'style.json', {'heatmap': {'heatorient': 'horizontal'}})

    with pytest.raises(ValueError) as raised:
        toolsTG.load_style_file(path)
    message = str(raised.value)
    assert 'heatorient' in message
    assert '--config' in message
    assert 'flags' in message


def test_a_style_key_in_a_config_file_says_it_belongs_in_the_style_file(tmp_path):
    path = _write(tmp_path, 'config.json',
                  {'name': 'x', 'flags': {'figsize': [4, 4]}})

    with pytest.raises(ValueError) as raised:
        _load_config(path)
    message = str(raised.value)
    assert 'figsize' in message
    assert '--style' in message


def test_a_genuine_typo_gets_a_near_match_rather_than_a_file_pointer(tmp_path):
    """`marker_sizes` belongs in neither file under that spelling, so the useful answer is
    the spelling that does exist."""
    path = _write(tmp_path, 'style.json', {'volcano': {'marker_sizes': 12}})

    with pytest.raises(ValueError) as raised:
        toolsTG.load_style_file(path)
    assert 'marker_size' in str(raised.value)


def test_the_original_validation_detail_is_kept(tmp_path):
    """The added sentence is guidance, not a replacement: the pydantic report still names the
    exact location, which matters when several keys are wrong at once."""
    path = _write(tmp_path, 'style.json', {'heatmap': {'heatorient': 'horizontal'}})

    with pytest.raises(ValueError) as raised:
        toolsTG.load_style_file(path)
    assert 'heatmap.heatorient' in str(raised.value)


# --- template ordering ------------------------------------------------------------------

import pathlib


@pytest.mark.parametrize('name', ['config.template.json', 'style.template.json'])
def test_every_template_object_is_alphabetically_ordered(name):
    """The templates are read as JSON, which cares nothing for order -- this is for the
    person reading the file and for the diff when a key is added. Regenerating a template by
    appending new keys leaves them at the bottom, which is how `corrmask` and `heatorient`
    ended up out of order."""
    template = json.loads(pathlib.Path('src/trnagraph/assets').joinpath(name).read_text())

    def check(obj, path='<root>'):
        if not isinstance(obj, dict):
            return
        keys = list(obj)
        assert keys == sorted(keys), f'{name}: {path} is not alphabetical: {keys}'
        for key, value in obj.items():
            check(value, f'{path}.{key}')

    check(template)


def test_a_malformed_file_is_a_usage_error_not_a_crash(tmp_path):
    """A file the user wrote is their input, so a mistake in it belongs on the terminal as a
    message rather than as a traceback -- the same treatment a typo'd flag gets. It still
    subclasses ValueError, so anything catching that keeps working."""
    path = _write(tmp_path, 'style.json', {'heatmap': {'heatorient': 'horizontal'}})

    with pytest.raises(toolsTG.UsageError):
        toolsTG.load_style_file(path)
    assert issubclass(toolsTG.UsageError, ValueError)
