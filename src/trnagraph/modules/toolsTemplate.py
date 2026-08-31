#!/usr/bin/env python3
'''
Write a blank, fully-enumerated copy of a tRNAgraph JSON config into the working directory.

The style file's surface is large -- a colors block, six gradient roles, a categorical
palette, and presentation settings that can be set globally or per graph type -- and JSON
carries no comments, so there is nowhere in the format itself to document what a file may
contain. The templates are that documentation: every key the schema accepts, present and set
to null.

Emitting them through a command rather than telling users to copy the packaged file is
deliberate. The assets live inside the installed package (`<site-packages>/trnagraph/assets/`),
which a user who ran a plain `pip install` cannot reasonably be expected to find, and a
previous attempt to locate assets by walking directories up from `__file__` resolved to a
path that does not exist under a non-editable install. toolsTG.assets_dir() resolves through
importlib.resources instead, which is correct for both layouts.

Every value being null matters beyond readability: a template is a no-op if passed straight
to `--style`, so nothing in it is silently load-bearing, and a unit test asserts exactly that
alongside the template's agreement with the schema.
'''

import logging
import os
import shutil

from . import toolsTG

# Template name -> (packaged asset filename, what it configures). The `config` template is
# added alongside the graph-flags work; keeping this a registry means the command's shape
# does not change when it arrives.
TEMPLATES = {
    'style': ('style.template.json', 'colors, gradients and presentation settings (--style)'),
}


class TemplateWriter:
    '''Copy blank config templates out of the installed package into a target directory.'''

    def __init__(self, args):
        self.args = args
        self.logger = logging.getLogger(__name__)

    def selected(self):
        '''
        Which templates to write: the ones named by flags, or all of them when none are named.

        Writing everything on a bare `trnagraph tools template` is the discoverable default --
        a user who does not yet know what the files are called cannot be expected to name one.
        '''
        chosen = [name for name in TEMPLATES if getattr(self.args, name, False)]
        return chosen or list(TEMPLATES)

    def run(self):
        output = os.path.abspath(getattr(self.args, 'output', None) or '.')
        os.makedirs(output, exist_ok=True)
        written = []
        for name in self.selected():
            asset_name, description = TEMPLATES[name]
            source = os.path.join(toolsTG.assets_dir(), asset_name)
            if not os.path.isfile(source):
                raise FileNotFoundError(
                    f'Packaged template {asset_name} is missing from {toolsTG.assets_dir()}. '
                    f'This indicates a broken install rather than a usage error.'
                )
            destination = os.path.join(output, asset_name)
            # Refusing by default rather than overwriting: the working pattern for these files
            # is emit, edit, then re-emit after an upgrade, and silently replacing an edited
            # style file would destroy real work with no way back.
            if os.path.exists(destination) and not getattr(self.args, 'overwrite', False):
                raise FileExistsError(
                    f'{destination} already exists. Pass --overwrite to replace it, or write '
                    f'to a different directory with -o.'
                )
            shutil.copyfile(source, destination)
            self.logger.info(f'Wrote {destination} -- {description}')
            written.append(destination)
        return written
