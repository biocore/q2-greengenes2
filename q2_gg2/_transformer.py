# ----------------------------------------------------------------------------
# Copyright (c) 2022-, Greengenes2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from .plugin_setup import plugin
from ._format import ReferenceMapFmt

import pandas as pd


@plugin.register_transformer
def _1(ff: ReferenceMapFmt) -> pd.DataFrame:
    return pd.read_csv(str(ff), sep='\t')
