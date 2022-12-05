# ----------------------------------------------------------------------------
# Copyright (c) 2022-, Greengenes2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from ._methods import (taxonomy_from_table, taxonomy_from_features,
                       filter_features, relabel, compute_effect_size)
from ._version import get_versions
__version__ = get_versions()['version']
__all__ = ['taxonomy_from_table', 'taxonomy_from_features',
           'filter_features', 'relabel', 'compute_effect_size']
