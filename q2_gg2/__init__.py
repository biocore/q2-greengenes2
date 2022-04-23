# ----------------------------------------------------------------------------
# Copyright (c) 2022-, Greengenes2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from ._methods import taxonomy_from_table, taxonomy_from_features
from . import _version
__version__ = _version.get_versions()['version']
__all__ = ['taxonomy_from_table', 'taxonomy_from_features']
