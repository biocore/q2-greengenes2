# ----------------------------------------------------------------------------
# Copyright (c) 2022-, Greengenes2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from .plugin_setup import plugin
from ._format import ReferenceMapFmt, CladeAssessmentFmt, ASVAssessmentFmt
from ._methods import CladeAssessment, ASVAssessment
import json
import pandas as pd


@plugin.register_transformer
def _1(ff: ReferenceMapFmt) -> pd.DataFrame:
    return pd.read_csv(str(ff), sep='\t')


@plugin.register_transformer
def _2(ff: CladeAssessmentFmt) -> CladeAssessment:
    d = json.loads(str(ff))
    return CladeAssessment(d)


@plugin.register_transformer
def _3(data: CladeAssessment) -> CladeAssessmentFmt:
    ff = CladeAssessmentFmt()
    with open(str(ff), 'w') as fp:
        fp.write(json.dumps(data, indent=2))
    return ff


@plugin.register_transformer
def _4(ff: ASVAssessmentFmt) -> ASVAssessment:
    d = json.loads(str(ff))
    return ASVAssessment(d)


@plugin.register_transformer
def _5(data: ASVAssessment) -> ASVAssessmentFmt:
    ff = ASVAssessmentFmt()
    with open(str(ff), 'w') as fp:
        fp.write(json.dumps(data, indent=2))
    return ff
