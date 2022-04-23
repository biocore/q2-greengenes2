# ----------------------------------------------------------------------------
# Copyright (c) 2022-, Greengenes2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.plugin.model as model
from qiime2.plugin import ValidationError


# adapted from q2-demux ErrorCorrectionDetailsFmt
class ReferenceMapFmt(model.TextFileFormat):
    METADATA_COLUMNS = ['id', 'md5', 'sequence']

    def _validate_(self, level):
        line = open(str(self)).readline()
        if len(line.strip()) == 0:
            raise ValidationError("Failed to locate header.")

        header = set(line.strip().split('\t'))
        for column in sorted(self.METADATA_COLUMNS):
            if column not in header:
                raise ValidationError(f"{column} is not a column")


ReferenceMapDirFmt = model.SingleFileDirectoryFormat(
    'ReferenceMapDirFmt',
    'details.tsv',
    ReferenceMapFmt)
