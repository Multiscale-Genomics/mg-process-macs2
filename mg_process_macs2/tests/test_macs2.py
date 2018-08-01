"""
.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

from __future__ import print_function

import os.path
import pytest

from basic_modules.metadata import Metadata
from mg_process_macs2.tool.macs2 import Macs2


@pytest.mark.chipseq
def test_macs2():
    """
    Function to test MACS2
    """

    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "bam": resource_path + "macs2.Human.DRR000150.22_aln_filtered.bam"
    }

    output_files = {
        "narrow_peak": resource_path + "macs2.Human.DRR000150.22_peaks.narrowPeak",
        "summits": resource_path + "macs2.Human.DRR000150.22_peaks.summits.bed",
        "broad_peak": resource_path + "macs2.Human.DRR000150.22_peaks.broadPeak",
        "gapped_peak": resource_path + "macs2.Human.DRR000150.22_peaks.gappedPeak"
    }

    metadata = {
        "bam": Metadata(
            "data_chipseq", "fastq", [], None,
            {'assembly': 'test'}),
    }

    macs_handle = Macs2()
    macs_handle.run(input_files, metadata, output_files)

    assert os.path.isfile(resource_path + "macs2.Human.DRR000150.22_peaks.narrowPeak") is True
    assert os.path.getsize(resource_path + "macs2.Human.DRR000150.22_peaks.narrowPeak") > 0
    assert os.path.isfile(resource_path + "macs2.Human.DRR000150.22_peaks.summits.bed") is True
    assert os.path.getsize(resource_path + "macs2.Human.DRR000150.22_peaks.summits.bed") > 0

    os.remove(resource_path + "macs2.Human.DRR000150.22_peaks.narrowPeak")
    os.remove(resource_path + "macs2.Human.DRR000150.22_peaks.summits.bed")
    os.remove(resource_path + "macs2.Human.DRR000150.22_aln_filtered.bam.bai")
    os.remove(resource_path + "macs2.Human.DRR000150.22_aln_filtered.bam_tmp.bai")
    os.remove(resource_path + "macs2.Human.DRR000150.22_aln_filtered_peaks.narrowPeak")
    os.remove(resource_path + "macs2.Human.DRR000150.22_aln_filtered_peaks.xls")
    os.remove(resource_path + "macs2.Human.DRR000150.22_aln_filtered_summits.bed")
