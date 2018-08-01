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

import os
import shlex
import subprocess

from utils import logger

from basic_modules.metadata import Metadata
from basic_modules.tool import Tool
from mg_common.tool.bam_utils import bamUtilsTask


# ------------------------------------------------------------------------------

class Macs2(Tool):
    """
    Tool for peak calling for ChIP-seq data
    """

    def __init__(self, configuration=None):
        """
        Init function
        """
        logger.info("MACS2 Peak Caller")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @staticmethod
    def macs2_peak_calling(  # pylint: disable=too-many-locals,too-many-statements,too-many-statements,too-many-arguments
            name, bam_file, bai_file,
            narrowpeak, summits_bed, broadpeak, gappedpeak,
            bam_file_bgd=None, bai_file_bgd=None):
        """
        Function to run MACS2 for peak calling on aligned sequence files and
        normalised against a provided background set of alignments.

        Parameters
        ----------
        name : str
            Name to be used to identify the files
        bam_file : str
            Location of the aligned FASTQ files as a bam file
        bai_file : str
            Location of the bam index file
        narrowpeak : str
            Location of the output narrowpeak file
        summits_bed : str
            Location of the output summits bed file
        broadpeak : str
            Location of the output broadpeak file
        gappedpeak : str
            Location of the output gappedpeak file
        bam_file_bgd : str
            Location of the aligned FASTQ files as a bam file representing
            background values for the cell
        bai_file_bgd : str
            Location of the background bam index file

        Returns
        -------
        narrowPeak : file
            BED6+4 file - ideal for transcription factor binding site
            identification
        summitPeak : file
            BED4+1 file - Contains the peak summit locations for everypeak
        broadPeak : file
            BED6+3 file - ideal for histone binding site identification
        gappedPeak : file
            BED12+3 file - Contains a merged set of the broad and narrow peak
            files

        Definitions defined for each of these files have come from the MACS2
        documentation described in the docs at https://github.com/taoliu/MACS
        """
        od_list = bam_file.split("/")
        output_dir = "/".join(od_list[0:-1])

        with open(narrowpeak, "w") as f_out:
            f_out.write("")
        with open(summits_bed, "w") as f_out:
            f_out.write("")
        with open(broadpeak, "w") as f_out:
            f_out.write("")
        with open(gappedpeak, "w") as f_out:
            f_out.write("")

        from mg_common.tool.bam_utils import bamUtils

        bam_utils_handle = bamUtils()

        command_param = [
            'macs2 callpeak --nomodel', '-t', bam_file, '-n', name
        ]
        if bam_file_bgd is not None:
            bgd_command = '-c ' + bam_file_bgd
            command_param.append(bgd_command)

        command_param.append('--outdir ' + output_dir)
        command_line = ' '.join(command_param)

        if int(bam_utils_handle.bam_count_reads(bam_file, aligned=True)) > 0:
            try:
                args = shlex.split(command_line)
                process = subprocess.Popen(args)
                process.wait()
            except (IOError, OSError) as msg:
                logger.fatal("I/O error({0}): {1}\n{2}".format(
                    msg.errno, msg.strerror, command_line))
                return False

            if process.returncode is not 0:
                logger.fatal("MACS2 ERROR", process.returncode)
                return False

            logger.info('Process Results 1:', process)

        logger.info('LIST DIR 1:', os.listdir(output_dir))

        out_suffix = ['peaks.narrowPeak', 'peaks.broadPeak', 'peaks.gappedPeak', 'summits.bed']
        for f_suf in out_suffix:
            output_tmp = output_dir + '/' + name + '_' + f_suf
            logger.info(output_tmp, os.path.isfile(output_tmp))
            if os.path.isfile(output_tmp) is True and os.path.getsize(output_tmp) > 0:
                if f_suf == 'peaks.narrowPeak':
                    with open(narrowpeak, "wb") as f_out:
                        with open(output_tmp, "rb") as f_in:
                            f_out.write(f_in.read())
                elif f_suf == 'summits.bed':
                    with open(summits_bed, "wb") as f_out:
                        with open(output_tmp, "rb") as f_in:
                            f_out.write(f_in.read())
                elif f_suf == 'peaks.broadPeak':
                    with open(broadpeak, "wb") as f_out:
                        with open(output_tmp, "rb") as f_in:
                            f_out.write(f_in.read())
                elif f_suf == 'peaks.gappedPeak':
                    with open(gappedpeak, "wb") as f_out:
                        with open(output_tmp, "rb") as f_in:
                            f_out.write(f_in.read())

        return True

    def run(self, input_files, input_metadata, output_files):  # pylint: disable=too-many-locals,too-many-statements
        """
        The main function to run MACS 2 for peak calling over a given BAM file
        and matching background BAM file.

        Parameters
        ----------
        input_files : list
            List of input bam file locations where 0 is the bam data file and 1
            is the matching background bam file
        metadata : dict


        Returns
        -------
        output_files : list
            List of locations for the output files.
        output_metadata : list
            List of matching metadata dict objects

        """
        root_name = input_files['bam'].split("/")
        root_name[-1] = root_name[-1].replace('.bam', '')
        name = root_name[-1]

        # input and output share most metadata
        output_bed_types = {
            'narrow_peak': "bed4+1",
            'summits': "bed6+4",
            'broad_peak': "bed6+3",
            'gapped_peak': "bed12+3"
        }

        bam_utils_handle = bamUtilsTask()
        bam_utils_handle.bam_index(
            input_files['bam'],
            input_files['bam'] + '.bai'
        )
        if 'bam_bg' in input_files:
            bam_utils_handle.bam_index(
                input_files['bam_bg'],
                input_files['bam_bg'] + '.bai'
            )

        if 'bam_bg' in input_files:
            result = self.macs2_peak_calling(
                name,
                str(input_files['bam']), str(input_files['bam']) + '.bai',
                str(output_files['narrow_peak']),
                str(output_files['summits']),
                str(output_files['broad_peak']),
                str(output_files['gapped_peak']),
                str(input_files['bam_bg']), str(input_files['bam_bg']) + '.bai'
            )
        else:
            result = self.macs2_peak_calling(
                name,
                str(input_files['bam']), str(input_files['bam']) + '.bai',
                str(output_files['narrow_peak']),
                str(output_files['summits']),
                str(output_files['broad_peak']),
                str(output_files['gapped_peak'])
            )

        if result is False:
            logger.fatal("MACS2: Something went wrong with the peak calling")

        output_files_created = {}
        output_metadata = {}
        for result_file in output_files:
            if (
                    os.path.isfile(output_files[result_file]) is True
                    and os.path.getsize(output_files[result_file]) > 0
            ):
                output_files_created[result_file] = output_files[result_file]

                sources = [input_metadata["bam"].file_path]
                if 'bam_bg' in input_files:
                    sources.append(input_metadata["bam_bg"].file_path)

                output_metadata[result_file] = Metadata(
                    data_type="data_chip_seq",
                    file_type="BED",
                    file_path=output_files[result_file],
                    sources=sources,
                    taxon_id=input_metadata["bam"].taxon_id,
                    meta_data={
                        "assembly": input_metadata["bam"].meta_data["assembly"],
                        "tool": "macs2",
                        "bed_type": output_bed_types[result_file]
                    }
                )
            else:
                os.remove(output_files[result_file])

        logger.info('MACS2: GENERATED FILES: ', ' '.join(output_files))

        return (output_files_created, output_metadata)

# ------------------------------------------------------------------------------
