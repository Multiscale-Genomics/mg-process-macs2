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

Requirements and Installation
=============================

Requirements
------------

Software
^^^^^^^^

- Python 2.7.10+
- Cython 0.25+
- HTSlib
- SAMtools

.. code-block:: none
   :linenos:

   sudo apt-get install make build-essential zlib1g-dev libbz2-dev liblzma-dev curl pigz

   cd ${HOME}/lib
   git clone https://github.com/samtools/htslib.git
   cd htslib
   autoheader
   autoconf
   ./configure --prefix=${HOME}/lib/htslib
   make
   make install

   cd ${HOME}/lib
   git clone https://github.com/samtools/samtools.git
   cd samtools
   autoheader
   autoconf -Wno-syntax
   ./configure --prefix=${HOME}/lib/samtools
   make
   make install

The following will then need to be on your $PATH:

.. code-block:: none
   :linenos:

   cd ${HOME}/bin

   ln -s ${HOME}/lib/htslib/bin/bgzip bgzip
   ln -s ${HOME}/lib/htslib/bin/htsfile htsfile
   ln -s ${HOME}/lib/htslib/bin/tabix tabix

   ln -s ${HOME}/lib/samtools/bin/ace2sam ace2sam
   ln -s ${HOME}/lib/samtools/bin/blast2sam.pl blast2sam.pl
   ln -s ${HOME}/lib/samtools/bin/bowtie2sam.pl bowtie2sam.pl
   ln -s ${HOME}/lib/samtools/bin/export2sam.pl export2sam.pl
   ln -s ${HOME}/lib/samtools/bin/interpolate_sam.pl interpolate_sam.pl
   ln -s ${HOME}/lib/samtools/bin/maq2sam-long maq2sam-long
   ln -s ${HOME}/lib/samtools/bin/maq2sam-short maq2sam-short
   ln -s ${HOME}/lib/samtools/bin/md5fa md5fa
   ln -s ${HOME}/lib/samtools/bin/md5sum-lite md5sum-lite
   ln -s ${HOME}/lib/samtools/bin/novo2sam.pl novo2sam.pl
   ln -s ${HOME}/lib/samtools/bin/plot-bamstats plot-bamstats
   ln -s ${HOME}/lib/samtools/bin/psl2sam.pl psl2sam.pl
   ln -s ${HOME}/lib/samtools/bin/sam2vcf.pl sam2vcf.pl
   ln -s ${HOME}/lib/samtools/bin/samtools samtools
   ln -s ${HOME}/lib/samtools/bin/samtools.pl samtools.pl
   ln -s ${HOME}/lib/samtools/bin/seq_cache_populate.pl seq_cache_populate.pl
   ln -s ${HOME}/lib/samtools/bin/soap2sam.pl soap2sam.pl
   ln -s ${HOME}/lib/samtools/bin/varfilter.py varfilter.py
   ln -s ${HOME}/lib/samtools/bin/wgsim wgsim
   ln -s ${HOME}/lib/samtools/bin/wgsim_eval.pl wgsim_eval.pl
   ln -s ${HOME}/lib/samtools/bin/zoom2sam.pl zoom2sam.pl


Python Modules
^^^^^^^^^^^^^^

- mg-tool-api
- mg-common
- pylint
- pytest
- macs2

Installation
------------
Directly from GitHub:

.. code-block:: none
   :linenos:

   git clone https://github.com/Multiscale-Genomics/mg-process-macs2.git

Using pip:

.. code-block:: none
   :linenos:

   pip install git+https://github.com/Multiscale-Genomics/mg-process-macs2.git