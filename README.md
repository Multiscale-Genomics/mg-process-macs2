# mg-process-macs2

[![Documentation Status](https://readthedocs.org/projects/mg-process-macs2/badge/?version=latest)](https://mg-process-macs2.readthedocs.io/en/latest/?badge=latest) [![Build Status](https://travis-ci.org/Multiscale-Genomics/mg-process-macs2.svg?branch=master)](https://travis-ci.org/Multiscale-Genomics/mg-process-macs2) [![Code Health](https://landscape.io/github/Multiscale-Genomics/mg-process-macs2/master/landscape.svg?style=flat)](https://landscape.io/github/Multiscale-Genomics/mg-process-macs2/master)

Example pipelines file that is ready to run in the VRE matching the code in the HowTo documentation.

This repo structure workflows and tools can be forked and used as the base template for new tools and workflows. It should have all of the base functionality and is set up for unit macs2ing and with pylint to ensure code clarity.

# Requirements
- pyenv and pyenv-virtualenv
- Python 2.7.12
- Python Modules:
  - pylint
  - pytest
  - mg-tool-api

Installation
------------

Directly from GitHub:

```
cd ${HOME}/code

git clone https://github.com/Multiscale-Genomics/mg-process-macs2.git

cd mg-process-macs2
```

Create the Python environment

```
pyenv-virtualenv 2.7.12 mg-process-macs2
pyenv activate mg-process-macs2
pip install -e .
pip install -r requirements.txt
```
