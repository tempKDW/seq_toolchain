#!/bin/bash
python -m venv ~/.virtualenv/test
source ~/.virtualenv/test/bin/activate
pip install -U pip
pip install biopython