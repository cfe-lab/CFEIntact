#! /bin/sh

python3 -m pip install .[dev]
python3 -m pip install .[test]
python3 -m pip uninstall -y cfeintact
