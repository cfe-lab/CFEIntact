
import os
import sys
from unittest.mock import patch
import cfeintact.main as main

def test_basic_run(tmp_path, request):
    output_dir = os.path.join(tmp_path, "wd")
    pwd = request.fspath.dirname
    input = os.path.join(pwd, "data-small.fasta")

    testargs = ["cfeintact", "check", "--subtype", "B", "--output", output_dir, "--output-json", input]

    with patch.object(sys, 'argv', testargs):
        try: main.cli()
        except SystemExit as e: assert e.code == 0

    assert os.path.exists(os.path.join(output_dir, "defects.json"))


def test_version_run(tmp_path, request):
    testargs = ["cfeintact", "version"]

    with patch.object(sys, 'argv', testargs):
        try: main.cli()
        except SystemExit: pass
