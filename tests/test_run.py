
import os
import sys
from unittest.mock import patch
import cfeintact.main as main

def test_basic_run(tmp_path, request):
    working_dir = os.path.join(tmp_path, "wd")
    pwd = request.fspath.dirname
    input = os.path.join(pwd, "data-small.fasta")

    testargs = ["cfeintact", "check", "--subtype", "B", "--working-folder", working_dir, "--output-json", input]

    with patch.object(sys, 'argv', testargs):
        try: main.cli()
        except SystemExit as e: assert e.code == 0

    assert os.path.exists(os.path.join(working_dir, "errors.json"))
