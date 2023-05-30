
import pytest
import requests
import filecmp
import os
import intact.intact as main

@pytest.fixture(scope="session")
def data_file(tmpdir_factory):
    return "tests/small.fasta"

def test_end_to_end(tmp_path, data_file):
    main.intact(
        working_dir=tmp_path,
        input_file=data_file,
        subtype="B",
        include_packaging_signal=True,
        include_rre=True,
        check_major_splice_donor_site=True,
        run_hypermut=True,
        check_long_deletion=True,
        include_small_orfs=False,
        )

    result = filecmp.dircmp(tmp_path, "./tests/expected-results-small")
    assert result.left_list == result.right_list
    assert result.diff_files == []
    assert result.common_funny == []
    assert result.common == result.right_list
