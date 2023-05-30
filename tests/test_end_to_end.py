
import pytest
import requests
import filecmp
import os
import intact.intact as main


@pytest.fixture(scope="session")
def small_data_file(tmpdir_factory):
    return "tests/small.fasta"


@pytest.fixture(scope="session")
def large_data_file(tmpdir_factory):
    return "tests/data.fasta"


def run_end_to_end(tmp_path, data_file, expected_dir):
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

    result = filecmp.dircmp(tmp_path, expected_dir)
    assert result.left_list == result.right_list
    assert result.diff_files == []
    assert result.common_funny == []
    assert result.common == result.right_list


def test_small(tmp_path, small_data_file):
    run_end_to_end(tmp_path, small_data_file, "./tests/expected-results-small")


@pytest.mark.slow
def test_large(tmp_path, large_data_file):
    run_end_to_end(tmp_path, large_data_file, "./tests/expected-results-large")

