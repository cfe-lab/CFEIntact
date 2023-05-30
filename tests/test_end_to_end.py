
import pytest
import requests
import filecmp
import os
import intact.intact as main

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


def test_small(tmp_path):
    run_end_to_end(tmp_path, "tests/data-small.fasta", "./tests/expected-results-small")


@pytest.mark.slow
def test_large(tmp_path):
    run_end_to_end(tmp_path, "tests/data-large.fasta", "./tests/expected-results-large")


@pytest.fixture(scope="session")
def huge_data_file(tmpdir_factory):
    tmp_dir = tmpdir_factory.mktemp('downloaded_files')
    tmp_tar_path = tmp_dir.join('test_file.tar')
    # response = requests.get("https://github.com/cfe-lab/HIVIntact/releases/download/test-data/huge-data-file.tar")
    # with open(tmp_path, 'wb') as f:
        # f.write(response.content)

    # return (tmp_path, "tmp1")
    return ("/home/doname/my/wd/hiv_sequences/all-psd-subtype-b-long.fasta", "haha")


@pytest.mark.slow
@pytest.mark.overnight
def test_huge(tmp_path, huge_data_file):
    file_path, expected_dir_path = huge_data_file
    run_end_to_end(tmp_path, file_path, expected_dir_path)

