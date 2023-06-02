
import pytest
import requests
import filecmp
import os
import tarfile
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


def test_small(tmp_path, request):
    pwd = request.fspath.dirname
    run_end_to_end(tmp_path,
                   os.path.join(pwd, "data-small.fasta"),
                   os.path.join(pwd, "expected-results-small"))


@pytest.mark.slow
def test_large(tmp_path, request):
    pwd = request.fspath.dirname
    run_end_to_end(tmp_path,
                   os.path.join(pwd, "data-large.fasta"),
                   os.path.join(pwd, "expected-results-large"))


@pytest.fixture(scope="session")
def huge_data_file(tmpdir_factory):
    tmp_dir = tmpdir_factory.mktemp("downloaded_files")
    tmp_tar_path = tmp_dir.join("test_file.tar.gz")
    data_file = tmp_dir.join("hivintact_data", "all-psd-subtype-b-long.fasta")
    expected_dir = tmp_dir.join("hivintact_data", "expected_output")

    response = requests.get("https://github.com/cfe-lab/HIVIntact/releases/download/test-data/hivintact_data.tar.gz")
    with open(tmp_tar_path, "wb") as f:
        f.write(response.content)

    with tarfile.open(tmp_tar_path, "r:gz") as tar:
        tar.extractall(path=tmp_dir)

    return (data_file, expected_dir)


@pytest.mark.slow
@pytest.mark.overnight
def test_huge(tmp_path, huge_data_file):
    file_path, expected_dir_path = huge_data_file
    run_end_to_end(tmp_path, file_path, expected_dir_path)

