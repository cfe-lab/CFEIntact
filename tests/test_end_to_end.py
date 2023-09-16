
import pytest
import requests
import filecmp
import os
import tarfile
import intact.intact as main

def run_end_to_end(tmp_path, data_file, expected_dir, output_csv):
    main.intact(
        working_dir=tmp_path,
        input_file=data_file,
        subtype="all",
        include_packaging_signal=True,
        include_rre=True,
        check_major_splice_donor_site=True,
        run_hypermut=True,
        check_long_deletion=True,
        check_nonhiv=True,
        check_scramble=True,
        check_internal_inversion=True,
        include_small_orfs=True,
        output_csv=output_csv,
        )

    result = filecmp.dircmp(tmp_path, expected_dir)
    assert result.left_list == result.right_list
    assert result.diff_files == []
    assert result.common_funny == []
    assert result.common == result.right_list


def test_single(tmp_path, request):
    pwd = request.fspath.dirname
    run_end_to_end(tmp_path,
                   os.path.join(pwd, "data-single.fasta"),
                   os.path.join(pwd, "expected-results-single"),
                   output_csv=False)


def test_single_csv(tmp_path, request):
    pwd = request.fspath.dirname
    run_end_to_end(tmp_path,
                   os.path.join(pwd, "data-single.fasta"),
                   os.path.join(pwd, "expected-results-single-csv"),
                   output_csv=True)


def test_small(tmp_path, request):
    pwd = request.fspath.dirname
    run_end_to_end(tmp_path,
                   os.path.join(pwd, "data-small.fasta"),
                   os.path.join(pwd, "expected-results-small"),
                   output_csv=False)


def test_small_csv(tmp_path, request):
    pwd = request.fspath.dirname
    run_end_to_end(tmp_path,
                   os.path.join(pwd, "data-small.fasta"),
                   os.path.join(pwd, "expected-results-small-csv"),
                   output_csv=True)


@pytest.mark.slow
def test_large(tmp_path, request):
    pwd = request.fspath.dirname
    run_end_to_end(tmp_path,
                   os.path.join(pwd, "data-large.fasta"),
                   os.path.join(pwd, "expected-results-large"),
                   output_csv=False)


@pytest.mark.slow
def test_large_csv(tmp_path, request):
    pwd = request.fspath.dirname
    run_end_to_end(tmp_path,
                   os.path.join(pwd, "data-large.fasta"),
                   os.path.join(pwd, "expected-results-large-csv"),
                   output_csv=True)


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
    run_end_to_end(tmp_path, file_path, expected_dir_path, output_csv=False)

