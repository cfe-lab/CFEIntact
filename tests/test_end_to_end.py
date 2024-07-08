
import pytest
import requests
import filecmp
import os
import tarfile
import cfeintact.intact as main


def check_outputs(tmp_path, expected_dir):
    result = filecmp.dircmp(tmp_path, expected_dir)
    assert result.left_list == result.right_list
    assert result.diff_files == []
    assert result.common_funny == []
    assert result.common == result.right_list


def run_end_to_end(tmp_path, data_file, expected_dir, subtype, output_csv):
    main.check(
        working_dir=tmp_path,
        input_file=data_file,
        subtype=subtype,
        check_packaging_signal=True,
        check_rre=True,
        check_major_splice_donor_site=True,
        check_hypermut=True,
        check_long_deletion=True,
        check_nonhiv=True,
        check_scramble=True,
        check_internal_inversion=True,
        check_unknown_nucleotides=True,
        check_small_orfs=True,
        check_distance=True,
        output_csv=output_csv,
    )

    check_outputs(tmp_path, expected_dir)


def test_single(tmp_path, request):
    pwd = request.fspath.dirname
    run_end_to_end(tmp_path,
                   os.path.join(pwd, "data-single.fasta"),
                   os.path.join(pwd, "expected-results-single"),
                   subtype="all",
                   output_csv=False)


def test_single_csv(tmp_path, request):
    pwd = request.fspath.dirname
    run_end_to_end(tmp_path,
                   os.path.join(pwd, "data-single.fasta"),
                   os.path.join(pwd, "expected-results-single-csv"),
                   subtype="all",
                   output_csv=True)


def test_single_hxb2(tmp_path, request):
    pwd = request.fspath.dirname
    run_end_to_end(tmp_path,
                   os.path.join(pwd, "data-single.fasta"),
                   os.path.join(pwd, "expected-results-single-hxb2"),
                   subtype="HXB2",
                   output_csv=False)


def test_single_noblast(tmp_path, request):
    pwd = request.fspath.dirname
    data_file = os.path.join(pwd, "data-single.fasta")
    expected_dir = os.path.join(pwd, "expected-results-single-noblast")
    subtype = "HXB2"
    output_csv = True

    main.check(
        working_dir=tmp_path,
        input_file=data_file,
        subtype=subtype,
        check_packaging_signal=True,
        check_rre=True,
        check_major_splice_donor_site=True,
        check_hypermut=True,
        check_long_deletion=True,
        check_nonhiv=False,
        check_scramble=False,
        check_internal_inversion=False,
        check_unknown_nucleotides=True,
        check_small_orfs=True,
        check_distance=True,
        output_csv=output_csv,
    )

    check_outputs(tmp_path, expected_dir)


def test_small(tmp_path, request):
    pwd = request.fspath.dirname
    run_end_to_end(tmp_path,
                   os.path.join(pwd, "data-small.fasta"),
                   os.path.join(pwd, "expected-results-small"),
                   subtype="all",
                   output_csv=False)


def test_small_csv(tmp_path, request):
    pwd = request.fspath.dirname
    run_end_to_end(tmp_path,
                   os.path.join(pwd, "data-small.fasta"),
                   os.path.join(pwd, "expected-results-small-csv"),
                   subtype="all",
                   output_csv=True)


@pytest.mark.slow
def test_large(tmp_path, request):
    pwd = request.fspath.dirname
    run_end_to_end(tmp_path,
                   os.path.join(pwd, "data-large.fasta"),
                   os.path.join(pwd, "expected-results-large"),
                   subtype="all",
                   output_csv=False)


@pytest.mark.slow
def test_large_csv(tmp_path, request):
    pwd = request.fspath.dirname
    run_end_to_end(tmp_path,
                   os.path.join(pwd, "data-large.fasta"),
                   os.path.join(pwd, "expected-results-large-csv"),
                   subtype="all",
                   output_csv=True)


@pytest.mark.slow
def test_large_hxb2(tmp_path, request):
    pwd = request.fspath.dirname
    run_end_to_end(tmp_path,
                   os.path.join(pwd, "data-large.fasta"),
                   os.path.join(pwd, "expected-results-large-hxb2"),
                   subtype="HXB2",
                   output_csv=False)


@pytest.fixture(scope="session")
def huge_data_file(tmpdir_factory):
    tmp_dir = tmpdir_factory.mktemp("downloaded_files")
    tmp_tar_path = tmp_dir.join("test_file.tar.gz")
    data_file = tmp_dir.join("cfeintact_data", "all-psd-subtype-b-long.fasta")
    expected_dir = tmp_dir.join("cfeintact_data", "expected_output")

    response = requests.get("https://github.com/cfe-lab/CFEIntact/releases/download/test-data/cfeintact_data.tar.gz")
    with open(tmp_tar_path, "wb") as f:
        f.write(response.content)

    with tarfile.open(tmp_tar_path, "r:gz") as tar:
        tar.extractall(path=tmp_dir)

    return (data_file, expected_dir)


@pytest.mark.slow
@pytest.mark.overnight
def test_huge(tmp_path, huge_data_file):
    file_path, expected_dir_path = huge_data_file
    run_end_to_end(tmp_path, file_path, expected_dir_path, subtype="all", output_csv=False)


def test_edgy(tmp_path, request):
    pwd = request.fspath.dirname
    run_end_to_end(tmp_path,
                   os.path.join(pwd, "data-edgy.fasta"),
                   os.path.join(pwd, "expected-results-edgy"),
                   subtype="all",
                   output_csv=False)
