import json
import os
import pytest
import cfeintact.intact as main
import cfeintact.subtypes as st


def test_nl43_is_known_subtype():
    assert "NL43" in st.subtypes()


def _run_nl43_self_check(tmp_path):
    nl43_fasta = os.path.join(
        os.path.dirname(__file__), "..",
        "src", "cfeintact", "subtype_alignments", "NL43.fasta",
    )
    main.check(
        output_dir=tmp_path,
        input_file=nl43_fasta,
        subtype="NL43",
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
        output_csv=False,
    )
    with open(os.path.join(tmp_path, "holistic.json")) as f:
        holistic = json.load(f)
    with open(os.path.join(tmp_path, "regions.json")) as f:
        regions = json.load(f)
    with open(os.path.join(tmp_path, "defects.json")) as f:
        defects = json.load(f)
    return holistic, regions, defects


def _get_orf(regions, name):
    for orf in regions["HIV_ENVGFP"]:
        if orf["region"] == name:
            return orf
    return None


def test_nl43_self_check_is_intact(tmp_path):
    holistic, regions, defects = _run_nl43_self_check(tmp_path)
    assert "HIV_ENVGFP" in holistic
    assert "HIV_ENVGFP" in regions
    assert holistic["HIV_ENVGFP"]["intact"] is True, \
        f"NL43 against self should be intact, got defects: {defects}"


def test_nl43_self_check_no_defects(tmp_path):
    _, _, defects = _run_nl43_self_check(tmp_path)
    hiv_key = "HIV_ENVGFP"
    assert hiv_key in defects
    assert defects[hiv_key] == [], \
        f"Expected no defects, got: {defects[hiv_key]}"


def test_nl43_all_orfs_have_zero_distance(tmp_path):
    _, regions, _ = _run_nl43_self_check(tmp_path)
    hiv_key = "HIV_ENVGFP"
    orfs = regions[hiv_key]
    assert len(orfs) > 0, "Expected at least one ORF"
    for orf in orfs:
        msg = f"{orf['region']}: distance={orf['distance']} != 0"
        assert orf["distance"] == 0.0, msg


# --- Coordinate tests against GTF annotations ---
# GTF coordinates from HIV_NL43_EnvGFP.gtf.txt (1-based).
# All assertions convert to 0-based to match the tool's output.

def test_nl43_gag_coordinates(tmp_path):
    _, regions, _ = _run_nl43_self_check(tmp_path)
    orf = _get_orf(regions, "gag")
    assert orf is not None
    assert orf["start"] == 789   # GTF 790 -> 0-based
    assert orf["end"] == 2291    # GTF 2292 -> 0-based


@pytest.mark.xfail(reason="HXB2 pol maps to NL43 at 2085 (identical sequence), "
                          "but GTF annotates NL43 pol start at 2358 due to "
                          "different gag-pol organization between strains")
def test_nl43_pol_coordinates(tmp_path):
    _, regions, _ = _run_nl43_self_check(tmp_path)
    orf = _get_orf(regions, "pol")
    assert orf is not None
    assert orf["start"] == 2357   # GTF 2358 -> 0-based
    assert orf["end"] == 5095    # GTF 5096 -> 0-based


def test_nl43_vif_coordinates(tmp_path):
    _, regions, _ = _run_nl43_self_check(tmp_path)
    orf = _get_orf(regions, "vif")
    assert orf is not None
    assert orf["start"] == 5040   # GTF 5041 -> 0-based
    assert orf["end"] == 5618    # GTF 5619 -> 0-based


def test_nl43_vpr_coordinates(tmp_path):
    _, regions, _ = _run_nl43_self_check(tmp_path)
    orf = _get_orf(regions, "vpr")
    assert orf is not None
    assert orf["start"] == 5558   # GTF 5559 -> 0-based
    assert orf["end"] == 5848    # GTF 5849 -> 0-based


def test_nl43_tat_exon1_coordinates(tmp_path):
    _, regions, _ = _run_nl43_self_check(tmp_path)
    orf = _get_orf(regions, "tat_exon1")
    assert orf is not None
    assert orf["start"] == 5829   # GTF 5830 -> 0-based
    assert orf["end"] == 6043    # GTF 6044 -> 0-based


def test_nl43_rev_exon1_coordinates(tmp_path):
    _, regions, _ = _run_nl43_self_check(tmp_path)
    orf = _get_orf(regions, "rev_exon1")
    assert orf is not None
    assert orf["start"] == 5968   # GTF 5969 -> 0-based
    assert orf["end"] == 6043    # GTF 6044 -> 0-based


def test_nl43_vpu_coordinates(tmp_path):
    _, regions, _ = _run_nl43_self_check(tmp_path)
    orf = _get_orf(regions, "vpu")
    assert orf is not None
    assert orf["start"] == 6060   # GTF 6061 -> 0-based
    assert orf["end"] == 6305    # GTF 6306 -> 0-based


def test_nl43_env_coordinates(tmp_path):
    _, regions, _ = _run_nl43_self_check(tmp_path)
    orf = _get_orf(regions, "env")
    assert orf is not None
    assert orf["start"] == 6220   # GTF envgfp 6221 -> 0-based
    assert orf["end"] == 9504    # GTF envgfp 9505 -> 0-based


@pytest.mark.xfail(reason="HXB2 tat_exon2 (8377-8469) maps to ~9091-9183 via "
                          "coordinate mapping, but the GTF annotates it at "
                          "9089-9134 (off by ~2-45 bases)")
def test_nl43_tat_exon2_coordinates(tmp_path):
    _, regions, _ = _run_nl43_self_check(tmp_path)
    orf = _get_orf(regions, "tat_exon2")
    assert orf is not None
    assert orf["start"] == 9088   # GTF 9089 -> 0-based
    assert orf["end"] == 9133    # GTF 9134 -> 0-based


@pytest.mark.xfail(reason="HXB2 rev_exon2 (8378-8653) maps to ~9092-9367 via "
                          "coordinate mapping, but the GTF annotates it at "
                          "9371-9563 (off by ~279 bases)")
def test_nl43_rev_exon2_coordinates(tmp_path):
    _, regions, _ = _run_nl43_self_check(tmp_path)
    orf = _get_orf(regions, "rev_exon2")
    assert orf is not None
    assert orf["start"] == 9370   # GTF 9371 -> 0-based
    assert orf["end"] == 9562    # GTF 9563 -> 0-based


def test_nl43_nef_coordinates(tmp_path):
    _, regions, _ = _run_nl43_self_check(tmp_path)
    orf = _get_orf(regions, "nef")
    assert orf is not None
    assert orf["start"] == 9506   # GTF 9507 -> 0-based
    assert orf["end"] == 10126   # GTF 10127 -> 0-based
