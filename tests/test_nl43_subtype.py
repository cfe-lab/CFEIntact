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
