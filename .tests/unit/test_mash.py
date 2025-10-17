import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from mash_f import (
    open_report,
    process_mash_line,
    get_first_non_phage_hit,
    parse_report,
    contamination_call,
    write_report,
)

test_data_dir = os.path.join(os.path.dirname(__file__), "test_data")


def test_marc_sample(tmp_path):
    sample_report = os.path.join(
        test_data_dir, "marc.bacteremia.2902_sorted_winning.tab"
    )
    sample_name, top_lines = open_report(sample_report)
    parsed_report = parse_report(top_lines)
    mash_dict = contamination_call(parsed_report)
    report = write_report(tmp_path / "marc_test.out", sample_name, mash_dict)


def test_open_report():
    sample_report = os.path.join(test_data_dir, "sample_sorted_winning.tab")
    sample_name, top_lines = open_report(sample_report)
    assert sample_name == "sample"
    assert len(top_lines) == 20


def test_process_mash_line():
    example_line = "0.998502\t969/1000\t111\t0\tGCF_000781075.1_ASM78107v1_genomic.fna.gz\t[54 seqs] NZ_JSOE01000001.1 Escherichia coli strain upec-121 upec-121_ctg_315, whole genome shotgun sequence [...]"
    species, median_multiplicity, identity, hits = process_mash_line(example_line)
    assert species == "Escherichia coli"
    assert median_multiplicity == 111.0
    assert identity == 0.998502
    assert hits == 969


def test_parse_report():
    sample_report = os.path.join(test_data_dir, "sample_sorted_winning.tab")
    with open(sample_report, "r") as f:
        top_lines = f.readlines()[:20]
    target_set = parse_report(top_lines)
    assert isinstance(target_set, set)
    assert len(target_set) > 0


def test_get_first_non_phage_hit():
    sample_report = os.path.join(test_data_dir, "sample_sorted_winning.tab")
    with open(sample_report, "r") as f:
        lines = f.readlines()[:3]
    phage_line = lines[0].rstrip() + " phage\n"
    test_lines = [phage_line] + lines[1:]
    result, idx = get_first_non_phage_hit(test_lines)
    expected = process_mash_line(lines[1])
    assert result == expected
    assert idx == 1


def test_contamination_call():
    example_dict = {"Escherichia coli", "Proteus penneri"}
    mash_dict = contamination_call(example_dict)
    assert isinstance(mash_dict, dict)
    assert mash_dict == {"Contaminated": "Escherichia coli Proteus penneri"}


def test_write_report(tmp_path):
    sample_output = tmp_path / "sample_mash.out"
    sample_name = "sample"
    mash_dict = {"Contaminated": "Escherichia coli Proteus penneri"}
    report = write_report(sample_output, sample_name, mash_dict)

    with open(report, "r") as f:
        content = f.read()
    assert "sample\tContaminated\tEscherichia coli Proteus penneri" in content
