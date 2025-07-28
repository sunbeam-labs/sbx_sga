import os
from pathlib import Path


def parse_filelines(lines, sample_name):
    parsed_lines = []
    for line in lines[1:]:
        if not line.strip():
            continue
        line = line.rstrip().split("\t")
        line.insert(0, sample_name)
        parsed_lines.append("\t".join(line))
    return parsed_lines


def write_virus_report(output_file, lines, write_header=False):
    mode = "w" if write_header else "a"
    with open(output_file, mode) as virus_op:
        if write_header:
            virus_op.write(
                "genome\tgene\tstart\tend\tlength\tstrand\tgc_content\tgenetic_code\trbs_motif\tmarkerevalue\tbitscore\tuscg\tplasmid_hallmark\tvirus_hallmark\ttaxid\ttaxname\tannotation_conjscan\tannotation_amr\tannotation_accessions\tannotation_description\n"
            )
        for line in lines:
            virus_op.write(f"{line}\n")


def write_plasmid_report(output_file, lines, write_header=False):
    mode = "w" if write_header else "a"
    with open(output_file, mode) as plasmid_op:
        if write_header:
            plasmid_op.write(
                "gene\tstart\tend\tlength\tstrand\tgc_content\tgenetic_code\trbs_motif\tmarkerevalue\tbitscore\tuscg\tplasmid_hallmark\tvirus_hallmark\ttaxid\ttaxname\tannotation_conjscan\tannotation_amr\tannotation_accessions\tannotation_description\n"
            )
        for line in lines:
            plasmid_op.write(f"{line}\n")


def write_report(op_file, report_paths, op_type):
    first = True
    for report_path in report_paths:
        sample_path = Path(report_path)
        sample_name = sample_path.parts[-3]

        with open(report_path, "r") as opened_f:
            report_lines = opened_f.readlines()
            formatted_lines = parse_filelines(report_lines, sample_name)

            if op_type == "virus":
                write_virus_report(op_file, formatted_lines, write_header=first)
            elif op_type == "plasmid":
                write_plasmid_report(op_file, formatted_lines, write_header=first)

        first = False
