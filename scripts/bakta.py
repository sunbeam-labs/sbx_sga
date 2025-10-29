import csv
import traceback
from typing import Callable

def parse_file(filelines, log: Callable[[str], None] = print):
    log(f"[bakta_f] Parsing {len(filelines)} raw lines from Bakta report")
    parsed_dict = {}
    if len(filelines) != 0:
        for line in filelines:
            line = line.rstrip().split(":")
            try:
                key = line[0]
                value = line[1].strip()
            except Exception:
                continue
            parsed_dict[key] = value
            log(f"[bakta_f] Captured entry key={key!r} value={value!r}")
    return parsed_dict


def filter_keys(parsed_dict, log: Callable[[str], None] = print):
    filtered = {key: value for key, value in parsed_dict.items() if value != ""}
    log(
        "[bakta_f] Filtered parsed entries: "
        f"kept {len(filtered)} of {len(parsed_dict)} keys with non-empty values"
    )
    return filtered


def write_to_report(report_fp, output_fp, log: Callable[[str], None] = print):
    log(f"[bakta_f] Opening Bakta report at {report_fp}")
    with open(report_fp, "r") as f_in:
        lines = f_in.readlines()

    parsed_dict = parse_file(lines, log)
    filtered = filter_keys(parsed_dict, log)

    with open(output_fp, "w") as op:
        writer = csv.writer(op, delimiter="\t")
        writer.writerow(filtered.keys())
        writer.writerow(filtered.values())
    log(f"[bakta_f] Wrote Bakta summary with {len(filtered)} keys to {output_fp}")


if "snakemake" in globals():
    log_fp = snakemake.log[0] # type: ignore
    with open(log_fp, "w") as log_file:

        def log(message: str) -> None:
            log_file.write(f"[bakta.py] {message}\n")
            log_file.flush()

        report = snakemake.input[0] # type: ignore
        output = snakemake.output[0] # type: ignore

        try:
            log(f"Starting write_to_report with report={report} -> output={output}")
            write_to_report(report, output, log)
            log(f"Finished writing parsed Bakta summary to {output}")
        except Exception as error:
            log(f"Encountered error: {error}")
            log(traceback.format_exc())
            raise
