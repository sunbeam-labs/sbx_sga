import csv
import os
import traceback
from typing import Callable, Iterable


def write_report(
    writer,
    report_reader: Iterable[list],
    sample_name: str,
    log: Callable[[str], None] = print,
):
    log(f"[concat_files_f] Writing rows for sample {sample_name}")
    for row in report_reader:
        row.insert(0, sample_name)
        writer.writerow(row)


def summarize_all(
    report_paths,
    out_fp,
    folder_suffix="",
    header=True,
    log: Callable[[str], None] = print,
):
    log(
        "[concat_files_f] Preparing to summarize "
        f"{len(report_paths)} reports to {out_fp} (header={header})"
    )
    first_non_empty = True
    header_first = []
    with open(out_fp, "w") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        for report_path in report_paths:
            log(f"[concat_files_f] Processing report {report_path}")
            if os.path.getsize(report_path) < 5:
                log(
                    f"[concat_files_f] Skipping {report_path} because file is nearly empty"
                )
                continue
            sample_name = os.path.basename(os.path.dirname(report_path))
            sample_name = sample_name.removesuffix(folder_suffix)
            with open(report_path, "r") as in_f:
                report_reader = csv.reader(in_f, delimiter="\t")
                if header:
                    header_line = next(report_reader)
                    header_line.insert(0, "SampleID")
                    if first_non_empty:
                        log(f"[concat_files_f] Setting header to {header_line}")
                        header_first = header_line
                        writer.writerow(header_first)
                    else:
                        if header_line != header_first:
                            log(
                                f"[concat_files_f] Header mismatch for sample {sample_name}: {header_line}"
                            )
                            log(f"[concat_files_f] Expected header: {header_first}")
                            raise ValueError(
                                f"Headers in sample {sample_name} doesn't match the first file"
                            )
                write_report(writer, report_reader, sample_name, log)
            first_non_empty = False
    log(f"[concat_files_f] Finished writing combined report to {out_fp}")


if "snakemake" in globals():
    log_fp = snakemake.log[0]  # type: ignore
    with open(log_fp, "w") as log_file:

        def log(message: str) -> None:
            log_file.write(f"[concat_files.py] {message}\n")
            log_file.flush()

        input_files = list(snakemake.input) # type: ignore
        output_path = snakemake.output[0] # type: ignore
        suffix = snakemake.params.suffix # type: ignore
        header = snakemake.params.header # type: ignore

        try:
            log(
                "Starting summarize_all with "
                f"{len(input_files)} files, suffix={suffix!r}, header={header} -> {output_path}"
            )

            summarize_all(
                input_files,
                output_path,
                suffix,
                header,
                log,
            )

            log(f"Finished summarizing files into {output_path}")
        except Exception as error:
            log(f"Encountered error: {error}")
            log(traceback.format_exc())
            raise
