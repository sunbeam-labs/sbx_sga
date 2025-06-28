import sys
from scripts.mash_f import open_report, parse_report, contamination_call, write_report


sorted_report = snakemake.input[0]
output = snakemake.output[0]

sample, filelines = open_report(sorted_report)

if len(filelines) == 0:
    empty_dict = {"":""}
    write_report(output, sample, empty_dict)
else:
    parsed_report = parse_report(filelines)
    mash_dict = contamination_call(parsed_report)
    write_report(output, sample, mash_dict)
