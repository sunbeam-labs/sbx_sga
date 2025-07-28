from scripts.genomad_f import (
    parse_filelines,
    write_plasmid_report,
    write_virus_report,
    write_report,
)


plasmid_reports = snakemake.input.plasmid_summary
plasmid_output = snakemake.output.plasmid_summary_report
write_report(plasmid_output, plasmid_reports, op_type="plasmid")

virus_reports = snakemake.input.virus_summary
virus_output = snakemake.output.virus_summary_report
write_report(virus_output, virus_reports, op_type="virus")
