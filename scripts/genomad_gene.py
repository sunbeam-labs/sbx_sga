from scripts.genomad_gene_f import (
    parse_filelines,
    write_plasmid_report,
    write_virus_report,
    write_report,
)

plasmid_reports = snakemake.input.plasmid_gene_summary
plasmid_output = snakemake.output.plasmid_gene_report
write_report(plasmid_output, plasmid_reports, op_type="plasmid")

virus_reports = snakemake.input.virus_gene_summary
virus_output = snakemake.output.virus_gene_report
write_report(virus_output, virus_reports, op_type="virus")
