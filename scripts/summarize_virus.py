import logging
import pandas as pd
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


if "snakemake" in globals():
    log_fp = snakemake.log[0]  # type: ignore
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        handlers=[logging.FileHandler(log_fp), logging.StreamHandler()],
    )
    logger = logging.getLogger(__name__)
    try:
        logger.debug("Starting summary script")
        from scripts.map import tools_to_model
        from scripts.parse import parse_all_outputs, parse_tsv

        parsers = {
            "genomad_plasmid_summary": parse_tsv,
            "genomad_virus_summary": parse_tsv,
            "genomad_plasmid_genes": parse_tsv,
            "genomad_virus_genes": parse_tsv,
        }

        outputs: dict[str, list[Path]] = {
            "genomad_plasmid_summary": [Path(fp) for fp in snakemake.input if "plasmid_summary" in fp],  # type: ignore
            "genomad_virus_summary": [Path(fp) for fp in snakemake.input if "virus_summary" in fp],  # type: ignore
            "genomad_plasmid_genes": [Path(fp) for fp in snakemake.input if "plasmid_genes" in fp],  # type: ignore
            "genomad_virus_genes": [Path(fp) for fp in snakemake.input if "virus_genes" in fp],  # type: ignore
        }

        tool_reports = {Path(fp).stem: Path(fp) for fp in snakemake.output.tool_reports}  # type: ignore

        virus_summary_output = Path(snakemake.output.virus)  # type: ignore

        logger.debug("Parsing tool outputs", extra={"outputs": outputs})
        parsed_outputs = parse_all_outputs(outputs, parsers)

        logger.debug("Writing individual tool reports", extra={"tool_reports": tool_reports})
        for tool, df in parsed_outputs.items():
            df.to_csv(tool_reports[tool], sep="\t", index=False)
            logger.debug(
                "Wrote virus tool report",
                extra={"tool": tool, "dataframe": df.to_dict(orient="list")},
            )

        logger.debug(
            "Dummying final virus summary for now (TODO: Implement same as summarize_all.py)"
        )
        pd.DataFrame().to_csv(virus_summary_output, sep="\t", index=False)
        logger.debug("Finished writing final summaries")
    except Exception:
        logger.exception("Encountered error during virus summarization")
        raise
