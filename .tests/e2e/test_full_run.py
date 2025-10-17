import os
import pytest
import shutil
import subprocess as sp
from pathlib import Path


@pytest.fixture
def setup(tmp_path: Path):
    reads_fp = Path(".tests/data/reads/").resolve()
    snippy_fp = Path(".tests/data/Bfragilis.fasta").resolve()
    checkm_fp = Path(".tests/data/checkm_test.dmnd").resolve()
    mash_fp = Path(".tests/data/mash_test.msh").resolve()

    project_dir = tmp_path / "project/"

    sp.check_output(["sunbeam", "init", "--data_fp", reads_fp, project_dir])

    config_fp = project_dir / "sunbeam_config.yml"

    config_str = f"sbx_sga: {{mash_ref: '{mash_fp}'}}"
    sp.check_output(
        [
            "sunbeam",
            "config",
            "--modify",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    config_str = f"sbx_sga: {{checkm_ref: '{checkm_fp}'}}"
    sp.check_output(
        [
            "sunbeam",
            "config",
            "--modify",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    config_str = f"sbx_sga: {{bakta_ref: '{tmp_path}/bakta/db/'}}"
    Path(tmp_path / "bakta/db/").mkdir(parents=True, exist_ok=True)

    sp.check_output(
        [
            "sunbeam",
            "config",
            "--modify",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    genomad_fp = tmp_path / "genomad_db"
    config_str = f"sbx_sga: {{genomad_ref: '{genomad_fp.parent}'}}"
    genomad_fp.mkdir(parents=True, exist_ok=True)
    (genomad_fp / "version.txt").touch()

    sp.check_output(
        [
            "sunbeam",
            "config",
            "--modify",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    config_str = f"sbx_sga: {{snippy_ref: '{snippy_fp}'}}"
    sp.check_output(
        [
            "sunbeam",
            "config",
            "--modify",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    config_str = f"sbx_sga: {{sylph_ref: '{tmp_path}/gtdb-r220-c200-dbv1.syldb'}}"
    (tmp_path / "gtdb-r220-c200-dbv1.syldb").touch()
    sp.check_output(
        [
            "sunbeam",
            "config",
            "--modify",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    yield tmp_path, project_dir


@pytest.fixture
def run_sunbeam(setup):
    tmp_path, project_dir = setup
    output_fp = project_dir / "sunbeam_output"
    log_fp = output_fp / "logs"
    stats_fp = project_dir / "stats"

    sp.run(
        [
            "sunbeam",
            "run",
            "--profile",
            str(project_dir),
            "--skip",
            "decontam",
            "--exclude",
            "all",
            "-n",
        ],
        check=True,
    )

    sbx_proc = sp.run(
        [
            "sunbeam",
            "run",
            "--profile",
            str(project_dir),
            "clean_sga",
            "all_sga_virus",
            "all_sga_snippy",
            "--directory",
            str(tmp_path),
            "--include",
            "sbx_sga",
            "--show-failed-logs",
            "-n",
        ],
        capture_output=True,
        text=True,
    )

    print("STDOUT: ", sbx_proc.stdout)
    print("STDERR: ", sbx_proc.stderr)

    if os.getenv("GITHUB_ACTIONS") == "true":
        try:
            shutil.copytree(log_fp, "logs/")
            shutil.copytree(stats_fp, "stats/")
        except FileNotFoundError:
            print("No logs or stats directory found.")

    output_fp = project_dir / "sunbeam_output"
    benchmarks_fp = project_dir / "stats/"

    yield output_fp, benchmarks_fp, sbx_proc


def test_dry_run(run_sunbeam):
    output_fp, benchmarks_fp, proc = run_sunbeam

    assert proc.returncode == 0, f"Sunbeam dryrun failed with error: {proc.stderr}"
