import argparse
import datetime
import json
import pathlib
import shutil
import time

from reader import Ini

import moalmanac

# ruff: noqa: E501
METADATA = {
    "patient_id": "example",
    "reported_tumor_type": "MEL",
    "stage": "Metastatic",
    "description": "Test patient for development runs",
    "purity": 0.85,
    "ploidy": 4.02,
    "WGD": True,
    "microsatellite_status": "msih",
    "tumor_mutational_burden": 10.1,
}

INPUTS_EMPTY = {
    "snv_handle": "",
    "indel_handle": "",
    "bases_covered_handle": "",
    "called_cn_handle": "",
    "cnv_handle": "",
    "fusion_handle": "",
    "germline_handle": "",
    "validation_handle": "",
    "mutational_signatures_path": "",
    "disable_matchmaking": False,
}

INPUTS_FULL = {
    "snv_handle": "../example_data/example_patient.capture.somatic.snvs.maf",
    "indel_handle": "../example_data/example_patient.capture.somatic.indels.maf",
    "bases_covered_handle": "../example_data/example_patient.capture.somatic.coverage.txt",
    "called_cn_handle": "../example_data/example_patient.capture.somatic.called.cna.txt",
    "cnv_handle": "../example_data/example_patient.capture.somatic.seg.annotated",
    "fusion_handle": "../example_data/example_patient.rna.star.fusions.txt",
    "germline_handle": "../example_data/example_patient.capture.germline.maf",
    "validation_handle": "../example_data/example_patient.rna.somatic.snvs.maf",
    "mutational_signatures_path": "../example_data/example_patient.capture.sbs_contributions.txt",
    "disable_matchmaking": False,
}

# (name, inputs, function toggle value)
RUNS = [
    ("full-all-enabled", INPUTS_FULL, "on"),
    ("empty-all-enabled", INPUTS_EMPTY, "on"),
    ("full-all-disabled", INPUTS_FULL, "off"),
    ("empty-all-disabled", INPUTS_EMPTY, "off"),
]


def read_ini_files():
    """
    Reads the configuration and datasource path files from the working directory.

    Returns:
        tuple[configparser.ConfigParser, dict, dict]: The runtime configuration, the annotation datasource paths, and the preclinical datasource paths.
    """
    config_ini = Ini.read(
        "config.ini",
        extended_interpolation=False,
        convert_to_dictionary=False,
    )
    db_paths = Ini.read(
        "annotation-databases.ini",
        extended_interpolation=True,
        convert_to_dictionary=True,
        resolve_paths=True,
    )["paths"]
    preclinical_db_paths = Ini.read(
        "preclinical-databases.ini",
        extended_interpolation=True,
        convert_to_dictionary=True,
        resolve_paths=True,
    )["paths"]
    return config_ini, db_paths, preclinical_db_paths


def main(label, output_root, date, overwrite=False):
    """
    Runs the four example executions and writes a manifest of their output folders.
    The runs are full and empty inputs, each with all function toggles on and off.
    Output folder paths in the manifest are relative to the manifest's folder.

    Args:
        label (str): Either "pre" or "post", appended to each output folder name.
        output_root (str): Folder in which output folders and the manifest are created.
        date (str): ISO 8601 date prefixed to each output folder and the manifest name.
        overwrite (bool): Delete existing, non-empty output folders before running
            (default: False).

    Raises:
        SystemExit: If an output folder already contains files and overwrite is False,
            since files left from an earlier run would hide missing outputs.
    """
    output_root = pathlib.Path(output_root)
    folder_names = {
        name: f"{date}-example-outputs-{name}-{label}" for name, _, _ in RUNS
    }

    existing = [
        output_root / folder
        for folder in folder_names.values()
        if (output_root / folder).is_dir() and any((output_root / folder).iterdir())
    ]
    if existing and not overwrite:
        listed = "\n  ".join(str(folder) for folder in existing)
        raise SystemExit(
            f"Output folders already contain files:\n  {listed}\n"
            "Remove them, choose another --date or --output-root, or pass --overwrite.",
        )
    for folder in existing:
        shutil.rmtree(folder)

    config_ini, db_paths, preclinical_db_paths = read_ini_files()
    manifest = {"label": label, "date": date, "runs": {}}

    for name, inputs, toggle in RUNS:
        for key, _ in config_ini.items("function_toggle"):
            config_ini["function_toggle"][key] = toggle

        output_directory = output_root / folder_names[name]
        output_directory.mkdir(parents=True, exist_ok=True)

        start_time = time.time()
        moalmanac.main(
            patient=METADATA,
            inputs=inputs,
            output_folder=str(output_directory),
            config=config_ini,
            dbs=db_paths,
            dbs_preclinical=preclinical_db_paths,
        )
        elapsed = round(time.time() - start_time, 4)
        n_files = len(list(output_directory.iterdir()))
        print(f"[{name}] {elapsed} seconds, {n_files} files -> {output_directory}")
        manifest["runs"][name] = folder_names[name]

    manifest_path = output_root / f"{date}-regression-{label}.json"
    manifest_path.write_text(json.dumps(manifest, indent=4) + "\n")
    print(f"Wrote manifest {manifest_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Molecular Oncology Almanac regression runs",
        description="Run the four example executions used for output regression tests",
    )
    parser.add_argument(
        "--label",
        choices=["pre", "post"],
        default="pre",
    )
    parser.add_argument(
        "--output-root",
        default=".",
    )
    parser.add_argument(
        "--date",
        default=datetime.date.today().isoformat(),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Delete existing output folders for this date and label before running",
    )
    args = parser.parse_args()
    main(
        args.label,
        args.output_root,
        args.date,
        args.overwrite,
    )
