import importlib.metadata
import json
import logging
import sys

DIRECT_DEPENDENCIES = [
    "moalmanac",
    "flask",
    "Frozen-Flask",
    "ipykernel",
    "jupyter",
    "matplotlib",
    "numpy",
    "pandas",
    "requests",
    "scikit-learn",
    "scipy",
]


def collect_environment_metadata():
    packages = {}
    for pkg in DIRECT_DEPENDENCIES:
        try:
            packages[pkg] = importlib.metadata.version(pkg)
        except importlib.metadata.PackageNotFoundError:
            packages[pkg] = "not installed"
    return {"python": sys.version, "packages": packages}


class Logger:
    @staticmethod
    def generate_filename(folder, label):
        return f"{folder}/{label}.log"

    @classmethod
    def setup(cls, output_folder, file_prefix, config):
        filename = cls.generate_filename(folder=output_folder, label=file_prefix)
        logger = logging.getLogger()
        if logger.hasHandlers():
            logger.handlers.clear()
        logging.basicConfig(
            filename=filename,
            filemode="w",
            format="%(asctime)s - %(levelname)s - %(message)s",
            level=config["logging"]["level"],
        )

    @classmethod
    def shutdown(cls):
        logging.shutdown()


class Messages:
    @staticmethod
    def dataframe_size(label, dataframe, add_line_break=False):
        message = f"{label}: {dataframe.shape[0]} rows and {dataframe.shape[1]} columns"
        if add_line_break:
            message = f"{message}\n"
        logging.info(message)

    @staticmethod
    def general(message, add_line_break=False):
        if add_line_break:
            message = f"{message}\n"
        logging.info(message)

    @staticmethod
    def header(label):
        message = f"--- {label} ---"
        logging.info(message)

    @staticmethod
    def inputs(label, dictionary):
        message = f"Input data, {label}:\n%s"
        logging.info(message, json.dumps(dictionary, indent=4))

    @staticmethod
    def environment(metadata):
        logging.info("Environment metadata:\n%s", json.dumps(metadata, indent=4))

    @staticmethod
    def start():
        message = "Starting execution of Molecular Oncology Almanac"
        logging.info(message)
