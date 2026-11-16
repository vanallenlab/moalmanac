import pathlib

import pytest
import reader

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
MOALMANAC_DIR = REPO_ROOT / "moalmanac"


@pytest.fixture(scope="session")
def repo_root():
    return REPO_ROOT


@pytest.fixture(scope="session")
def datasources_dir():
    return REPO_ROOT / "datasources"


@pytest.fixture(scope="session", autouse=True)
def _chdir_to_moalmanac():
    """Run tests from `moalmanac/`, since the package resolves paths relative to it."""
    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.chdir(MOALMANAC_DIR)
        yield


@pytest.fixture(scope="session")
def config_ini():
    return reader.Ini.read(MOALMANAC_DIR / "config.ini")
