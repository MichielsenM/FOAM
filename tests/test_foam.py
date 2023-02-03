import tomli
from pathlib import Path
import foam

def test_versions_are_in_sync():
    """Checks if the pyproject.toml and package.__init__.py __version__ are in sync."""

    path = Path(__file__).resolve().parents[1] / "pyproject.toml"
    with open(path, mode="rb") as toml_file:
        pyproject = tomli.load(toml_file)
        pyproject_version = pyproject["tool"]["poetry"]["version"]

    package_init_version = foam.__version__
    assert package_init_version == pyproject_version

test_versions_are_in_sync()
