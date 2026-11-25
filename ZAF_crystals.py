"""CIF-backed crystal models and a small, per-user crystal library for ZAF.

The original CIF is retained so that importing a structure never depends on
round-tripping through a reduced representation of its atoms or symmetry.
Pymatgen is imported only when a CIF is actually parsed.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import sys
import tempfile
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import TypeAlias
from uuid import uuid4


Matrix3: TypeAlias = tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]
Site: TypeAlias = tuple[float, float, float, float]


@dataclass(frozen=True)
class CrystalModel:
    """A full CIF cell in ZAF's normalized, three-index coordinate system.

    ``direct_matrix`` has Cartesian components in its rows and the a, b, c
    cell vectors in its columns. All lengths are divided by the cell's a
    length, while ``a_nm`` retains the physical scale. A site is fractional
    x, y, z followed by its occupancy-weighted atomic-number sum.
    """

    name: str
    crystal_system: str
    is_hexagonal: bool
    a_nm: float
    direct_matrix: Matrix3
    sites: tuple[Site, ...]
    symmetry_rotations: tuple[Matrix3, ...]
    cif_path: Path

    @property
    def lattice_system(self) -> str:
        """Alias matching CrysDiS terminology."""
        return self.crystal_system


def _clean_name(value: str) -> str:
    """Keep display names readable and harmless in labels and library keys."""
    name = re.sub(r"[\\/\x00-\x1f\x7f]+", "_", str(value)).strip(" ._")
    name = re.sub(r"\s+", " ", name)
    return name[:100].rstrip(" ._") or "Custom crystal"


def _matrix_tuple(matrix: object, *, integral: bool = False) -> Matrix3:
    rows = []
    for i in range(3):
        row = []
        for j in range(3):
            value = float(matrix[i][j])  # type: ignore[index]
            if integral and abs(value - round(value)) < 1e-8:
                value = float(round(value))
            elif abs(value) < 1e-12:
                value = 0.0
            row.append(value)
        rows.append(tuple(row))
    return tuple(rows)  # type: ignore[return-value]


def from_cif(path: str | Path, name: str | None = None) -> CrystalModel:
    """Parse a CIF's full input cell without primitive-cell reindexing.

    The conventional-cell indices in the CIF are therefore the indices ZAF
    shows to the user. A symmetry analysis provides candidate-axis rotations;
    if it is unavailable, the identity still permits a correct unsymmetrized
    search.
    """
    cif_path = Path(path).expanduser().resolve()
    if not cif_path.is_file():
        raise ValueError(f"CIF file does not exist: {cif_path}")
    try:
        from pymatgen.io.cif import CifParser
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    except ImportError as exc:
        raise ValueError("pymatgen is required to import CIF files") from exc

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            structure = CifParser(str(cif_path)).parse_structures(primitive=False)[0]
    except (IndexError, OSError, ValueError, TypeError, KeyError, RuntimeError) as exc:
        raise ValueError(f"Could not parse CIF file {cif_path.name}: {exc}") from exc

    lattice = structure.lattice
    a_angstrom = float(lattice.a)
    if not 0.0 < a_angstrom < float("inf"):
        raise ValueError("CIF has an invalid a lattice parameter")
    # Pymatgen stores cell vectors as rows; ZAF's basis convention uses columns.
    direct_matrix = _matrix_tuple(lattice.matrix.T / a_angstrom)

    sites: list[Site] = []
    for site in structure.sites:
        weight = sum(float(occupancy) * float(specie.Z) for specie, occupancy in site.species.items())
        frac = site.frac_coords
        sites.append((float(frac[0] % 1.0), float(frac[1] % 1.0), float(frac[2] % 1.0), weight))
    if not sites or not any(site[3] > 0.0 for site in sites):
        raise ValueError("CIF contains no occupied atomic sites")

    identity: Matrix3 = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    crystal_system = "triclinic"
    rotations: set[Matrix3] = {identity}
    try:
        analyzer = SpacegroupAnalyzer(structure, symprec=0.01)
        crystal_system = str(analyzer.get_crystal_system()).lower()
        for operation in analyzer.get_symmetry_operations(cartesian=False):
            rotations.add(_matrix_tuple(operation.rotation_matrix, integral=True))
    except (ValueError, TypeError, RuntimeError):
        pass

    return CrystalModel(
        name=_clean_name(name if name is not None else cif_path.stem),
        crystal_system=crystal_system,
        is_hexagonal=crystal_system == "hexagonal",
        a_nm=a_angstrom * 0.1,
        direct_matrix=direct_matrix,
        sites=tuple(sites),
        symmetry_rotations=tuple(sorted(rotations)),
        cif_path=cif_path,
    )


def user_library_dir() -> Path:
    """Return ZAF's normal per-user application-data location."""
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Application Support" / "ZAF"
    if os.name == "nt":
        base = os.environ.get("APPDATA")
        return (Path(base) if base else Path.home() / "AppData" / "Roaming") / "ZAF"
    base = os.environ.get("XDG_DATA_HOME")
    return (Path(base) if base else Path.home() / ".local" / "share") / "ZAF"


class CrystalLibrary:
    """Persist original CIF copies and a small name-to-file index.

    Imported names are made unique rather than replacing an existing entry.
    Existing CIFs remain untouched if indexing a new import fails.
    """

    def __init__(self, root: str | Path | None = None) -> None:
        self.root = Path(root).expanduser() if root is not None else user_library_dir()
        self.cif_dir = self.root / "crystals"
        self.index_path = self.root / "crystals.json"

    def _records(self) -> list[dict[str, str]]:
        if not self.index_path.exists():
            return []
        try:
            data = json.loads(self.index_path.read_text(encoding="utf-8"))
            records = data["crystals"]
            if not isinstance(records, list):
                raise ValueError("'crystals' must be a list")
        except (OSError, KeyError, TypeError, json.JSONDecodeError, ValueError) as exc:
            raise ValueError(f"Could not read ZAF crystal library: {exc}") from exc
        valid: list[dict[str, str]] = []
        for record in records:
            if not isinstance(record, dict):
                continue
            name, filename = record.get("name"), record.get("file")
            if not isinstance(name, str) or not isinstance(filename, str):
                continue
            if not re.fullmatch(r"[0-9a-f]{32}\.cif", filename):
                continue
            valid.append({"name": _clean_name(name), "file": filename})
        return valid

    def _write_records(self, records: list[dict[str, str]]) -> None:
        self.root.mkdir(parents=True, exist_ok=True)
        payload = {"version": 1, "crystals": records}
        temporary: Path | None = None
        try:
            with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=self.root, prefix=".crystals-", suffix=".json", delete=False) as stream:
                temporary = Path(stream.name)
                json.dump(payload, stream, indent=2, ensure_ascii=False)
                stream.write("\n")
            os.replace(temporary, self.index_path)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)

    def list_models(self) -> list[CrystalModel]:
        """List all readable saved structures, sorted by display name."""
        models: list[CrystalModel] = []
        for record in self._records():
            path = self.cif_dir / record["file"]
            if not path.is_file():
                continue
            try:
                models.append(from_cif(path, name=record["name"]))
            except ValueError:
                continue
        return sorted(models, key=lambda model: model.name.casefold())

    def import_cif(self, path: str | Path) -> CrystalModel:
        """Validate and save a CIF, returning a model with its durable path."""
        original = from_cif(path)
        records = self._records()
        existing = {record["name"].casefold() for record in records}
        name = original.name
        suffix = 2
        while name.casefold() in existing:
            suffix_text = f" ({suffix})"
            name = _clean_name(original.name[: 100 - len(suffix_text)].rstrip() + suffix_text)
            suffix += 1

        self.cif_dir.mkdir(parents=True, exist_ok=True)
        saved_path = self.cif_dir / f"{uuid4().hex}.cif"
        try:
            shutil.copyfile(original.cif_path, saved_path)
            # Parse the saved copy before adding it to the index. This also
            # catches any interrupted or incomplete copy.
            saved = from_cif(saved_path, name=name)
            self._write_records([*records, {"name": name, "file": saved_path.name}])
        except Exception:
            saved_path.unlink(missing_ok=True)
            raise
        return saved

    def delete(self, name: str) -> bool:
        """Remove one saved custom crystal by its exact display name."""
        records = self._records()
        removed = [record for record in records if record["name"] == name]
        if not removed:
            return False
        self._write_records([record for record in records if record["name"] != name])
        for record in removed:
            (self.cif_dir / record["file"]).unlink(missing_ok=True)
        return True

    def load(self, name: str) -> CrystalModel:
        """Load one saved structure by display name."""
        for record in self._records():
            if record["name"] == name:
                return from_cif(self.cif_dir / record["file"], name=name)
        raise KeyError(name)
