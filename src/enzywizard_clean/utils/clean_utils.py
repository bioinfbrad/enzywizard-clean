from __future__ import annotations
from Bio.PDB.Atom import Atom
from typing import Dict, List, Tuple, Optional, Union, Any
from Bio.PDB.Residue import Residue
from ..resources.aa_resources import AA3_STANDARD, modres
from pdbfixer import PDBFixer
from ..utils.logging_utils import Logger
from ..utils.sequence_utils import normalize_aa_name_to_one_letter
from copy import deepcopy
from openmm.app import Topology
from openmm.vec3 import Vec3
from openmm import unit

def standardize_resname(resname: str) -> str:
    resname = resname.strip().upper()
    if resname in AA3_STANDARD:
        return resname
    if resname in modres:
        return modres[resname].strip().upper()
    return resname

def normalize_atom_name(atom_name: str) -> str:
    return atom_name.strip().upper()

def is_hydrogen_atom(atom: Atom) -> bool:
    element = atom.element.strip().upper() if atom.element is not None else ""
    name = normalize_atom_name(atom.get_name())

    if element == "H":
        return True
    if name.startswith("H"):
        return True
    return False

def choose_atom_altloc(atom_list: List[Atom]) -> Atom:
    # Prefer blank altloc
    for a in atom_list:
        if (a.get_altloc() or " ").strip() == "":
            return a

    # Else pick highest occupancy
    best = atom_list[0]
    best_occ = best.get_occupancy()
    best_occ = best_occ if best_occ is not None else -1.0
    for a in atom_list[1:]:
        occ = a.get_occupancy()
        occ = occ if occ is not None else -1.0
        if occ > best_occ:
            best, best_occ = a, occ
    return best

def get_residue_heavy_atom_name_set(res: Residue) -> set[str]:
    heavy_atom_name_set: set[str] = set()

    for atom in res.get_atoms():
        if is_hydrogen_atom(atom):
            continue
        heavy_atom_name_set.add(normalize_atom_name(atom.get_name()))

    return heavy_atom_name_set

def clone_atom(atom: Atom, *, new_coord=None) -> Atom:
    coord = new_coord if new_coord is not None else atom.get_coord()
    normalized_name = normalize_atom_name(atom.get_name())

    return Atom(
        name=normalized_name,
        coord=coord,
        bfactor=atom.get_bfactor(),
        occupancy=atom.get_occupancy(),
        altloc=" ",
        fullname=f"{normalized_name:>4}",
        serial_number=atom.get_serial_number(),
        element=atom.element.strip().upper() if atom.element is not None else atom.element,
    )


def get_single_chain_from_pdbfixer(fixer: PDBFixer,logger: Logger):
    topology = fixer.topology

    if topology is None:
        logger.print(f"[ERROR] No topology found in PDBFixer")
        return None

    for chain in topology.chains():
        return chain

    logger.print(f"[ERROR] No chain found in PDBFixer topology")
    return None

def is_water_residue_name(resname: str) -> bool:
    return resname.strip().upper() in {"HOH", "WAT", "H2O"}

def is_protein_residue_name(resname: str) -> bool:
    return resname.strip().upper() in AA3_STANDARD

def count_single_chain_heterogens_from_pdbfixer_chain(chain) -> int:
    count = 0
    for res in chain.residues():
        if not is_protein_residue_name(res.name):
            count += 1
    return count

def count_hydrogen_atoms_in_fixer(fixer: PDBFixer) -> int:
    count = 0
    for atom in fixer.topology.atoms():
        element = atom.element
        if element is not None and element.symbol == "H":
            count += 1
    return count

def get_aa1_from_resname(resname: str) -> str:
    return normalize_aa_name_to_one_letter(resname)

def get_single_chain_protein_residue_info_from_pdbfixer_chain(chain) -> List[Dict[str, Any]]:
    residue_info_list: List[Dict[str, Any]] = []

    for res in chain.residues():
        resname = res.name.strip().upper()
        if not is_protein_residue_name(resname):
            continue

        try:
            aa_index = int(str(res.id).strip())
        except Exception:
            aa_index = -1

        residue_info_list.append(
            {
                "aa_index": aa_index,
                "aa_name": get_aa1_from_resname(resname),
            }
        )

    return residue_info_list

def build_missing_atom_stats_from_fixer(fixer: PDBFixer) -> Tuple[int, int]:
    fixed_residue_count = 0
    added_heavy_atom_count = 0

    missing_atoms_dict = getattr(fixer, "missingAtoms", None)
    missing_terminals_dict = getattr(fixer, "missingTerminals", None)

    if missing_atoms_dict is None:
        missing_atoms_dict = {}
    if missing_terminals_dict is None:
        missing_terminals_dict = {}

    all_residues = set(missing_atoms_dict.keys()) | set(missing_terminals_dict.keys())

    for res in all_residues:
        atom_num = 0

        atom_list = missing_atoms_dict.get(res, [])
        atom_num += len(atom_list)

        terminal_list = missing_terminals_dict.get(res, [])
        atom_num += len(terminal_list)

        if atom_num > 0:
            fixed_residue_count += 1
            added_heavy_atom_count += atom_num

    return fixed_residue_count, added_heavy_atom_count

def build_nonstandard_residue_count_from_fixer(fixer: PDBFixer) -> int:
    nonstandard_residues = getattr(fixer, "nonstandardResidues", None)
    if nonstandard_residues is None:
        return 0
    return len(nonstandard_residues)

def renumber_single_chain_fixer_residues(
    fixer: PDBFixer,
    logger: Logger
) -> PDBFixer | None:
    old_topology = fixer.topology
    old_positions = fixer.positions

    old_chain_list = list(old_topology.chains())
    if len(old_chain_list) == 0:
        logger.print(f"[ERROR] No chain found in fixer topology")
        return None

    if len(old_chain_list) != 1:
        logger.print(f"[ERROR] Fixer topology is not single-chain before renumbering")
        return None

    old_chain = old_chain_list[0]

    try:
        new_topology = Topology()
        new_chain = new_topology.addChain(id="A")

        atom_mapping: Dict[Any, Any] = {}
        new_positions_raw: List[Vec3] = []

        new_resseq = 0

        for old_res in old_chain.residues():
            new_resseq += 1

            new_res = new_topology.addResidue(
                name=old_res.name,
                chain=new_chain,
                id=str(new_resseq),
                insertionCode=" "
            )

            for old_atom in old_res.atoms():
                new_atom = new_topology.addAtom(
                    name=old_atom.name,
                    element=old_atom.element,
                    residue=new_res,
                    id=str(old_atom.id) if old_atom.id is not None else None,
                    formalCharge=old_atom.formalCharge
                )
                atom_mapping[old_atom] = new_atom
                new_positions_raw.append(old_positions[old_atom.index].value_in_unit(unit.nanometer))

        for bond in old_topology.bonds():
            atom1 = bond[0]
            atom2 = bond[1]

            new_atom1 = atom_mapping.get(atom1)
            new_atom2 = atom_mapping.get(atom2)

            if new_atom1 is None or new_atom2 is None:
                logger.print(f"[ERROR] Failed to map bonded atoms during renumbering")
                return None

            new_topology.addBond(
                new_atom1,
                new_atom2,
                type=bond.type,
                order=bond.order
            )

        new_fixer = deepcopy(fixer)
        new_fixer.topology = new_topology
        new_fixer.positions = unit.Quantity(new_positions_raw, unit.nanometer)

        return new_fixer

    except Exception as e:
        logger.print(f"[ERROR] Exception in renumbering fixer residues: {e}")
        return None

def postprocess_clean_report_to_schema(raw_report: Dict[str, Any],logger: Logger) -> Dict[str, Any] | None:
    """
    Convert the raw EnzyWizard-Clean report to the current FAIR-style JSON report schema.
    """

    if not isinstance(raw_report, dict):
        logger.print("[ERROR] Invalid raw clean report: report is not a dictionary")
        return None

    raw_report_type = raw_report.get("output_type")
    if raw_report_type != "enzywizard_clean":
        logger.print(
            f"[ERROR] Invalid raw clean report type: expected 'enzywizard_clean', got '{raw_report_type}'"
        )
        return None

    raw_mapping = raw_report.get("amino_acid_mapping_old_to_new")
    if not isinstance(raw_mapping, list):
        logger.print("[ERROR] Invalid raw clean report: 'amino_acid_mapping_old_to_new' is missing or not a list")
        return None

    raw_statistics = raw_report.get("clean_statistics")
    if not isinstance(raw_statistics, dict):
        logger.print("[ERROR] Invalid raw clean report: 'clean_statistics' is missing or not a dictionary")
        return None

    residue_mapping_old_to_new: List[Dict[str, Dict[str, Any]]] = []

    for i, mapping_item in enumerate(raw_mapping):
        if not isinstance(mapping_item, dict):
            logger.print(f"[ERROR] Invalid residue mapping item at index {i}: item is not a dictionary")
            return None

        old_residue = mapping_item.get("old_residue")
        new_residue = mapping_item.get("new_residue")

        if not isinstance(old_residue, dict):
            logger.print(f"[ERROR] Invalid old_residue at mapping index {i}: old_residue is not a dictionary")
            return None

        if not isinstance(new_residue, dict):
            logger.print(f"[ERROR] Invalid new_residue at mapping index {i}: new_residue is not a dictionary")
            return None

        old_residue_index = old_residue.get("aa_id")
        old_residue_name = old_residue.get("aa_name")
        old_hydrogen_atom_count = old_residue.get("hydrogen_atom_count")

        new_residue_index = new_residue.get("aa_id")
        new_residue_name = new_residue.get("aa_name")
        new_hydrogen_atom_count = new_residue.get("hydrogen_atom_count")

        if not isinstance(old_residue_index, int):
            logger.print(f"[ERROR] Invalid old residue index at mapping index {i}: {old_residue_index}")
            return None

        if not isinstance(new_residue_index, int):
            logger.print(f"[ERROR] Invalid new residue index at mapping index {i}: {new_residue_index}")
            return None

        if not isinstance(old_residue_name, str) or len(old_residue_name) != 1:
            logger.print(f"[ERROR] Invalid old residue name at mapping index {i}: {old_residue_name}")
            return None

        if not isinstance(new_residue_name, str) or len(new_residue_name) != 1:
            logger.print(f"[ERROR] Invalid new residue name at mapping index {i}: {new_residue_name}")
            return None

        if not isinstance(old_hydrogen_atom_count, int):
            logger.print(
                f"[ERROR] Invalid old hydrogen atom count at mapping index {i}: {old_hydrogen_atom_count}"
            )
            return None

        if not isinstance(new_hydrogen_atom_count, int):
            logger.print(
                f"[ERROR] Invalid new hydrogen atom count at mapping index {i}: {new_hydrogen_atom_count}"
            )
            return None

        residue_mapping_old_to_new.append(
            {
                "old_residue": {
                    "residue_index": old_residue_index,
                    "residue_name": old_residue_name,
                    "hydrogen_atom_count": old_hydrogen_atom_count,
                },
                "new_residue": {
                    "residue_index": new_residue_index,
                    "residue_name": new_residue_name,
                    "hydrogen_atom_count": new_hydrogen_atom_count,
                },
            }
        )

    required_raw_statistics_keys = [
        "removed_heterogen",
        "changed_resname",
        "fixed_residues",
        "added_heavy_atoms",
        "added_hydrogen_atoms",
        "kept_residues",
    ]

    for key in required_raw_statistics_keys:
        if key not in raw_statistics:
            logger.print(f"[ERROR] Missing raw clean statistic field: {key}")
            return None

        if not isinstance(raw_statistics[key], int):
            logger.print(f"[ERROR] Invalid raw clean statistic field '{key}': {raw_statistics[key]}")
            return None

    clean_statistics = {
        "removed_heterogen_count": raw_statistics["removed_heterogen"],
        "standardized_residue_name_count": raw_statistics["changed_resname"],
        "repaired_residue_count": raw_statistics["fixed_residues"],
        "added_heavy_atom_count": raw_statistics["added_heavy_atoms"],
        "added_hydrogen_atom_count": raw_statistics["added_hydrogen_atoms"],
        "retained_residue_count": raw_statistics["kept_residues"],
    }

    return {
        "report_type": "enzywizard_clean",
        "residue_mapping_old_to_new": residue_mapping_old_to_new,
        "clean_statistics": clean_statistics,
    }