from __future__ import annotations
from ..utils.logging_utils import Logger
from pathlib import Path
from ..utils.IO_utils import file_exists, get_stem, check_filename_length, load_protein_structure, write_cif ,write_pdb, write_json_from_dict_inline_leaf_lists, write_fasta, load_pdbfixer, write_cif_from_pdbfixer,write_pdb_from_pdbfixer
from ..algorithms.clean_algorithms import clean_pdbfixer_to_single_chain_A, generate_clean_report, check_cleaned_structure
from ..utils.common_utils import get_optimized_filename

def run_clean_service(input_path: str | Path,output_dir: str | Path,add_H: bool = True,pH: float = 7.0,force_field_file: str = "charmm36.xml") -> bool:
    # ---- logger ----
    input_path = Path(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger = Logger(output_dir)
    logger.print(f"[INFO] Clean processing started: {input_path}")

    # ---- check input ----
    if not (0.0 <= pH <= 14.0):
        logger.print(f"[ERROR] Invalid pH value: {pH}. Must be between 0 and 14.")
        return False

    if not file_exists(input_path):
        logger.print(f"[ERROR] Input not found: {input_path}")
        return False

    # ---- get name ----
    name = get_stem(input_path)
    if not check_filename_length(name, logger):
        return False
    logger.print(f"[INFO] Protein name resolved: {name}")


    # ---- load fixer ----
    fixer = load_pdbfixer(input_path, logger)
    if fixer is None:
        logger.print(f"[ERROR] Failed to load PDBFixer: {input_path}")
        return False

    logger.print("[INFO] PDBFixer loaded")

    # ---- run algorithm ----
    logger.print("[INFO] Cleaning started")
    clean_result = clean_pdbfixer_to_single_chain_A(fixer=fixer,add_H=add_H,logger=logger,pH=pH,force_field_file=force_field_file)
    if clean_result is None:
        return False

    cleaned_fixer, mapping_old_to_new, stats = clean_result

    # ---- save cif/pdb ----
    cleaned_cif_path = output_dir / get_optimized_filename(f"cleaned_{name}.cif")
    cleaned_pdb_path = output_dir / get_optimized_filename(f"cleaned_{name}.pdb")


    write_cif_from_pdbfixer(cleaned_fixer, cleaned_cif_path)
    logger.print(f"[INFO] Cleaned CIF saved: {cleaned_cif_path}")

    write_pdb_from_pdbfixer(cleaned_fixer, cleaned_pdb_path)
    logger.print(f"[INFO] Cleaned PDB saved: {cleaned_pdb_path}")

    # ---- check structure ----
    structure = load_protein_structure(input_path, name, logger)
    if structure is None:
        logger.print(f"[ERROR] Failed to load original structure: {input_path}")
        return False

    cleaned_structure = load_protein_structure(cleaned_cif_path, name, logger)
    if cleaned_structure is None:
        logger.print(f"[ERROR] Failed to load cleaned structure: {cleaned_cif_path}")
        return False

    if not check_cleaned_structure(cleaned_structure, logger):
        return False

    logger.print(f"[INFO] Structure checked")

    # ---- save results ----
    cleaned_fasta_path = output_dir / get_optimized_filename(f"cleaned_{name}.fasta")
    json_report_path = output_dir / get_optimized_filename(f"clean_report_{name}.json")

    if write_fasta(cleaned_structure, cleaned_fasta_path, logger):
        logger.print(f"[INFO] Cleaned FASTA saved: {cleaned_fasta_path}")
    else:
        return False

    report = generate_clean_report(structure=structure,cleaned_structure=cleaned_structure,mapping_old_to_new=mapping_old_to_new,stats=stats,logger=logger)
    if report is None:
        return False

    write_json_from_dict_inline_leaf_lists(report, json_report_path)
    logger.print(f"[INFO] Report JSON saved: {json_report_path}")

    logger.print("[INFO] Clean processing finished")

    return True
