[![DOI](https://zenodo.org/badge/1219033403.svg)](https://doi.org/10.5281/zenodo.19709411)

# EnzyWizard-Clean


EnzyWizard-Clean is a command-line tool for cleaning protein structures, generating 
multi-format protein files (CIF, PDB, and FASTA), and providing a detailed traceable
cleaning report. It standardizes residue names, removes problematic residues 
(non-standard residues, residues with missing backbone atoms, residues with 
missing required heavy atoms, residues with unexpected heavy atoms, residues 
with invalid occupancy), repairs residue order by renumbering residues continuously,
converts the structure into a cleaned protein chain, optionally adds hydrogens 
using OpenMM, and outputs cleaned structure files together with a JSON report 
summarizing residue mapping and cleaning statistics.



# example usage:

Example command:

enzywizard-clean -i examples/input/3GP6.cif -o examples/output/



# input parameters:

-i, --input_path
Required.
Path to the input protein structure file in CIF or PDB format.

-o, --output_dir
Required.
Path to the output directory for saving cleaned structure files and the JSON report.

--no_add_H
Optional flag.
Disable hydrogen addition using OpenMM.
By default, hydrogens are added.

--pH
Optional.
pH value used for hydrogen addition.
Default: 7.0
Valid range: 0.0 to 14.0


# output content:

The program outputs the following files into the output directory:

1. Cleaned structure files
   - cleaned_{name}.cif
   - cleaned_{name}.pdb
   - cleaned_{name}.fasta

2. A JSON report
   - clean_report_{name}.json

   The JSON report contains:

   - "output_type"
     A string identifying the report type:
     "enzywizard_clean"

   - "amino_acid_mapping_old_to_new"
     A list describing how residues in the original structure correspond to
     residues in the cleaned structure.

     Each entry contains:
     - "old_residue"
       Information for the original residue before cleaning:
       - "aa_id": original residue index
       - "aa_name": original residue one-letter amino acid code
       - "hydrogen_atom_count": number of hydrogen atoms found in the original residue

     - "new_residue"
       Information for the cleaned residue after cleaning:
       - "aa_id": new residue index after continuous renumbering
       - "aa_name": cleaned residue one-letter amino acid code
       - "hydrogen_atom_count": number of hydrogen atoms in the cleaned residue

     This mapping helps users track:
     - which residues were kept,
     - how residue numbering changed,
     - whether residue names were standardized,
     - and how hydrogen content changed after optional hydrogen addition.

   - "clean_statistics"
     A dictionary summarizing the overall cleaning process.

     It includes:
     - "changed_resname"
       Number of residues whose names were standardized using the MODRES mapping.

     - "removed_nonstd"
       Number of non-standard residues removed.

     - "removed_missing_bb"
       Number of residues removed because required backbone atoms
       (N, CA, C) were missing.

     - "removed_missing_heavy_atoms"
       Number of residues removed because one or more required heavy atoms
       were missing.

     - "removed_unexpected_heavy_atoms"
       Number of residues removed because unexpected heavy atoms were present.

     - "removed_bad_occ"
       Number of residues removed because selected atoms had invalid occupancy.

     - "removed_inscodes"
       Number of residues whose original insertion codes were present and then removed
       during residue renumbering.

     - "kept_residues"
       Number of residues kept in the final cleaned structure.


# Process:

This command processes the input protein structure as follows:

1. Load the input structure
   - Read the CIF or PDB file using Biopython (Bio.PDB).
   - Resolve the protein name from the input filename.

2. Validate basic input conditions
   - Check that the input file exists.
   - Check that the pH value is within the valid range.

3. Clean the structure (Biopython-based processing)
   - Extract a single chain using Biopython structure utilities.
   - Standardize residue names using the MODRES mapping.
   - Remove residues with missing backbone atoms (N, CA, C).
   - Remove residues with missing required heavy atoms.
   - Remove residues with unexpected heavy atoms.
   - Remove residues with invalid occupancy values.
   - Remove insertion codes.
   - Repair discontinuous residue numbering by rebuilding residue indices.
   - Renumber all kept residues continuously starting from 1.
   - Rebuild the output structure as a single chain with chain ID A.

4. Optionally add hydrogens
   - Convert the cleaned Biopython structure into an OpenMM object.
   - Use OpenMM Modeller.addHydrogens() with a specified pH.
   - Convert the hydrogen-added structure back into a Biopython structure.

5. Validate the cleaned structure

6. Save outputs
   - Save the cleaned structure in CIF format (Biopython MMCIFIO).
   - Save the cleaned structure in PDB format (Biopython PDBIO).
   - Extract and save the cleaned amino acid sequence in FASTA format.
   - Generate and save the JSON report summarizing residue mapping and statistics.


# dependencies:

- Biopython
- OpenMM
- NumPy


# references:

- Biopython:
  https://biopython.org/

- OpenMM:
  https://openmm.org/

- wwPDB Chemical Component Dictionary / MODRES-related residue standardization resource:
  https://www.wwpdb.org/data/ccd

- Rosetta structure preparation overview:
  https://docs.rosettacommons.org/docs/latest/rosetta_basics/preparation/preparing-structures
