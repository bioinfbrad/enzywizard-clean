
# Command: enzywizard-clean


EnzyWizard-Clean is a command-line tool for cleaning protein structures, generating 
multi-format protein files (CIF, PDB, and FASTA), and providing a detailed traceable
cleaning report. Using PDBFixer APIs, it removes non-protein components (including ligands and water),
replaces non-standard residues with standard ones, repairs missing
heavy atoms, optionally adds hydrogens, and converts the structure 
into a continuous protein chain with standardized residue numbering.
The tool outputs cleaned structure files together with a JSON report summarizing 
residue mapping and cleaning statistics.


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
     - how residue numbering changed,
     - whether residue names were standardized,
     - and how hydrogen content changed after optional hydrogen addition.

   - "clean_statistics"
     A dictionary summarizing the overall cleaning process.

     It includes:
     - "removed_heterogen"
       Number of non-protein residues removed (including ligands and water).

     - "changed_resname"
       Number of residues identified as non-standard and replaced by PDBFixer.

     - "fixed_residues"
       Number of residues where missing atoms were repaired.

     - "added_heavy_atoms"
       Total number of heavy atoms added during missing atom reconstruction.

     - "added_hydrogen_atoms"
       Number of hydrogen atoms added (if hydrogen addition is enabled).

     - "kept_residues"
       Number of residues kept in the final cleaned structure.


# Process:

This command processes the input protein structure as follows:

1. Load the input structure
   - Read the CIF or PDB file using PDBFixer.
   - Resolve the protein name from the input filename.

2. Validate basic input conditions
   - Check that the input file exists.
   - Check that the pH value is within the valid range.

3. Clean the structure (PDBFixer-based processing)
   - Keep only the first chain.
   - Identify non-standard residues.
   - Remove all heterogens (including water and ligands).
   - Replace non-standard residues with standard residues.
   - Disable missing residue reconstruction (do not add missing residues).
   - Detect missing atoms.
   - Add missing heavy atoms.
   - Optionally add hydrogens using OpenMM ForceField with specified pH.
   - Check for invalid coordinates (e.g., NaN values).

4. Renumber the structure
   - Rebuild the topology to ensure:
     - single chain (chain ID A),
     - continuous residue numbering starting from 1.

5. Validate the cleaned structure
   - Ensure single model and single chain.
   - Ensure no hetero residues remain.
   - Ensure residue names are standardized.
   - Ensure residue numbering is continuous.
   - Ensure all required backbone and heavy atoms are present.

6. Save outputs
   - Save the cleaned structure in CIF format.
   - Save the cleaned structure in PDB format.
   - Extract and save the cleaned amino acid sequence in FASTA format.
   - Generate and save the JSON report summarizing residue mapping and statistics.


# dependencies:

- Biopython
- OpenMM
- PDBFixer
- NumPy


# references:

- Biopython:
  https://biopython.org/

- OpenMM:
  https://openmm.org/

- PDBFixer:
  https://github.com/openmm/pdbfixer

- wwPDB Chemical Component Dictionary:
  https://www.wwpdb.org/data/ccd
