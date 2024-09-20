from Bio import PDB
from Bio.PDB.Atom import Atom as BioAtom
import copy
import warnings
import numpy as np
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import os
from sys import stdout

def add_methyl_to_uracil(residue):
    if 'C7' in residue:
        return

    try:
        C4 = residue['C4']
        C5 = residue['C5']
        C6 = residue['C6']
    except KeyError:
        print(f"Warning: Required atoms not found in residue {residue.id}. Skipping methyl addition.")
        return

    C5_vector = C5.coord - C4.coord
    C6_vector = C6.coord - C5.coord
    methyl_vector = C5_vector + C6_vector
    methyl_vector = methyl_vector / np.linalg.norm(methyl_vector)
    C7_coord = C5.coord + methyl_vector * 1.5  # 假设C-C键长为1.5埃
    C7 = BioAtom("C7", C7_coord, 0.0, 1.0, " ", "C7", None, element="C")
    residue.add(C7)

    # 删除H5原子
    if 'H5' in residue:
        residue.detach_child('H5')

    # 添加三个氢原子，基于四面体结构
    tetrahedral_angle = 109.5 * (np.pi / 180)  # 角度转换为弧度
    bond_length = 1.09  # 假设C-H键长为1.09埃

    # 计算氢原子的位置
    H1_vector = np.array([bond_length, 0.0, 0.0])
    H2_vector = bond_length * np.array([np.cos(tetrahedral_angle), np.sin(tetrahedral_angle), 0.0])
    H3_vector = bond_length * np.array([np.cos(tetrahedral_angle), -np.sin(tetrahedral_angle), 0.0])

    H1_coord = C7_coord + H1_vector
    H2_coord = C7_coord + H2_vector
    H3_coord = C7_coord + H3_vector

    H1 = BioAtom("H71", H1_coord, 0.0, 1.0, " ", "H71", None, element="H")
    H2 = BioAtom("H72", H2_coord, 0.0, 1.0, " ", "H72", None, element="H")
    H3 = BioAtom("H73", H3_coord, 0.0, 1.0, " ", "H73", None, element="H")

    residue.add(H1)
    residue.add(H2)
    residue.add(H3)

def process_rna_to_dna(input_file):
    warnings.simplefilter('ignore', PDBConstructionWarning)
    parser = PDB.PDBParser(QUIET=True)
    temp_cleaned_path = ''
    try:
        structure = parser.get_structure("RNA", input_file)
    except Exception as e:
        print(f"Error parsing the file: {e}")
        print("Trying to parse the file line by line...")
        
        with open(input_file, 'r') as f:
            lines = f.readlines()
        
        valid_lines = [line for line in lines if line.startswith("ATOM") and len(line.split()) >= 11]
        input_name = os.path.splitext(os.path.basename(input_file))[0]
        temp_cleaned_path = f"temp_cleaned_{input_name}.pdb" 
        with open(temp_cleaned_path, 'w') as f:
            f.writelines(valid_lines)
        
        structure = parser.get_structure("RNA", temp_cleaned_path)

    new_structure = copy.deepcopy(structure)

    for model in new_structure:
        for chain in model:
            for residue in chain:
                if residue.resname == "A":
                    residue.resname = "DA"
                elif residue.resname == "G":
                    residue.resname = "DG"
                elif residue.resname == "C":
                    residue.resname = "DC"
                elif residue.resname == "U":
                    residue.resname = "DT"
                    add_methyl_to_uracil(residue)
                
                # Remove O2' and HO2'
                atoms_to_remove = ["O2'", "HO2'"]
                for atom_name in atoms_to_remove:
                    if atom_name in residue:
                        residue.detach_child(atom_name)

                # Add H2''
                if "C2'" in residue:
                    C2_prime = residue["C2'"]
                    H2_double_prime_coord = C2_prime.coord + np.array([0.9572, 0.0, 0.0])
                    H2_double_prime = BioAtom("H2''", H2_double_prime_coord, 0.0, 1.0, " ", "H2''", None, element="H")
                    residue.add(H2_double_prime)

                # Renumber atoms
                atoms = list(residue.get_atoms())
                for i, atom in enumerate(atoms, start=1):
                    atom.serial_number = i

    # # Use input_file instead of input_pdb_path
    # output_file = os.path.splitext(input_file)[0] + "_ss.pdb"
    # io = PDB.PDBIO()
    # io.set_structure(new_structure)
    # io.save(output_file)
    # Generate the output file name
    input_dir = os.path.dirname(input_file)
    input_filename = os.path.basename(input_file)
    if input_filename.startswith("d_"):
        output_filename = "ss" + input_filename[2:]  # Remove "d_" prefix
    else:
        output_filename = "ss" + input_filename

    output_file = os.path.join(input_dir, output_filename)
    output_file = os.path.splitext(output_file)[0] + ".pdb"  # Ensure .pdb extension

    # Normalize the path
    output_file = os.path.normpath(output_file)

    io = PDB.PDBIO()
    io.set_structure(new_structure)
    io.save(output_file)

    if os.path.exists(temp_cleaned_path):
        os.remove(temp_cleaned_path)

    return "./" + output_file.replace("\\", "/").lstrip("./")

# Example use (commented out)
# if __name__ == "__main__":
#     input_pdb_path = "aptamer_pdb/dna_ROGERRFODIWHRFIFSHIHOWOYPATISIPESIQ.pdb"
#     output_pdb_path = process_rna_to_dna(input_pdb_path)
#     print(f"Converted ssDNA file saved as: {output_pdb_path}")