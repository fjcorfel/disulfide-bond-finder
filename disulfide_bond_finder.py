import sys
import math
from Bio.PDB import PDBParser, Structure, Residue, vectors


def get_file_name() -> str:
    """
    Reads the PDB file name from de command line or input 
    and check if it is valid.

    Returns:
        str: PDB file name.
    """
    
    if len(sys.argv) == 1:
        file_name = input("Enter file name (PDB): ").strip()
    elif len(sys.argv) == 2:
        file_name = sys.argv[1]
    else:
        sys.exit("Too many arguments used.\nUsage: python pdb_diS.py <example_protein.pdb>")
        
    if not file_name.endswith(".pdb"):
        sys.exit("Invalid file extension.")
    
    return file_name


def read_pdb_file(file_name: str) -> Structure:
    """
    Returns the structure contained in a PDB file.

    Args:
        file_name (str): PDB file name.

    Returns:
        Structure: Protein structure.
    """
    
    parser = PDBParser(QUIET=True)
    
    try:
        structure = parser.get_structure("protein", file_name)
    except FileNotFoundError:
        sys.exit(f"File '{file_name}' not found.")
    except Exception as e:
        sys.exit(f"Error parsing PDB file: {e}")
    
    return structure


def discard_residue(res: Residue, model_type: str) -> bool:
    """
    Determines if a residue should be discarded based on B-factor and pLDDT,
    depending on model type (experimental or AlphaFold model).

    Args:
        res (Residue): Residue to check.
        model_type (str): Type of model ('alphafold' or 'experimental')

    Returns:
        bool: True if residue should be discarded, False otherwise.
    """
    
    threshold = 50 if model_type == "alphafold" else 30
    
    for atom in res:

        if ((model_type == "alphafold" and atom.get_bfactor() < threshold) or
            (model_type == "experimental" and atom.get_bfactor() > threshold)):
            return True
    
    return False
    

def get_cysteines(structure: Structure) -> list:
    """
    Extracts CYS residues from the main structure.

    Args:
        structure (Structure): Protein structure read from a PDB file.

    Returns:
        list: List of CYS residues in the protein structure.
    """
    
    structure_name = structure.header.get("name")
    
    if "alphafold" in structure_name.strip().lower():
        model_type = "alphafold"
    else:
        model_type = "experimental"

    cysteines = []
    
    for model in structure:
        for chain in model:
            for res in chain:
                if (res.get_resname() == "CYS" and 
                    not discard_residue(res, model_type)):
                    cysteines.append(res)
                    
    return cysteines


def calculate_cysteines_distance(cys1: Residue, cys2: Residue) -> float:
    """
    Calculates the Euclidean distance  (angstrom) between two CYS sulfur atoms.

    Args:
        res1 (Residue): First CYS residue.
        res2 (Residue): Second CYS residue.

    Returns:
        float: Distance between CYS sulfur atoms.
    """
    
    return cys1["SG"] - cys2["SG"]


def calculate_dihedral_angle(cys1: Residue, cys2: Residue) -> float:
    """
    Calculates the dihedral angle (degrees) between two CYS residues.

    Args:
        cys1 (Residue): First CYS residue.
        cys2 (Residue): Second CYS residue.

    Returns:
        float: Dihedral angle between CYS residues.
    """

    atom1 = cys1["CB"]
    atom2 = cys1["SG"]
    atom3 = cys2["SG"]
    atom4 = cys2["CB"]

    angle_radians = vectors.calc_dihedral(
        atom1.get_vector(), atom2.get_vector(),
        atom3.get_vector(), atom4.get_vector()
        )
    
    angle_degrees = math.degrees(angle_radians)

    return abs(angle_degrees)


def find_disulfide_bonds(cysteines: list) -> list[tuple]:
    """
    Finds potential disulfide bonds between cysteine residues.

    Args:
        cysteines (list): List of CYS residues.

    Returns:
        list[tuple]: List of tuples containing CYS pairs [0-1], distance [2] and
        dihedral angle [3].
    """
    
    disulfide_bonds = []
    
    for i in range(len(cysteines)):
        cys1 = cysteines[i]
        
        for j in range(i + 1, len(cysteines)):
            cys2 = cysteines[j]
        
            distance = calculate_cysteines_distance(cys1, cys2)
            if 1.5 <= distance <= 2.5:
                    
                angle = calculate_dihedral_angle(cys1, cys2)
                if 84 <= angle <= 96:
                    
                    disulfide_bonds.append((cys1, cys2, distance, angle))
    
    return disulfide_bonds
    
    
def main():
    
    print("-------------------------")
    print("PDB Disulfide Bond Finder")
    print("-------------------------\n")
    
    file_name = get_file_name()
    structure = read_pdb_file(file_name)
    cysteines = get_cysteines(structure)
    disulfide_bonds = find_disulfide_bonds(cysteines)

    if disulfide_bonds:    
        print(f"Potential disulfide bonds found in '{file_name}': ")
        for bond in disulfide_bonds:
            print(f"* CYS {bond[0].id[1]} - CYS {bond[1].id[1]} | Distance: {bond[2]:.2f} Å | Angle: {bond[3]:.2f}º")
    else:
        print(f"Potential disulfide bonds NOT found in '{file_name}'.")


if __name__ == "__main__":
    main()
