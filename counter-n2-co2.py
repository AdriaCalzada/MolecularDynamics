import numpy as np
from scipy.spatial import cKDTree
import time

# Open the input file and read the number of atoms
with open("lammpstrj.xyz") as f:
    lines = f.readlines()
n_atoms = int(lines[0])
total_molecules = 128

def atom_index():
    nitrogen_index = 6
    carbon_index = 4
    return nitrogen_index, carbon_index


# Function to read the atom positions from a list of lines in an XYZ file
def read_xyz_file(frame_lines):
    data = np.array([line.split() for line in frame_lines], dtype=float)
    atoms = data[:, 0].astype(int)
    positions = data[:, 1:]
    return atoms, positions


def find_co2_molecules(atoms, carbon_index):
    carbon_mask = atoms == carbon_index

    carbon_atoms = np.where(carbon_mask)

    co2 = []
    co2.append(carbon_atoms)

    return co2, carbon_mask


def find_n2_molecules(atoms, nitrogen_index):
    # Create a boolean mask for nitrogen atoms
    nitrogen_mask = atoms == nitrogen_index
    
    # Get the indices of nitrogen atoms
    nitrogen_atoms = np.where(nitrogen_mask)[0]

    n2 = []
    n2.append(nitrogen_atoms)
    
    return n2, nitrogen_mask


# Function to find the nearest carbon atom to a given CO2 molecule
def find_nearest_carbon(co2_positions, material_carbon_positions, tree):
    # Use cKDTree to find the nearest neighbors
    distances, indices = tree.query(co2_positions, k=1)

    # Get positions of the nearest carbon atoms
    nearest_carbon_positions = material_carbon_positions[indices]

    return nearest_carbon_positions


# Function to find the nearest carbon atom to a given N2 molecule
def find_nearest_material_carbon_to_n2(n2_positions, material_carbon_positions, tree):
    # Use cKDTree to find the nearest neighbors
    distances, indices = tree.query(n2_positions, k=1)

    # Get positions of the nearest carbon atoms
    nearest_carbon_positions = material_carbon_positions[indices]

    return nearest_carbon_positions


def count_crossing_co2_molecules(co2_molecules, positions, tree, material_carbon_positions, carbon_mask):
    # Extract the positions of the nitrogen atoms in the N2 molecules
    carbon_positions = positions[carbon_mask]

    # Find the nearest carbon atom to each N2 molecule using the pre-built cKDTree
    nearest_carbon_positions = find_nearest_carbon(carbon_positions, material_carbon_positions, tree)

    # Count how many N2 molecules are crossing the material by comparing z-coordinates
    n_crossing = np.sum(carbon_positions[:, 2] < nearest_carbon_positions[:, 2])

    return n_crossing


def count_crossing_n2_molecules(n2_molecules, positions, tree, material_nitrogen_positions, nitrogen_mask):
    # Extract the positions of the nitrogen atoms in the N2 molecules
    n2_nitrogen_positions = positions[nitrogen_mask]

    # Find the nearest carbon atom to each N2 molecule using the pre-built cKDTree
    nearest_carbon_positions = find_nearest_material_carbon_to_n2(n2_nitrogen_positions, material_nitrogen_positions, tree)

    # Count how many N2 molecules are crossing the material by comparing z-coordinates
    n_crossing = np.sum(n2_nitrogen_positions[:, 2] < nearest_carbon_positions[:, 2])

    return n_crossing


def main():
    start= time.time()

    nitrogen_index, carbon_index = atom_index()
    
    with open(input("Enter the name of the output:\n"), mode='w') as file_object:

        file_object.write("Frame N\u2082 CO\u2082 \n")

        print("--------------------------\n The program has started")
    
        for frame in range(16001):

            frame_lines = lines[frame*(n_atoms+2)+2 : (frame+1)*(n_atoms+2)]
            atoms, positions = read_xyz_file(frame_lines)
            
            material_carbon_positions = positions[atoms == 1]
            tree = cKDTree(material_carbon_positions)

            n2_molecules, nitrogen_mask = find_n2_molecules(atoms, nitrogen_index)
            n_n2_crossings = count_crossing_n2_molecules(n2_molecules, positions, tree, material_carbon_positions, nitrogen_mask)
            co2_molecules, carbon_mask= find_co2_molecules(atoms, carbon_index)
            n_co2_crossings = count_crossing_co2_molecules(co2_molecules, positions, tree, material_carbon_positions, carbon_mask)

            print('Frame',frame,':\n',n_n2_crossings,'N\u2082', n_co2_crossings,' CO\u2082')

            # Calculate pressure on the right side of the membrane
        

            file_object.write(f"{frame} {n_n2_crossings} {n_co2_crossings}\n")
             
    end=time.time()
    print(' The program has finished\n--------------------------')
    print(f'Time: {end-start}')


main()
