from rdkit import Chem
import numpy as np
import os
import glob
import warnings
warnings.filterwarnings('ignore', 'Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.')
import csv
import re
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

def compute_mol2_centroid(mol):
    conformer = mol.GetConformer()
    centroid = rdMolTransforms.ComputeCentroid(conformer)
    return centroid

def calculate_distances(ref_mol_file, mol_file):
    # 读取参考分子
    ref_mol = Chem.MolFromMol2File(ref_mol_file, sanitize=False)

    # 读取待处理分子
    mols = Mol2MolSupplier(mol_file)

    mols.append(ref_mol)

    center_coord = []
    distances = []
    for mol in mols:
        if mol is None:
            # 跳过当前循环，对距离矩阵赋值为-1
            distances.append(-1)
            continue
        
        coord = []
        for i, atom in enumerate(mol.GetAtoms()):
            atomic_number = atom.GetAtomicNum()
            # 如果是氢原子，则跳过当前循环
            if atomic_number == 1:
                continue
            coord.append([mol.GetConformer().GetAtomPosition(i).x,
                          mol.GetConformer().GetAtomPosition(i).y,
                          mol.GetConformer().GetAtomPosition(i).z])
        center_coord.append(np.average(coord, 0))
    
    ref_center = center_coord[-1]
    
    for i in center_coord[:-1]:
        if i is None:
            # 如果中心坐标为空，跳过当前循环，对距离矩阵赋值为-1
            distances.append(-1)
            continue
        distances.append(np.sqrt(np.sum((ref_center - i) ** 2)))
    
    return distances


def compute_distances(ref_mol_file, mol_file):
    ref_mol = Chem.MolFromMol2File(ref_mol_file, sanitize=False)
    mols = Mol2MolSupplier(mol_file)
    mols.append(ref_mol)

    coord_list = []
    distances = []

    for mol in mols:
        corrdmol = compute_mol2_centroid(mol)
        coord = [corrdmol.x, corrdmol.y, corrdmol.z]
        coord_list.append(coord)

    ref_center = np.array(coord_list[-1])
    
    for coord in coord_list[:-1]:
        coord_array = np.array(coord)
        dist = np.sqrt(np.sum((ref_center - coord_array) ** 2))
        distances.append(dist)
    
    return distances

def Mol2MolSupplier(file=None, sanitize=False):
    mols = []
    with open(file, 'r') as f:
        doc = [line for line in f.readlines()]

    start = [index for (index, p) in enumerate(doc) if '@<TRIPOS>MOLECULE' in p]
    start.append(len(doc) + 1)

    interval = list(zip(start[:-1], start[1:]))
    for i in interval:
        block = ",".join(doc[i[0]:i[1]]).replace(',', '')
        m = Chem.MolFromMol2Block(block, sanitize=sanitize)
        mols.append(m)
    return mols
    
def set_dataset(folder_path_ligand, folder_path_dock):
    subfolders_ligand = [f.path for f in os.scandir(folder_path_ligand) if f.is_dir()]
    subfolders_dock = [f.path for f in os.scandir(folder_path_dock) if f.is_dir()]

    ligand_dict = {}
    dock_dict = {}

    for folder in subfolders_ligand:
        for file_path in glob.glob(os.path.join(folder, '*.mol2')):
            filename = os.path.basename(file_path)
            ligand_dict[filename] = file_path

    for folder in subfolders_dock:
        for file_path in glob.glob(os.path.join(folder, '*.mol2')):
            filename = os.path.basename(file_path)
            dock_dict[filename] = file_path
    
    return ligand_dict, dock_dict



def calculate_selected_distances(ligand_dict, dock_dict, method, output_filename):
    results = []
    
    with open(output_filename, 'w', newline='') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(['Ligand', 'Distance1', 'Distance2', 'Distance3', 'Distance4', 'Distance5', 'Distance6', 'Distance7', 'Distance8', 'Distance9',...])  # 写入表头，包括列名
        
        for ligand_filename in ligand_dict:
            if ligand_filename in dock_dict:
                try:
                    ligand_l = ligand_dict[ligand_filename]
                    dock_l = dock_dict[ligand_filename]
                    ligand_name = os.path.basename(ligand_l)  # 获取ligand文件的名称
                    
                    if method == "calculate_distances":
                        dist = calculate_distances(ligand_l, dock_l)
                    elif method == "compute_distances":
                        dist = compute_distances(ligand_l, dock_l)
                    
                    results.append(dist)
                    # 将 [ligand_name] + dist 写入 CSV 文件，以逗号分隔
                    writer.writerow([ligand_name] + dist)
                    
                except AttributeError:
                    continue

'''
folder_path_ligand = r'/mnt/c/Users/86173/rmsd/origin/origin'
folder_path_dock = r'/mnt/c/Users/86173/rmsd/output_origin_pdbbind_site_specific.tar/output_origin_pdbbind_site_specific'
'''

folder_path_ligand = 'C:/Users/lenovo/origin/'
folder_path_dock = 'C:/Users/lenovo/output_origin_pdbbind_site_specific.tar/output_origin_pdbbind_site_specific'

ligand_dict, dock_dict = set_dataset(folder_path_ligand, folder_path_dock)

calculate_selected_distances(ligand_dict, dock_dict, "compute_distances",'results_spec_dict_rdkit0.csv')