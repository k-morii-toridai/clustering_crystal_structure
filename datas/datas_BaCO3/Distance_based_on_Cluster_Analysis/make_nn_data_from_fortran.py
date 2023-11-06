import os
import re
import pickle
import shutil
from collections import defaultdict
import numpy as np
import subprocess
import glob
from pymatgen.core.structure import IStructure
from pymatgen.io.vasp.inputs import Poscar

def read_nnlist(filepath:str) -> dict[list]:
    """
    Read nnlist file.
    """
    nnfile = open(filepath).readlines()
    nnlist = defaultdict(list)
    for i in nnfile:
        info = re.split('\s+',i)[1:-1]
        center_isite = int(info[0])
        neighbor_isite = int(info[1])
        distance = float(info[2])
        coord = [float(info[3]),float(info[4]),float(info[5])]
        neighbor_isite_info = [int(info[6]),int(info[7]),int(info[8])]
        #add nnlist
        nnlist[center_isite].append((neighbor_isite,distance,coord,neighbor_isite_info))
    return dict(nnlist)

def sort_and_trim(nnlist: dict[list], siteinfo: dict, max_neib: dict) -> dict:
    """
    Sort and trim neighbor lists.
    """
    sorted_nnlist = {}
    for isite, atom in siteinfo.items():
        nnlist[isite].sort(key=lambda x: x[1])
        sorted_nnlist[isite] = nnlist[isite][1:max_neib[atom] + 1]
    return sorted_nnlist

def nnlist2nn_data(nnlist:dict,siteinfo:dict) -> dict():
    #set nn_data
    nn_data=defaultdict(list)
    #center_info=[siteinfo[1],np.nan,0,0,0]
    for isite,neibs in nnlist.items():
        center_info = [siteinfo[isite],None,0,0,0]
        nn_data[isite].append(center_info)
        for neib_isite,_,coord,_ in neibs:
            nn_data_=[neib_isite,siteinfo[neib_isite],*coord]
            nn_data[isite].append(nn_data_)
    return dict(nn_data)

def read_siteinfo_from_poscar(filename:str) -> dict:
    file = open(filename,'r').readlines()#正規表現を使ってdirectの行を探す
    re_direct = re.compile(r'Direct',re.IGNORECASE)
    for i, line in enumerate(file):
        if re_direct.search(line):
            Direct_index = i
            break
    siteinfo = file[Direct_index-2:Direct_index]
    siteinfo = [re.split('\s+',i.strip()) for i in siteinfo]
    siteinfo[1] = [int(i) for i in siteinfo[1]]

    total_atom_num = sum(siteinfo[1])
    #各定義域をリストに格納
    
    atom_num_list = [0]
    for i in siteinfo[1]:
        atom_num_list.append(atom_num_list[-1]+i)
    atom_num_list.pop(0)
    def get_siteinfo(isite):
        for atom_num in atom_num_list:
            if isite <= atom_num:
                return siteinfo[0][atom_num_list.index(atom_num)]
        return None
    return {i:get_siteinfo(i) for i in range(1,total_atom_num+1)}

def run_make_nnlist_out(fileaddress='POSCAR', rmax=1.0):
    """
    Generate nnlist using an external program.
    """
    subprocess.run(f'make_nnlist.out {fileaddress} {rmax} >> {fileaddress}.log', shell=True)

def make_poscar(ciffile,algorithm='pymatgen'):
    """
    Generate a POSCAR file from a CIF file.
    """
    if algorithm == 'pymatgen':
        poscar = Poscar(IStructure.from_file(ciffile))
        poscar.write_file('POSCAR')
    elif algorithm == 'cif2cell':
        subprocess.run(f'cif2cell {ciffile} -p vasp --vasp-format=5', shell=True)

class CIFDataProcessor:
    def __init__(self,max_neib: dict = {'Si': 4, 'O': 2}, algorithm: str = 'cif2cell'):
        self.max_neib = max_neib
        self.algorithm = algorithm

    def make_nn_data_from_cifdirs(self,dirpath,outputpath='result'):
        """
        Generate nn_data.pickle from CIF directories.
        """
        self.dirpath = dirpath
        self.outputpath = outputpath
        self.cwd = os.getcwd()
        result_dir_path = os.path.join(self.outputpath, self.dirpath)

        # Make output directory
        if not os.path.isdir(self.outputpath):
            os.makedirs(result_dir_path)
        else:
            if not os.path.isdir(result_dir_path):
                os.makedirs(result_dir_path)

        ciffilelist = glob.glob(os.path.join(self.dirpath, '*.cif'))

        for cif_path in ciffilelist:
            self.process_cif_file(cif_path, result_dir_path)

    def process_cif_file(self, cif_path, result_dir_path):
        basename = os.path.basename(cif_path)
        cifnum = basename.replace('.cif', '')
        result_path = 'nb_{}.pickle'.format(cifnum)

        cifdataout = os.path.join(result_dir_path, cifnum)
        if not os.path.isdir(cifdataout):
            os.mkdir(cifdataout)
        shutil.copyfile(cif_path, os.path.join(cifdataout, basename))

        os.chdir(cifdataout)
        try:
            make_poscar(basename, algorithm=self.algorithm)  # Use the instance variable for algorithm
            run_make_nnlist_out('POSCAR', 5)
            siteinfo = read_siteinfo_from_poscar('POSCAR')
            nnlist = read_nnlist('POSCAR.nnlist', siteinfo, self.max_neib)
            nn_data = nnlist2nn_data(nnlist, siteinfo)

            # Save nn_data
            self.save_nn_data(result_path, nn_data)
        except Exception as e:
            print(f"Error processing {cif_path}: {str(e)}")
        os.chdir(self.cwd)

    def save_nn_data(self, result_path, nn_data):
        if os.path.isfile(result_path):
            os.remove(result_path)
        with open(result_path, 'wb') as f:
            pickle.dump(nn_data, f)
