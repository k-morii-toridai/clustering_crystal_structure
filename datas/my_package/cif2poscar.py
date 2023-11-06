import os
from pymatgen.io.cif import CifParser


def cif2poscar_mk_dir(cif_file):
    
    
    def mk_cif_num_folder(cif_file):
        """
        To create a new dirctory which name is CIF file number, Use thie func().

        param1: example: cif_file='1507756.cif'
        created: a directory which name is a CIF file number
        """
        cif_file_number = cif_file.split(".")[0]
        os.mkdir(cif_file_number)
    
    
    def cif2poscar(cif_file):
        """
        To convert a cif to POSCAR file, Use this func().

        param1: example: cif_file='1507756.cif'
        created: a POSCAR file
        """
        parser = CifParser(cif_file)
        structure = parser.get_structures()[0]
        # make cif file number
        cif_file_number = cif_file.split(".")[0]
        # Createしたフォルダに，POSCARファイルとして書き出し
        structure.to(fmt="poscar", filename=f"{cif_file_number}/POSCAR") 
    
    
    # Create new folder
    mk_cif_num_folder(cif_file)
    # convert cif to POSCAR
    cif2poscar(cif_file)
