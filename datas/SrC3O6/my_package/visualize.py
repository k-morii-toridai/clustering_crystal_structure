import subprocess


def vesta(POSCAR='./POSCAR'):
    """
    for visualizing POSCAR file (or cif file) by using VESTA
    
    pram1:POSAR:POSCAR file path
    
    """
    VESTA = '/home/morii-k/vesta/VESTA-gtk3/VESTA'
    POSCAR = POSCAR
    myprocess = subprocess.Popen([VESTA, '-open', POSCAR])
    