#!/usr/bin/python3
import sys,subprocess

base = '/mnt/ssd_elecom_black/cif/' # cif base directory, change here
fname = sys.argv[1]
dir = base + fname[0] + '/' + fname[1:3] + '/' + fname[3:5] + '/' + fname
# cmd = "open -a VESTA " + dir
# subprocess.call(cmd.split(" "))
VESTA = '/home/morii-k/vesta/VESTA-gtk3/VESTA'
subprocess.Popen([VESTA, '-open', dir])


