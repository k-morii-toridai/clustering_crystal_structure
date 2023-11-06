import os
import re
import sys
import pickle
import argparse
import subprocess
import numpy as np
from tqdm import tqdm

def change_dir(dirpath):
	os.chdir(dirpath)

def create_dataset(datafile):
	with open(datafile) as file1:
		file1data = file1.readlines()
	
	nearest_neighbor = {}
	
	switch1 = False
	switch2 = False
	index_num = 'nodata'
	for index,i in enumerate(file1data):
		if '**isite= ' in i:
			isite_num = int(i.split()[1])
			index_num = index + 1
			switch1 = True
		
			
		if switch1:
			if index_num == index:
				isite_atom = i.split()[0]
				isite_xyz = re.findall('\[(.*)\]',i)[0].split()
				isite_data_list = [isite_atom,float(isite_xyz[0]),float(isite_xyz[1]),float(isite_xyz[2])] 
		
			if 'ce_fraction= ' in i:
				ce_fraction = i.split()[1]
				isite_data_list.insert(1,float(ce_fraction))
			
			if 'species:coords:indexend\n' == i:
				switch1 = False
				switch2 = False
			
			if switch2:
				firstNN_data = i.split(' : ')
				firstNN_site = int(firstNN_data[2].split('\n')[0])
				firstNN_xyz = re.findall('\[(.*)\]',firstNN_data[1])[0].split()
				nearest_neighbor[isite_num].append([firstNN_site,firstNN_data[0],float(firstNN_xyz[0]),float(firstNN_xyz[1]),float(firstNN_xyz[2]),])
			
			if 'species:coords:index\n' == i:
				nearest_neighbor[isite_num] = [isite_data_list]
				switch2 = True

	#print(nearest_neighbor[0])	
	return nearest_neighbor


def main():
	cwd = os.getcwd()
	parser = argparse.ArgumentParser()
	parser.add_argument('--output1', default='result')
	parser.add_argument('--output2', default='cod')
	parser.add_argument('-e','--explanation', default=False)
	args = parser.parse_args()

	if args.explanation:
		print('''cifdir : ['result/cod', 'result/cod/1000007', 'result/cod/1000017',...''')
		print('''cifdir_i_file : 'result/cod/1000007/1000007.txt' ''')
		print('''isite_data_list : ['Ca1', -1.07652858, 6.22898225, 3.78771229]''')
		print('''''')
		sys.exit()
	cifdir_ = subprocess.getoutput("find {0} -type d | sort".format(args.output1 + "/" + args.output2))
	cifdir = cifdir_.split('\n')
	del cifdir[0]
	print('make nn_data.pickle')
	for i in tqdm(cifdir):
		cifdir_i_file = subprocess.getoutput("find {0}/*.txt".format(i))
		try:
			cifdir_o_file = create_dataset(cifdir_i_file)
		except FileNotFoundError:
			subprocess.getoutput("rm -rf {0}".format(i))
			continue
		
		#print(cifdir_o_file)
		
		#cifnum = re.findall('\/(\w+)',i)[0]
		cifnum=os.path.basename(i)
		#print(cifnum)
		with open(("{0}" + "/" + "nb_{1}.dat").format(i,cifnum),"w") as file2:
			for j in cifdir_o_file.keys():
				print('isite = {0}  :  ce_fraction = {1}'.format(j,cifdir_o_file[j][0][1]),file=file2)
				for num,k in enumerate(cifdir_o_file[j]):
					if num == 0:
						print([cifdir_o_file[j][0][0],cifdir_o_file[j][0][2],cifdir_o_file[j][0][3],cifdir_o_file[j][0][4]],file=file2)
					else:
						print(k,file=file2)
				print(file=file2)
		
		with open(("{0}" + "/" + "nb_{1}.pickle").format(i,cifnum),"wb") as fwb:
			pickle.dump(cifdir_o_file,fwb)
	print('end nn_data.pickle')
if __name__ == '__main__':
	main()
