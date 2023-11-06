import os
import argparse
import subprocess
import numpy as np
from pymatgen.io.cif import CifParser
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import MultiWeightsChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import MultiWeightsChemenvStrategy
from tqdm import tqdm

def adjacent_table1(cif_file,checktag):
	### Setup the local geometry finder ###
	lgf = LocalGeometryFinder()
	lgf.setup_parameters(centering_type='centroid', include_central_site_in_centroid=True)

	try:
		parser = CifParser(cif_file)
		try:
			struct = parser.get_structures()[0]
		except:
			pass
			checktag = False
	except:
		pass

	try:
		lgf.setup_structure(structure=struct)
	except:
		pass
		#print('error lgf.setup_structure(structure)')
		#checktag = ' error lgf.setup_structure(structure) '
		checktag = False
	
	return checktag


def adjacent_table2(cif_file,cifnumname,cifdataout):
	### Setup the local geometry finder ###
	lgf = LocalGeometryFinder()
	lgf.setup_parameters(centering_type='centroid', include_central_site_in_centroid=True)
	
	cif_text=subprocess.getoutput('cat {cif}'.format(cif=cif_file))
	try:
		space_group_IT_number = cif_text.split('_space_group_IT_number')[1].split('\n')[0].strip()
	except:
		space_group_IT_number = None
	
	parser = CifParser(cif_file)
	struct = parser.get_structures()[0]
	lgf.setup_structure(structure=struct)
	
	### Get the StructureEnvironments ###
	try:    
		se = lgf.compute_structure_environments(maximum_distance_factor=1.41,only_cations=False)
	except:
		print('error lgf.compute_structure_environments(maximum_distance_factor=1.41,only_cations=True)')
		#return 0
	
	### Get lightstructure emvironment ###
	try:
		###  strategy  ###   
		strategy = MultiWeightsChemenvStrategy.stats_article_weights_parameters()
		lse = LightStructureEnvironments.from_structure_environments(strategy=strategy,structure_environments=se)
	except:
		print('error LightStructureEnvironments.from_structure_environments(strategy=strategy,structure_environments=se)')
		return 0
	
	### print ###
	with open( cifdataout + "/" + "{0}.txt".format(cifnumname),"w") as file2:
		print(lse.structure,file=file2)
		print(file=file2)
		print('space_group_IT_number',space_group_IT_number,file=file2)
		print('\n*** lse ***\n',file=file2)
		for isite in range(len(lse.structure)):
			print("\n**isite=",isite,file=file2)
			
			#print(se.structure[isite].species_and_occu,":",lse.structure.sites[isite].coords,file=file2)
			print(se.structure[isite].species,":",lse.structure.sites[isite].coords,file=file2)
			
			isite_coords=lse.structure.sites[isite].coords
			print(file=file2)
			if lse.coordination_environments[isite] is None:
				continue 
			for i in range(len(lse.coordination_environments[isite])):
				print('ce_symbol  =',lse.coordination_environments[isite][i]['ce_symbol'],file=file2)
				print('ce_fraction= ',lse.coordination_environments[isite][i]['ce_fraction'],file=file2)
				cn=len(lse.neighbors_sets[isite][i].all_nbs_sites_indices)
				_distance=[]
				print("species:coords:index",file=file2)
				for inbs in range(len(lse.neighbors_sets[isite][i].neighb_sites_and_indices)):      #neib loop                   
					coords=lse.neighbors_sets[isite][i].neighb_sites_and_indices[inbs]['site'].coords
					index=lse.neighbors_sets[isite][i].neighb_sites_and_indices[inbs]['index']
					coords_checkindex=lse.structure.sites[index].coords

					_distance.append(np.linalg.norm(isite_coords - coords))
					species=lse.structure.sites[index]._species
					print(species,":",coords,":",index,file=file2)
				print("species:coords:indexend",file=file2)
				if not set(se.neighbors_sets[isite][cn][0].distances)==set(_distance):
					print("!!!atention to nb coordination!!!",file=file2)
					print("        ",se.neighbors_sets[isite][cn][0].distances,file=file2)
					print("        ",_distance,file=file2)
				print(file=file2)
			print('END\n',file=file2)

	return 0


def main():
	cwd = os.getcwd()

	parser = argparse.ArgumentParser()
	parser.add_argument('--codpath', default='COD/O')
	parser.add_argument('--output1', default='result')
	parser.add_argument('--output2', default='cod')
	parser.add_argument('-rs','--restart', default=False)
	args = parser.parse_args()
	
	if not os.path.isdir(args.output1):
		os.makedirs(args.output1 + "/" + args.output2)
	else:
		if not os.path.isdir(args.output1 + "/" + args.output2):
			os.makedirs(args.output1 + "/" + args.output2)
	
	ciffile = subprocess.getoutput("find {0} -name '*.cif'|sort".format(args.codpath))
	ciffilelist = ciffile.split('\n')
	#print(ciffilelist)
	
	checker = args.restart
	print('make adjacent info by pymatgen')
	for i in tqdm(ciffilelist):
		
		#cifnum = re.findall('\/(\w+)\.cif',i)[0]
		cifnum=os.path.basename(i).replace('.cif','')
		cifdataout = os.path.join(args.output1,args.output2,cifnum)
		outfile = os.path.join(cifdataout,f"{cifnum}.txt")
		if os.path.isfile(outfile):
			continue
		if cifnum == args.restart or not checker:
			checker = True
		
		if checker:
			if adjacent_table1(i,True):
				if not os.path.isdir(cifdataout):
					os.mkdir(cifdataout)
				#subprocess.getoutput("cp {0} {1}".format(i,cifdataout))
			
				#print(cifnum)
				adjacent_table2(i,cifnum,cifdataout)
			
			else:
				print("{0} : error".format(cifnum))
	print('end adjacent info')
if __name__ == '__main__':
	main()


