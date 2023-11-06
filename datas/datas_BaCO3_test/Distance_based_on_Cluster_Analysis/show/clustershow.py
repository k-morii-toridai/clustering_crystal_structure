#!/usr/bin/env python3
import sys,re,argparse
import pandas as pd
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d


def clusterplot(cluster_df,title='cluster'):
	noods=list()
	for index,i in cluster_df.iterrows():
		if index==0:
			continue
		front_idx=i.loc['front_index']
		a=cluster_df.loc[front_idx].loc['x':'z']
		b=i.loc['x':'z']
		noods.append(([a.x,b.x],[a.y,b.y],[a.z,b.z]))
	fig = plt.figure(figsize = (12, 12))
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(cluster_df.x,cluster_df.y,cluster_df.z)
	ax.set_xlim(-5,5)
	ax.set_ylim(-5,5)
	ax.set_zlim(-5,5)
	for index,i in cluster_df.iterrows():
		text=i.atom+'_'+str(i.isite)
		ax.text(i.x,i.y,i.z,text)
	for nood in noods:
		line = art3d.Line3D(*nood)
		ax.add_line(line)
	fig.suptitle(title)
	plt.show()
	plt.close()
	


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('csvn')
	parser.add_argument('-e','--explanation', default=False)
	args = parser.parse_args()

	if args.explanation:
		print('''-csvn : /cluster_0_0.csv''')
		sys.exit()
	df=pd.read_csv(args.csvn,index_col=0)
	showname=re.split('/',args.csvn)[-1].replace('.csv','')
	clusterplot(df,title=showname)

if __name__ == '__main__':
	main()
