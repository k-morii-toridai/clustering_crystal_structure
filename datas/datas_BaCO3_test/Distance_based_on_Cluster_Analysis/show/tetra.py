#!/usr/bin/env python3
import sys,re,argparse
import pandas as pd
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d

def tetra(csv1,title='cluster'):
    x=csv1['x']
    y=csv1['y']
    z=csv1['z']
    isite=csv1['isite']
    atom=csv1['atom']
    
 

    #plt.style.use('ggplot')
    plt.rcParams["axes.facecolor"] = 'white'
    fig = plt.figure(figsize=(8, 8)) # 図の設定 
    ax = fig.add_subplot(projection='3d') # 3Dプロットの設定

    poly=list(zip(x.iloc[[5,6,7]],y.iloc[[5,6,7]],z.iloc[[5,6,7]]))
    ax.add_collection3d(art3d.Poly3DCollection([poly],color='#4169e133'))
    poly=list(zip(x.iloc[[6,7,8]],y.iloc[[6,7,8]],z.iloc[[6,7,8]]))
    ax.add_collection3d(art3d.Poly3DCollection([poly],color='#4169e133'))
    poly=list(zip(x.iloc[[5,7,8]],y.iloc[[5,7,8]],z.iloc[[5,7,8]]))
    ax.add_collection3d(art3d.Poly3DCollection([poly],color='#4169e133'))
    poly=list(zip(x.iloc[[5,6,8]],y.iloc[[5,6,8]],z.iloc[[5,6,8]]))
    ax.add_collection3d(art3d.Poly3DCollection([poly],color='#4169e133'))
    

    #for i in range(len(x)):
        #ax.scatter(x[i], y[i], z[i], label='(x, y, z)',c="blue") # 始点
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title(title)
    ax.set_xlim(-5, 5)
    ax.set_ylim(-5, 5)
    ax.set_zlim(-5, 5)
    for i in range(len(x)):
        text=str(atom[i])+'_'+str(isite[i])
        #text=str(isite[0])+'_'+str(isite[i])
		#text='ABW_'+str(isite[0])+'_'+str(isite[i])
        ax.text(x[i],y[i],z[i],text,size=8)

    noods=list()
    for index,i in csv1.iterrows():
        if index==0:
            continue
        front_idx=i.loc['front_index']
        a=csv1.loc[front_idx].loc['x':'z']
        b=i.loc['x':'z']
        noods.append(([a.x,b.x],[a.y,b.y],[a.z,b.z]))
    for nood in noods:
        line = art3d.Line3D(*nood)
        ax.add_line(line)

    plt.show()

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
    tetra(df,title=showname)

if __name__ == '__main__':
	main()
