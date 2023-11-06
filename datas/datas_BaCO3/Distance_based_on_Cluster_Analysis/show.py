import re,os
import pandas as pd
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d
from .distance import cal_distances
from copy import deepcopy
from glob import glob
import os
def change(csv2,csv1,show=True,save=False,hist=False,text=True):
    #print(os.getcwd())
    #sys.exit()
    csvaddress1=csv1
    csvaddress2=csv2
    csv1=pd.read_csv(csv1,index_col=0)
    csv2=pd.read_csv(csv2,index_col=0)
    x=csv1['x']
    y=csv1['y']
    z=csv1['z']
    isite=csv1['isite']
    atom=csv1['atom']
    
    u=csv2['x']
    v=csv2['y']
    w=csv2['z']
    isite2=csv2['isite']
    atom2=csv2['atom']
    
    val=cal_distance(csvaddress1,csvaddress2,values=True)
    u2=[]
    v2=[]
    w2=[]
    isite3=[]
    atom3=[]
    for i in range(len(val)):
        a=re.sub(r"\D","",val[i][0][-2:])
        u2.append(u[int(a)])
        v2.append(v[int(a)])
        w2.append(w[int(a)])
        isite3.append(isite2[int(a)])
        atom3.append(atom2[int(a)])

    #plt.style.use('ggplot')
    plt.rcParams["axes.facecolor"] = 'white'
    fig = plt.figure(figsize=(8, 8)) # 図の設定 
    ax = fig.add_subplot(projection='3d') # 3Dプロットの設定
    for i in range(len(x)):
        ax.quiver(x[i], y[i], z[i], u2[i]-x[i], v2[i]-y[i], w2[i]-z[i], arrow_length_ratio=0.1,color="red") # 矢印プロット
        ax.scatter(x[i], y[i], z[i], label='(x, y, z)',c="blue") # 始点
        ax.scatter(u2[i], v2[i], w2[i], label='(x+u, y+v, z+w)',c="green") # 終点
		#ax.set_xlabel('x')
		#ax.set_ylabel('y')
		#ax.set_zlabel('z')
        ax.set_title(os.path.basename(csvaddress1).replace('.csv','')+'(blue) to ' +os.path.basename(csvaddress2).replace('.csv','')+'(green)',size=10) # タイトル
        #ax.set_xlim(-5, 5)
        #ax.set_ylim(-5, 5)
        #ax.set_zlim(-5, 5)
    parlist=list()
    if text:
        for i in range(len(x)):
            text=str(atom[i])+'_'+str(isite[i])
            #text=str(isite[0])+'_'+str(isite[i])
            #text='ABW_'+str(isite[0])+'_'+str(isite[i])
            text2=str(atom3[i])+'_'+str(isite3[i])
            #text2=str(isite2[0])+'_'+str(isite3[i])
            #text2='ABW_'+str(isite2[0])+'_'+str(isite3[i])
            ax.text(x[i],y[i],z[i],text,size=8)
            ax.text(u2[i],v2[i],w2[i],text2,size=8)
            parlist.append(('{}_{}'.format(text,text2)))

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
    if show:
        plt.show()
    if save:
        plt.savefig('change.svg')
    plt.close()
    if hist:
        df=pd.Series(cal_distance(csvaddress1,csvaddress2,histgram=True),index=parlist)
        return df


def double_clusterplot(cluster_df1,cluster_df2,title='cluster.png',show=None,save=True):
    noods=list()
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(-5,5)
    ax.set_ylim(-5,5)
    ax.set_zlim(-5,5)
    for cluster_df in [cluster_df1,cluster_df2]:
        for index,i in cluster_df.iterrows():
            if index==0:
                continue
            front_idx=i.loc['front_index']
            a=cluster_df.loc[front_idx].loc['x':'z']
            b=i.loc['x':'z']
            noods.append(([a.x,b.x],[a.y,b.y],[a.z,b.z]))
        ax.scatter(cluster_df.x,cluster_df.y,cluster_df.z)
        for index,i in cluster_df.iterrows():
            text=i.atom+'_'+str(int(i.isite))
            ax.text(i.x,i.y,i.z,text)
        for nood in noods:
            line = art3d.Line3D(*nood)
            ax.add_line(line)
        fig.suptitle(title)
        if save:
            fig.savefig(title)
        if show:
            plt.show()
        plt.close()
from PIL import Image
class DrawGif():
    def __init__(self):
        self.image=list()
    def set_data(self,cluster_df1,cluster_df2):
        pngname='cluster.png'
        double_clusterplot(cluster_df1,cluster_df2,title=pngname)
        self.image.append(Image.open(pngname))
        os.remove(pngname)
        return
    def makegif(self,filename='clustergif.gif'):
        self.image[0].save(filename,save_all=True, append_images=self.image[1:],optimize=False, duration=500, loop=0)
        return

def emd_histgram(cifdir,database_address='database',show=False,save=True):
    distanlist=glob('{}/*distance'.format(cifdir))
    histdf=pd.DataFrame()
    for dataaddress in distanlist:
        data=pd.read_csv(dataaddress,index_col=0).iloc[0]
        clusteraddress_base='{}/{}_{}.csv'.format(database_address,data.isite_j,data.pattern_j)
        clusteraddress_cif='{}/{}_{}.csv'.format(cifdir,data.isite_i,data.pattern_i)
        d=cal_distance(csv_address1=clusteraddress_cif,csv_address2=clusteraddress_base,histgram=True)
        d.sort(reverse=True)
        histdf.loc[:,'{}_{}'.format(data.isite_i,data.isite_j)]=deepcopy(d)
    histdf=histdf.sort_index(axis=1)
    histdf.plot(kind='bar')
    if show:
        plt.show()
    if save:
        plt.savefig('{}/{}_move.png'.format(cifdir,data.isite_i))
    return histdf


def clusterplot(cluster_df,text=True,color=False):
    #set nood
    noods=list()
    for index,i in cluster_df.iterrows():
        if index==0:
            continue
        front_idx=i.loc['front_index']
        a=cluster_df.loc[front_idx].loc['x':'z']
        b=i.loc['x':'z']
        noods.append(([a.x,b.x],[a.y,b.y],[a.z,b.z]))
    
    fig = plt.figure(figsize = (6, 6))
    ax = fig.add_subplot(111, projection='3d')
    '''
    ax.set_xlim(-5,5)
    ax.set_ylim(-5,5)
    ax.set_zlim(-5,5)
    '''
    if text:
        for index,i in cluster_df.iterrows():
            text=i.atom+'_'+str(i.isite)
            ax.text(i.x,i.y,i.z,text)
    
    if not color is False:
        for atom,color_ in color.items():
            specific_atom=cluster_df[cluster_df.atom==atom].copy()
            ax.scatter(specific_atom.x,specific_atom.y,specific_atom.z,color=color_,label=atom)
    else:
        ax.scatter(cluster_df.x,cluster_df.y,cluster_df.z)

    for nood in noods:
        line = art3d.Line3D(*nood)
        ax.add_line(line)