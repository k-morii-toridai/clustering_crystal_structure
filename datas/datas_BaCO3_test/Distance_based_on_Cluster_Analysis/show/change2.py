import sys,re,argparse,os
import pandas as pd
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d
from ..distance_func import cal_distance

def change(csv1,csv2,title='cluster',title2='cluster'):
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
    
    val=cal_distance(title+'.csv',title2+'.csv',values=True)
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
        ax.set_title(title+'(blue) to '+title2+'(green)',size=10) # タイトル
        ax.set_xlim(-5, 5)
        ax.set_ylim(-5, 5)
        ax.set_zlim(-5, 5)
    for i in range(len(x)):
        text=str(atom[i])+'_'+str(isite[i])
        #text=str(isite[0])+'_'+str(isite[i])
		#text='ABW_'+str(isite[0])+'_'+str(isite[i])
        text2=str(atom3[i])+'_'+str(isite3[i])
		#text2=str(isite2[0])+'_'+str(isite3[i])
		#text2='ABW_'+str(isite2[0])+'_'+str(isite3[i])
        ax.text(x[i],y[i],z[i],text,size=8)
        ax.text(u2[i],v2[i],w2[i],text2,size=8)

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

#def main():
    #parser = argparse.ArgumentParser()
    #parser.add_argument('csvn')
    # parser.add_argument('csvn2')
    # parser.add_argument('-e','--explanation', default=False)
    # args = parser.parse_args()
    # if args.explanation:
      # print('''-csvn : /cluster_0_0.csv''')
      # sys.exit()
    # df=pd.read_csv(args.csvn,index_col=0)
    # df2=pd.read_csv(args.csvn2,index_col=0)
    # showname=re.split('/',args.csvn)[-1].replace('.csv','')
    # showname2=re.split('/',args.csvn2)[-1].replace('.csv','')
    # change(df,df2,title=showname,title2=showname2)

#if __name__ == '__main__':
	#main()