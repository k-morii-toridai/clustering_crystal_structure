import csv
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d
from distance_func import cal_distance
import re

csv1='result/testcif/ABW/ABW_0_0.csv'
csv2='result/testcif/ABW/ABW_2_0.csv'


xlist=[]
ylist=[]
zlist=[]
isitelist=[]
atomlist=[]
nextlist=[]
with open(csv1,newline='') as inputfile:
    for row in csv.reader(inputfile):
        if not row[0].isdecimal():
            continue
        xlist.append(float(row[4]))
        ylist.append(float(row[5]))
        zlist.append(float(row[6]))
        isitelist.append(row[2])
        atomlist.append(row[3])
        nextlist.append(row[1])

ulist=[]
vlist=[]
wlist=[]
isitelist2=[]
atomlist2=[]
with open(csv2,newline='') as inputfile:
    for row in csv.reader(inputfile):
        if not row[0].isdecimal():
            continue
        ulist.append(float(row[4]))
        vlist.append(float(row[5]))
        wlist.append(float(row[6]))
        isitelist2.append(int(row[2]))
        atomlist2.append(row[3])



val=cal_distance(csv1,csv2,values=True)

ulist2=[]
vlist2=[]
wlist2=[]
isitelist3=[]
atomlist3=[]
alist=[]
for i in range(len(val)):
    a=re.sub(r"\D","",val[i][0][-2:])
    alist.append(int(a))
    ulist2.append(ulist[int(a)])
    vlist2.append(vlist[int(a)])
    wlist2.append(wlist[int(a)])
    isitelist3.append(isitelist2[int(a)])
    atomlist3.append(atomlist[int(a)])

plt.style.use('ggplot')
plt.rcParams["axes.facecolor"] = 'white'
fig = plt.figure(figsize=(8, 8)) # 図の設定
ax = fig.add_subplot(projection='3d') # 3Dプロットの設定
for i in range(len(xlist)):
    ax.quiver(xlist[i], ylist[i], zlist[i], ulist2[i]-xlist[i], vlist2[i]-ylist[i], wlist2[i]-zlist[i], arrow_length_ratio=0.1) # 矢印プロット
    ax.scatter(xlist[i], ylist[i], zlist[i], label='(x, y, z)',c="blue") # 始点
    ax.scatter(ulist2[i], vlist2[i], wlist2[i], label='(x+u, y+v, z+w)',c="green") # 終点
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.set_zlabel('z')
ax.set_title(csv1+'(blue) to '+csv2+'(green)',size=10) # タイトル
#ax.legend() #　凡例
ax.set_xlim(-5, 5)
ax.set_ylim(-5, 5)
ax.set_zlim(-5, 5)
for i in range(len(xlist)):
    text=str(atomlist[i])+'_'+str(isitelist[i])
    #text=str(isitelist[0])+'_'+str(isitelist[i])
    #text='ABW_'+str(isitelist[0])+'_'+str(isitelist[i])
    text2=str(atomlist3[i])+'_'+str(isitelist3[i])
    #text2=str(isitelist2[0])+'_'+str(isitelist3[i])
    #text2='ABW_'+str(isitelist2[0])+'_'+str(isitelist3[i])
    ax.text(xlist[i],ylist[i],zlist[i],text,size=8)
    ax.text(ulist2[i],vlist2[i],wlist2[i],text2,size=8)
plt.show()