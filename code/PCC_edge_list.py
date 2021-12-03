import numpy as np
import csv
import sys
from scipy.stats import pearsonr
gr = open('all condition.csv','r')
gene_mean0 = list(csv.reader(gr,delimiter=','))
gr.close()

sd =int(sys.argv[1])*1000
sd0=(int(sys.argv[1])-1)*1000

gene=gene_mean0[5:]
        
gene_new=[]
gene_id=[]
gene_name=[]
gene_order=[]
gene_label=[]
for i in range(5,len(gene_mean0)):
    gene_order.append(int(gene_mean0[i][0]))
    gene_name.append(gene_mean0[i][1])

    #gene_id.append(x)
for i in range(len(gene)):
    gene[i].pop(0)
    gene[i].pop(0)
print(gene[4])
def findsameitem(xc0,y0):
    x =[]
    y =[]
    for index,item in enumerate(xc0):
        if item and y0[index]:
            x.append(float(item))
            y.append(float(y0[index]))
    return x,y
Ap=[]
Bonfe=0
basic=[]
apv=[]

if sd==10000:
    sd=len(gene)-1
for j in range(sd0,sd): 
    for i in range(j+1,len(gene)):
        x,y =findsameitem(gene[j],gene[i])
            #Corr =PearsonCorr(x,y)
        lenx=len(x)
        if lenx>4:
            corr,p_value =pearsonr(x,y)
            #Ap.append((gene_name[j],gene_label[j],gene_name[i],gene_label[i],corr,p_value,lenx))
            Ap.append((gene_order[j],gene_order[i],round(corr,3),lenx,p_value))
            apv.append((round(corr,3),lenx))
            
            if lenx<140:
                Bonfe+=1
lenAp=len(Ap)

boncorr1=[Bonfe,lenAp-Bonfe,lenAp-Bonfe]
boncorr=[0.05/ele for ele in boncorr1]


gi = open('pcc '+str(sd0)+'.csv','w',newline='')
cw = csv.writer(gi,delimiter=',')
cw.writerows(apv)
gi.close()

edl=[]
for i in range(lenAp):
    if float(Ap[i][4])<boncorr[int(Ap[i][3])//140]:# or (float(gene_mean[i][4])> 0.7 and float(gene_mean[i][5])<0.05):        
        if 0.99>float(Ap[i][2])>(1.34-1/(1.47+1.732*np.exp(-float(Ap[i][3])/38.5))):            
                #edl.append((gene_mean[i][0],gene_mean[i][2],(float(gene_mean[i][6]))*float(gene_mean[i][4])))
            edl.append((Ap[i][0],Ap[i][1],float(Ap[i][2]),int(Ap[i][3])))
# edl=[]
# for i in range(lenAp):
#     if float(Ap[i][4])<boncorr:# or (float(gene_mean[i][4])> 0.7 and float(gene_mean[i][5])<0.05):        
#         if 0.99>float(Ap[i][2])>(1.14-1/(1.9+4.3*np.exp(-float(Ap[i][3])/23.7))):            
#                 #edl.append((gene_mean[i][0],gene_mean[i][2],(float(gene_mean[i][6]))*float(gene_mean[i][4])))
#             edl.append((Ap[i][0],Ap[i][1],float(Ap[i][2]),int(Ap[i][3])))

gi = open('edges/edge'+str(sd0)+'.csv','w',newline='')
cw = csv.writer(gi,delimiter=',')
cw.writerows(edl)
gi.close()
edl=[]

# gi = open('results1_1/'+str(sd)+' removed edgelist.csv','w',newline='')
# cw = csv.writer(gi,delimiter=',')
# cw.writerows(edl)
# gi.close()
# edl=[]
