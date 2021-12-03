import sys
import numpy as np
import matplotlib.pyplot as plt
import csv
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import pairwise_distances
import seaborn as sn
import networkx as nx
from collections import defaultdict,Counter

plt.style.use('seaborn-whitegrid')

nw = open('new_network/LCC network rescaled.csv','r')
nwl = list(csv.reader(nw,delimiter=','))
nw.close()

nl=[]
for ele in nwl:
    if ele[0] not in nl:
        nl.append(ele[0])
    if ele[1] not in nl:
        nl.append(ele[1])
        
N = len(nl)
G=nx.Graph()
for ele in nl:
    G.add_node(ele)

for ele in nwl:
    G.add_edge(ele[0],ele[1],weight=round(float(ele[2]),2))
    
def rn(G,thres):
    for ele in list(G.nodes()):
        if G.degree(ele,weight='weight')<=thres:
            G.remove_node(ele)
    return G

numnode=[]
strth=[]
core=[['Id','core7464']] #74.52 74.64
while len(G.nodes())>0:
    deg=[]
    for ele in list(G.nodes()):
        deg.append(G.degree(ele,weight='weight'))
    strth.append(min(deg))
    numnode.append(len(G.nodes()))
    thr=min(deg)
    nl_1=len(G.nodes())
    G=rn(G,thr)
    if thr>74.52 and thr<74.65:
        for ele in list(G.nodes()):
            core.append([ele,'core7464'])
    while len(G.nodes())<nl_1:
        nl_1=len(G.nodes())
        G=rn(G,thr)
    
print(len(core))
gi = open('core74.64.csv','w',newline='')
cw = csv.writer(gi,delimiter=',')
cw.writerows(core)
gi.close()  
plt.scatter(strth,numnode,s=5,color='b',marker='^')
plt.xlabel('s-core threshold')
plt.ylabel('number of nodes')
plt.savefig('core.pdf')  
