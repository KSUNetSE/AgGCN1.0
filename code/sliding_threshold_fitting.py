from numpy import linalg as LA
import operator
from scipy.optimize import minimize
import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from sklearn.metrics import r2_score
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy import stats

gr = open('threshold 0.5% interval 10.csv','r')
gene_mean = list(csv.reader(gr,delimiter=','))
gr.close()
gene_mean =list(map(list,zip(*gene_mean)))
x0=[(int(gene_mean[0][i])+ int(gene_mean[1][i]))/2.0 for i in range(len(gene_mean[0]))]
y0=[float(ele) for ele in gene_mean[2]]
#print(x)
#print(y)


def Dis(alpha,eta,lam,beta):
    ssres=0
    sstot=0
    ssmod =0
    ybar=np.mean(y0)
    for i in range(len(x0)):
        ssres += (alpha-1/(eta+lam*np.exp(-x0[i]/beta))-y0[i])**2
        sstot +=(y0[i] -ybar)**2
        ssmod += (alpha-1/(eta+lam*np.exp(-x0[i]/beta))-ybar)**2
    p_value = 1-stats.f.cdf(ssmod/ssres*(len(y0)-4)/(4-1), 4-1, len(y0)-4)
    return 1-ssres/sstot,p_value

lam_series =np.linspace(1.6,6,201)
#print(lam_series)
beta_series =np.linspace(10,40,401)
alpha_series =np.linspace(1.1,1.5,21)
eta_series =np.linspace(1.3,2.2,401)
#L,B =np.meshgrid(lam_series,beta_series)
L,B =np.meshgrid(alpha_series,eta_series)
#Z,p_va =Dis(1.28,1.6,L,B)
Z,p_va =Dis(L,B,2,30)
fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
ax = fig.gca(projection='3d')

surf = ax.plot_surface(L, B, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=True)
#ax.scatter(L,B,Z,s=1)
ax.set_xlabel('lambda')
ax.set_ylabel('beta')
ax.set_zlabel('R2')
ax.zaxis.set_major_locator(LinearLocator(5))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
fig.colorbar(surf, shrink=0.5, aspect=10)
#ax.view_init(45, 45)
idlam=0
idbeta=0
maxr2=0
for i in range(len(Z)):
    #print(Z[i])
    for j in range(len(Z[i])):
        if Z[i][j]> maxr2:
            maxr2=Z[i][j]
            idlam=j
            idbeta=i

#print('R2 =',maxr2,'lambda =',lam_series[idlam],'beta =',beta_series[idbeta])
#R2,p_value=Dis(1.16,1.9,lam_series[idlam],beta_series[idbeta])

print('R2 =',maxr2,'alpha =',alpha_series[idlam],'eta =',eta_series[idbeta])
R2,p_value=Dis(alpha_series[idlam],eta_series[idbeta],4.3,24)

print('R2 =',R2,'p-value =',p_value)








