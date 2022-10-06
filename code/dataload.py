import numpy as np
import pandas as pd
import csv
from collections import defaultdict
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


class GeneCoExp():
    def __init__(self, corr='', paired_items='', filepath='expression.csv'):
        self.filepath = filepath
        self.N = 0
        self.M = 0
        self.corr = corr
        if corr and paired_items:
            print('reading correlation table...')
            self.rho = pd.read_csv(corr, delimiter=',', index_col=0)
            self.paired_items = pd.read_csv(paired_items, delimiter=',', index_col=0)
            self.rho.columns = [int(ele) for ele in self.rho.columns]
            self.paired_items.columns = [int(ele) for ele in self.paired_items.columns]
        
    def read_data(self, zero_removed=False, zscored=False,  rescaled=False, dropna=0, savefile=''):
        '''
        Parameters
        ----------
        zero_removed : required
            remove zeros from the table. The default is False.
        zscored : required
            z-score normalize each column. The default is False.
        rescaled : required
            log2 rescale each value. The default is False.
        savefile : optional
            save the normalized data. The default is ''.
        '''
        print('normalizing data...')
        self.data = pd.read_csv(self.filepath, delimiter=',', index_col=0)
        if self.N != len(self.data): self.N = len(self.data)
        if self.M != len(self.data.columns): self.M = len(self.data.columns)
        
        if not zero_removed:
            self.data.replace(0, np.nan, inplace=True)#remove 0s 

        if not zscored:  # zscore normalize the data
            if not rescaled: #rescale the expression data, if input is not rescaled
                for col in self.data.columns:
                    self.data[col] = np.log2(self.data[col])
                    
            for col in self.data.columns: # zscore normalize the data
                self.data[col] = (self.data[col] - self.data[col].mean())/self.data[col].std(ddof=0)    
                #print(self.data[col].mean(), self.data[col].var())
                assert abs(self.data[col].mean())<1e-1 and abs(self.data[col].var()-1.)<1e-1
        
        self.data = self.data.dropna(thresh = dropna, axis='columns')
        if savefile:
            self.data.to_csv('z-scored '+self.filepath, sep=',')
            print('normalized data are saved as: '+'z-scored '+self.filepath)
        
    def Calculate_PCC(self, min_periods=4, index2gene ='index2gene.csv', pcc_path='', paired_elements_path=''):
        '''
        Parameters
        ----------
        min_periods : int, optional
            Minimum number of observations required per pair of columns to have a valid result. The default is 4.
        gene2index : string, required
            corresponds each gene with a number. The default is 'gene2index.csv'.
        pcc_path : string.csv, optional
            save the PCC matrix table. The default is ''.
        paired_elements_path : string.csv, optional
            save the number of observations between each pair of columns. The default is ''.
        '''
        self.min_periods = min_periods #the min # of paired elements
        if self.corr:
            pass
        else:        
            #calculate PCC
            genelist = list(self.data.index)
            self.index2gene = {idx: ele for idx,ele in enumerate(genelist)}
            if index2gene:
                (pd.DataFrame.from_dict(data=self.index2gene, orient='index').to_csv(index2gene, header=False))
                print('index 2 gene table is saved as: '+ index2gene)
            self.data.index = list(range(len(genelist)))
            
            print('calculating PCC...')
            self.rho = self.data.T.astype('float16').corr(min_periods=self.min_periods)   #compute the correlation coefficient
            self.rho = self.rho.round(4)
            self.rho = pd.DataFrame(np.triu(self.rho.to_numpy(), 1))
            #save the PCC matrix and paired elements matrix
            if pcc_path:
                self.rho.to_csv(pcc_path, sep=',')
                print('PCC table is saved as: ' + pcc_path)
            print('determining the # of paired elements for every node pair...')
            #find the # of paired elements between every pair of genes
            val_df = self.data.notnull().astype('int').to_numpy() #replace non-empty with 1, and empty with 0    
            self.paired_items = np.dot(val_df, val_df.T) #(N, M) * (M, N) -> (N, N)
            self.paired_items = pd.DataFrame(np.triu(self.paired_items, 1))
    
            #save the paired elements matrix
            if paired_elements_path:
                self.paired_items.to_csv(paired_elements_path, sep =',')
                print('paired element table is saved as: ' + paired_elements_path)
        
    def calculate_threshold(self, bin_size=10, cutoff=0.005):
        '''
        Parameters
        ----------
        bin_size : 10 or 5, required
            divide the PCC into different intervals per the # of paired elements. The default is 10.
        cutoff : float, required
            choose the top fraction PCC to construct network. The default is 0.005.
        '''
        #self.rho, self.paired_items
        print('calculating threshold...')
        pair2pcc = defaultdict(list) #{4:[0.3, 0.5,...], 5:[0.1, 0.5, ...],...}
        triu = np.triu_indices(len(self.paired_items),1) # the upper triangular matrix
        pair_triu = list(self.paired_items.to_numpy()[triu]) #the upper triangular of paired elements matrix 
        pcc_triu =  list(self.rho.to_numpy()[triu]) #the upper triangluar of the PCC matrix
        for idx in range(len(pair_triu)): #
            if pair_triu[idx] >= self.min_periods:
                pair2pcc[pair_triu[idx]].append(pcc_triu[idx]) #
        
        pairs = sorted(list(pair2pcc.keys())) 
        print(pairs)
        bin2pair = defaultdict() # {4:7, 5:7,...,101:105, 102:105}
        for p in pairs:
            bin2pair[p] = ((p-1)//bin_size)*bin_size + (bin_size+1)/2.
            if p <=10:
                bin2pair[p] = (min(pairs)+10)/2                         
            if max(pairs)//bin_size==((p-1)//bin_size):
                bin2pair[p] = (((p-1)//bin_size)*bin_size + max(pairs)+1)/2.
        
        bin2pcc = defaultdict(list) # {7:[pcc], 15:[pcc], 25:[pcc]}
        for key in pairs:
            if len(pair2pcc[key])>0:
                bin2pcc[bin2pair[key]] += pair2pcc[key]

        self.bin2cut = {}
        for key,val in bin2pcc.items():
            cut_val = sorted(val)[int(len(val)*(1-cutoff))]
            if len(val)>50 and cut_val > 0.3:
                self.bin2cut[key] = cut_val
        print('threshold of each interval: ', self.bin2cut)
        
    def curve_fitting(self, alpha=1,eta=1.5,lam=2,beta=30, para_in= False, edgelist='edgelist.csv', thresholdcurve='threshold_curve.png'):
        '''
        Parameters
        ----------
        p0 : list, required
            the initial values of the 4 parameter to fit the curve. The default is [1,1.5,2,30].
        edgelist : string, required
            save the edge list. The default is 'edgelist.csv'.
        thresholdcurve : string, optional
            save the fitted threhold curve. The default is 'threshold_curve.png'.
        '''
        if len(self.bin2cut)<4:
            thresh = np.mean([val for key,val in self.bin2cut.items()])
            el = np.where(self.rho.to_numpy()>=thresh, self.rho, 0)
        else:
            p0=[alpha,eta,lam,beta]
            xdata=[]
            ydata=[]
            for key,val in self.bin2cut.items():
                xdata.append(key)
                ydata.append(val)
            xdata = np.array(xdata,dtype=np.float64)
            ydata = np.array(ydata,dtype=np.float64)
            popt, pcov =  curve_fit(f=func, xdata=xdata, ydata=ydata, p0=p0, maxfev=20000)
            
            if para_in:
                popt = np.array([alpha,eta,lam,beta])
                
            if thresholdcurve:
                plt.scatter(xdata, ydata)
                plt.plot(xdata, func(xdata, *popt))
                plt.savefig(thresholdcurve)
    
            thres_max = func(min(xdata), *popt)
            thres_min = func(max(xdata), *popt)
            print('parameters of the curve:')
            print('    alpha: ',popt[0])
            print('    eta: ',popt[1])
            print('    lambda: ',popt[2])
            print('    beta: ',popt[3])
            print('generating edge list...')
            thres_matrix = func(self.paired_items.to_numpy(), *popt)
            el = np.where(self.rho.to_numpy()>=thres_matrix, self.rho, 0)
            
            thres_matrix = thres_max + thres_min - thres_matrix
            el = np.multiply(el, thres_matrix)
        el = np.round(el,decimals = 3)
        
        el_index = np.nonzero(el)
        source =[self.index2gene[ele] for ele in el_index[0]]
        target = [self.index2gene[ele] for ele in el_index[1]]
        
        el = map(list, zip(*[source,target,el[el_index]]))

        el = pd.DataFrame(el, columns=['source','target','weight'])
        el.to_csv(edgelist, sep=',', index=False)
        print('edge list is saved as: ' + edgelist)


def func(x, alpha, eta, lam, beta):
    return alpha - 1/(eta+lam*np.exp(-x/beta))



