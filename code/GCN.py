'''***************************************************************
Caterina Scoglio and Kristin Michel
Kansas State University
Last Modified: Sept 2022
Copyright (c) 2022, Caterina Scoglio and Kristin Michel. All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted
********************************************************************'''
import dataload
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--data', type=str, default='anopheles.csv', help='A .csv matrix where rows represent nodes and columns represent experiments')
parser.add_argument('--zero_removed',  type=str, default=True, help='default:True, False: 0s will be removed from the matrix')
parser.add_argument('--zscored', type=str, default=True, help='default:True, False: z-score normalize each column')
parser.add_argument('--rescaled', type=str, default=True, help='default:True, False: rescale each value with log2')
parser.add_argument('--save_data', type=str, default=True, help='default:True, save the normalzied data')
parser.add_argument('--curve_params', type=str, default=False, help='default: False, use pre-determined parameters for the curve, otherwise the params will be optimized')
parser.add_argument('--alpha', type=float, default=1.0, help='default: 1.0, use this pre-determined alpha if curve_params is true')
parser.add_argument('--eta', type=float, default=1.5, help='default:1.5, use this pre-determined eta if curve_params is true')
parser.add_argument('--lam', type=float, default=2, help='default:2, use this pre-determined lam if curve_params is true')
parser.add_argument('--beta', type=float, default=30, help='default:30, use this pre-determined beta if curve_params is true')
parser.add_argument('--dropna', type=int, default=200, help='default:200, keep only the columns with at least 200 values')
parser.add_argument('--index2gene', type=str, default='', help='defalut:'', save a list which map each gene with an integer')
parser.add_argument('--pcc', type=str, default='PCC_matrix.csv', help='default: PCC_matrix.csv, save the PCC matrix if the input is not empty')
parser.add_argument('--paired_elements', type=str, default='paired_elements.csv', help='default: paired_elements.csv, save the paired_elements if the input is not empty')
parser.add_argument('--bin_size', type=int, default=10, help='default:10, binning the PCC into different intervals per paired element')
parser.add_argument('--cutoff', type=float, default=0.005, help='default:0.005, choose the top fraction of PCC as edges of the network')
parser.add_argument('--edgelist', type=str, default='edgelist.csv', help='default:edgelist.csv, save the edge list')
parser.add_argument('--thres_curve', type=str, default='threshold_curve.png', help='default:threshold_curve.png, save the fitted sliding threshold curve')


args = parser.parse_args()

if __name__ == "__main__":
    gene = dataload.GeneCoExp( filepath = args.data)
    
    gene.read_data( zero_removed = args.zero_removed, 
                   zscored = args.zscored, 
                   rescaled = args.rescaled, 
                   dropna = args.dropna,
                   savefile = args.save_data)
    
    gene.Calculate_PCC(index2gene =args.index2gene, 
                       pcc_path = args.pcc, 
                       paired_elements_path = args.paired_elements)
    
    gene.calculate_threshold(bin_size=args.bin_size, cutoff=args.cutoff)
    gene.curve_fitting(alpha=args.alpha,eta=args.eta,lam=args.lam,beta=args.beta, #assign initial values to the 4 parameters
                       para_in= args.curve_params, 
                       edgelist = args.edgelist, 
                       thresholdcurve = args.thres_curve)


