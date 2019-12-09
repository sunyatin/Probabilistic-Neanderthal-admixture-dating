# -*- coding: utf-8 -*-
"""
@author: Windows
# 21 oct 2019: created the script
"""

import numpy as np
import time
import argparse
import sys
from scipy.spatial.distance import cdist

## NOTE HERE THAT WE PRINT THE RIGHT BOUND ALSO, AS IN MOORJANI


np.seterr(divide='ignore', invalid='ignore')

parser = argparse.ArgumentParser(description='Calculate covariance in genotypes.')
parser.add_argument('-f', '--filePrefix', type=str, required=True, help='Prefix of the files to analyze.')
parser.add_argument('-p', '--targetPopulation', type=str, required=True, help='Name of the target population to analyze.')
parser.add_argument('-o', '--outfilePrefix', type=str, required=True, help='Prefix of the output file.')
parser.add_argument('-minD', '--minD', type=float, default=0.01, help='Minimum genetic distance.')
parser.add_argument('-maxD', '--maxD', type=float, default=20.0, help='Maximum genetic distance.')
parser.add_argument('-stepD', '--stepD', type=float, default=0.1, help='Bin size.')
parser.add_argument('--Morgans', action='store_true', default=False, help='Add this option if genetic distance are in Morgans (by default assumes centiMorgans).')
parser.add_argument('--chrom', type=int, required=False, default=None, help='Add this option to run the analysis on a particular chromosome.')

args = parser.parse_args()

import os
#os.chdir('/home/rtournebize/Desktop/postdoc/benpeter/covariance')
'''
input_prefix = 'Afon3'
target_popname = 'Afon3'
output_prefix = 'Afon3'
stepD_cM = 0.1
minD_cM = 0
maxD_cM = 10
input_distance_in_cM = False
chrom_to_analyze = cta = None
'''

input_prefix = args.filePrefix
target_popname = args.targetPopulation
output_prefix = args.outfilePrefix
stepD_cM = args.stepD
minD_cM = args.minD
maxD_cM = args.maxD
input_distance_in_cM = not args.Morgans
chrom_to_analyze = args.chrom

print('\n\n')
print(input_prefix)
print(output_prefix)
print('Target pop: '+target_popname)
print('D:'+str(minD_cM)+' '+str(maxD_cM)+' by '+str(stepD_cM))
print('Input distance in cM: '+str(input_distance_in_cM))
print('===============================================================')


###
#####
###

start = time.time()
print('>   '+input_prefix)
print('Analyzing only '+target_popname+'.')
corr_decimals = 1e-10

D_FULL = np.genfromtxt(input_prefix+'.snp', dtype=float, usecols = (1,2))
CHR = D_FULL[:,0].astype(int)
D_FULL = D_FULL[:,1]
if input_distance_in_cM is False:
    print('Converting Morgans input into centiMorgans.')
    D_FULL = 100.0 * D_FULL
POP = np.genfromtxt(input_prefix+'.ind', dtype=str, usecols = 2)
if len(POP.shape) == 0:
    POP = np.array(str(POP))

if chrom_to_analyze is None:
    chr_values = np.unique(CHR).tolist()
else:
    chr_values = np.array([chrom_to_analyze])
print('Chromosomes: '+' '.join([str(x) for x in chr_values]))

RR, first = [], True
for chrom in chr_values:
    foc = np.where(CHR == chrom)[0]

    r0, r1 = foc[0], foc[-1]  
    print('>> chrom:    '+str(chrom)+'   -~-   Range: '+str(r0)+' => '+str(r1))
    D = D_FULL[foc]

    del foc
    G = np.genfromtxt(input_prefix+'.geno', delimiter=[1]*len(np.atleast_1d(POP)), dtype=int,
                           skip_header = r0, max_rows = r1-r0+1)
    
    GL = np.genfromtxt(input_prefix+'.glhood', skip_header = r0, max_rows = r1-r0+1)
    
    if target_popname not in POP:
        sys.exit("ERROR! Pop1 not in list!")

    G = G[:,POP==target_popname]
    GL = GL[:,POP==target_popname]
    if G.shape[1] != 1:
        sys.stop('You did not provide a single sample specification.')

    nSNP_raw = G.shape[0]

    if maxD_cM is None:
        maxD_cM = D[len(D)]
    
    nSNP = G.shape[0]
    print('Raw nSNP: '+str(nSNP_raw))
    
    print('ASSUMES NO MISSING DATA')

    bins_left_bound = np.arange(minD_cM, maxD_cM, step=stepD_cM)
    n_bins = len(bins_left_bound)
    print('There will be '+str(n_bins)+' bins.')
    
    '''
    mD = ( (abs(D[:,None]-D) - minD_cM + corr_decimals) // stepD_cM).astype(np.int32)
    mD[np.triu_indices(mD.shape[0], 0)] = -1
    
    COV, N = [], []
    for bbin in range(len(bins_left_bound)):
        ind = np.where(mD==bbin)
        COV += [np.cov(G[ind[0],0], G[ind[1],0])[0][1]]
        N += [len(ind[0])]
    R = np.column_stack(([chrom]*len(bins_left_bound),
                   bins_left_bound,
                   bins_left_bound + stepD_cM/float(2),
                   np.array(COV),
                   np.array(N)))
    '''

    #BD = {key: [[],[]] for key in bins_left_bound}
    BINS1 = {key: [] for key in range(len(bins_left_bound))}
    BINS2 = {key: [] for key in range(len(bins_left_bound))}
    WEIGHTS = {key: [] for key in range(len(bins_left_bound))}

    for i in range(nSNP):

        gi = G[i,0]
        if gi == 9:
            sys.stop('there is a missing value!')
        gs = G[(i+1):nSNP,0]
        
        gli = GL[i,0]
        gls = GL[(i+1):nSNP,0]

        dd = D[(i+1):nSNP]
        d = D[i]
        d_ = abs(dd - d)
        bins = ( (d_ - minD_cM + corr_decimals) // stepD_cM).astype(np.int32)
        del d, d_, dd
        
        goods_Dref = np.where((0<=bins) & (bins<=(n_bins-1)))[0]
        bins = bins[goods_Dref]
        gs = gs[goods_Dref]
        gls = gls[goods_Dref]
        # check this in ASCEND!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        for bbin in range(len(bins_left_bound)):
            foc = np.where(bins==bbin)[0]
            #foc = foc[gs[foc]!=9]
            if len(foc) == 0:
                continue
            BINS1[bbin] += [gi]*len(foc)
            BINS2[bbin] += list(gs[foc])
            WEIGHTS[bbin] += list(gli+gls[foc])

    # compute covariance
    COV, COVW, N, W = [], [], [], []
    for bbin in range(len(bins_left_bound)):
        bins1 = np.array(BINS1[bbin])
        bins2 = np.array(BINS2[bbin])
        weights = np.power(10, np.array(WEIGHTS[bbin]))
        if sum(bins1)==0 or sum(bins2)==0:
            COV += [np.nan]
            COVW += [np.nan]
            N += [0]
            W += [np.sum(weights)]
        else:
            COV += [np.cov(bins1, bins2, ddof=None, bias=0)[0][1]]
            COVW += [np.cov(bins1, bins2, ddof=None, bias=0, aweights=weights)[0][1]]
            N += [len(bins1)]
            W += [np.sum(weights)]
    R = np.column_stack(([chrom]*len(bins_left_bound),
                   bins_left_bound,
                   bins_left_bound + stepD_cM,
                   np.array(COV),
                   np.array(COVW),
                   np.array(N),
                   np.array(W)))
    
    if first:
        RR = R
        first = False
    else:
        RR = np.concatenate((RR, R), axis=0)
    
    with open(output_prefix+'.out', 'w') as fout:
        fout.write('chrom\tbin.left.bound\tbin.right.bound\tcov\tweighted.cov\tn.pairs\tsum.log10.lhoods\n')
        for i in range(RR.shape[0]):
            fout.write(str(int(RR[i,0]))+'\t'+'\t'.join([str(x) for x in RR[i,1:]])+'\n')
    fout.close()
    
print('It took:   '+str(round((time.time()-start)/60,2))+' min\n')
