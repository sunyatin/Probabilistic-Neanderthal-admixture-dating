#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 12:42:34 2019

@author: rtournebize
"""

import numpy as np
import argparse
import sys
import random

parser = argparse.ArgumentParser(description='')
parser.add_argument('-f', '--input', type=str, required=True, help='File to convert.')
parser.add_argument('-o', '--outfilePrefix', type=str, required=True, help='Prefix of the output file.')
parser.add_argument('-s', '--snpFile', type=str, required=False, default=None, help='List of SNPs to retain.')

args = parser.parse_args()

input_prefix = args.input
output_prefix = args.outfilePrefix
snpFile = args.snpFile

if snpFile is not None:
    targetSNP = np.genfromtxt(snpFile, dtype=str, usecols = (1,3))
    targetSNP = set([targetSNP[x,0]+':'+targetSNP[x,1] for x in range(targetSNP.shape[0])])
    
    SNP = np.genfromtxt(input_prefix, dtype=str, usecols = (0,1,2), delimiter=',', skip_header = 1)
    
    print('Finding intersection')
    
    goodSNPs = [False]*SNP.shape[0]
    for i in range(SNP.shape[0]):
        snp = SNP[i,0]+':'+SNP[i,1]
        if snp in targetSNP:
            goodSNPs[i] = True
            
    print('Exporting')
    
    with open(output_prefix+'.snp', 'w') as fout:
        for i in range(SNP.shape[0]):
            if goodSNPs[i] is True:
                fout.write('. '+str(SNP[i,0])+' '+SNP[i,2]+' '+SNP[i,1]+' N N\n')
    with open(output_prefix+'_Morgans.snp', 'w') as fout:
        for i in range(SNP.shape[0]):
            if goodSNPs[i] is True:   
                fout.write('. '+str(SNP[i,0])+' '+str(float(SNP[i,2])*0.01)+' '+SNP[i,1]+' N N\n')
    
    print('Initial number of SNPs: '+str(SNP.shape[0]))
    print('Number of output SNPs: '+str(sum(goodSNPs)))
    
    del SNP
    
    GENO = np.genfromtxt(input_prefix, dtype=float, usecols = (5,6,7), delimiter=',', skip_header = 1)
    with open(output_prefix+'.geno', 'w') as fout:
        with open(output_prefix+'.glhood', 'w') as fout2:
            for i in range(GENO.shape[0]):
                if goodSNPs[i] is True:
                    gl = np.amax(GENO[i,:])
                    got = np.where(GENO[i,:] == gl)[0]
                    if len(got) != 1:
                        print('Two equally probable genotype, picking one random')
                        gt = int(random.choice(got))
                    else:
                        gt = int(got)
                    fout.write(str(gt)+'\n')
                    fout2.write(str(gl)+'\n')
    
    with open(output_prefix+'.ind', 'w') as fout:
        fout.write(output_prefix+' . '+output_prefix+'\n')

    
else:

    SNP = np.genfromtxt(input_prefix, dtype=str, usecols = (0,1,2), delimiter=',', skip_header = 1)
    
    print('Exporting')
    
    with open(output_prefix+'.snp', 'w') as fout:
        for i in range(SNP.shape[0]):
            fout.write('. '+str(SNP[i,0])+' '+SNP[i,2]+' '+SNP[i,1]+' N N\n')
    with open(output_prefix+'_Morgans.snp', 'w') as fout:
        for i in range(SNP.shape[0]):
            fout.write('. '+str(SNP[i,0])+' '+str(float(SNP[i,2])*0.01)+' '+SNP[i,1]+' N N\n')
    del SNP
    
    GENO = np.genfromtxt(input_prefix, dtype=float, usecols = (5,6,7), delimiter=',', skip_header = 1)
    with open(output_prefix+'.geno', 'w') as fout:
        with open(output_prefix+'.glhood', 'w') as fout2:
            for i in range(GENO.shape[0]):
                gl = np.amax(GENO[i,:])
                gt = int(np.where(GENO[i,:] == gl)[0])
                fout.write(str(gt)+'\n')
                fout2.write(str(gl)+'\n')
    
    with open(output_prefix+'.ind', 'w') as fout:
        fout.write(output_prefix+' . '+output_prefix+'\n')

