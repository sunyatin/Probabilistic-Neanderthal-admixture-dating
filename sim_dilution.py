#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: rtournebize
"""
import msprime
import numpy as np
import argparse
import sys

# 12092019: implements two demographic models

parser = argparse.ArgumentParser(description='')
parser.add_argument('-T', '--time', type=int, required=True, help='Admixture time between EU and BE')
parser.add_argument('-p', '--proportion', type=float, required=True, help='Admixture proportion between EU and BE')
parser.add_argument('-S', '--sampling', type=int, required=True, help='Sampling time of EU ancient samples')
parser.add_argument('-o', '--outfilePrefix', type=str, required=True, help='Output file prefix')
parser.add_argument('-m', '--model', type=str, required=False, default='qiaomei', help='Name of the demographic model')
args = parser.parse_args()

Ta_BE_EU = args.time
Ts = args.sampling
Aa_BE_EU = args.proportion
prefix = args.outfilePrefix
model = args.model

print('Using model: '+model)

######### Parameters

L = 50e6
recomb_rate = 1e-8
mut_rate = 1.2e-8
n_AFR = 10 # diploids
n_EU = 1 # diploids
n_NEA = 1 # diploids

nrep = 10

######### Script

SNP = open(prefix+'.snp', 'w')
GENO = open(prefix+'.geno', 'w')
GLHOOD = open(prefix+'.glhood', 'w')

SNP1 = open(prefix+'.1.snp', 'w')
GENO1 = open(prefix+'.1.geno', 'w')
GLHOOD1 = open(prefix+'.1.glhood', 'w')

samples = [msprime.Sample(population = 0, time = 0)]*n_AFR*2 # Africans
samples.extend([msprime.Sample(population = 1, time = Ts)]*n_EU*2) # Ancient European
#samples.extend([msprime.Sample(population = 2, time = 0)]) # Basal Eurasians
#samples.extend([msprime.Sample(population = 3, time = 0)]) # Neanderthal
samples.extend([msprime.Sample(population = 4, time = 2400)]*n_NEA*2) # Altai

print(samples)

migr_matrix = [[0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0]]

if model == 'qiaomei':
    
    mut_rate_ori = 1.5e-8
    ratio = mut_rate_ori / mut_rate
    
    sys.exit('the demo has not been rescaled')
    
    if False:

        pop_conf = [
        msprime.PopulationConfiguration(initial_size = 14000), # 0 = Africans
        msprime.PopulationConfiguration(initial_size = 33800), # 1 = Europeans
        msprime.PopulationConfiguration(initial_size = 12500), # 2 = Basal Eurasians
        msprime.PopulationConfiguration(initial_size = 2500), #  3 = Neanderthal
        msprime.PopulationConfiguration(initial_size = 2500)] #  4 = Altai
        
        g_E = 1/2000 * np.log(33800 / 1032)
        
        demo = [
                msprime.PopulationParametersChange(time = 0, growth_rate = g_E, population_id = 1), # EU
                msprime.MassMigration(time = Ta_BE_EU, source = 1, destination = 2, proportion = Aa_BE_EU), # EU-BE admix
                msprime.PopulationParametersChange(time = 2000, initial_size = 14000, growth_rate = 0, population_id = 1), # EU NOT SURE!!!!
                msprime.PopulationParametersChange(time = 2200, initial_size = 1860, population_id = 1), # EU
                msprime.MassMigration(time = 2201, source = 1, destination = 3, proportion = 0.03), # EU-NEA admix
                msprime.MassMigration(time = 3000, source = 1, destination = 0, proportion = 1), # EU-AFR split
                msprime.MassMigration(time = 3001, source = 2, destination = 0, proportion = 1), # BE-AFR split
                msprime.MassMigration(time = 4000, source = 4, destination = 3, proportion = 1), # Altai-Neanderthal split
                msprime.PopulationParametersChange(time = 6000, initial_size = 7300, population_id = 0), # AFR
                msprime.PopulationParametersChange(time = 12000, initial_size = 9800, population_id = 0), # AFR
                msprime.MassMigration(time = 12001, source = 3, destination = 0, proportion = 1)] # NEA-AFR split
    
elif model == 'harris':
    # Harris and Nielsen 2006
    
    mut_rate_ori = 2.5e-8
    ratio = mut_rate_ori / mut_rate
    
    pop_conf = [
    msprime.PopulationConfiguration(initial_size = 10000*ratio), # 0 = Africans
    msprime.PopulationConfiguration(initial_size = 20000*ratio), # 1 = Europeans
    msprime.PopulationConfiguration(initial_size = 10000*ratio), # 2 = Basal Eurasians
    msprime.PopulationConfiguration(initial_size = 1000*ratio), #  3 = Neanderthal
    msprime.PopulationConfiguration(initial_size = 1000*ratio)] #  4 = Altai
    
    g_E = 1/(1100*ratio) * np.log(20000 / 1032)
    
    demo = [
            msprime.PopulationParametersChange(time = 0, growth_rate = g_E, population_id = 1), # EU
            msprime.MassMigration(time = Ta_BE_EU, source = 1, destination = 2, proportion = Aa_BE_EU), # EU-BE admix
            msprime.PopulationParametersChange(time = 1100*ratio, initial_size = 1032*ratio, growth_rate = 0, population_id = 1), # EU
            msprime.PopulationParametersChange(time = 2000*ratio, initial_size = 10000*ratio, population_id = 1), # EU
            msprime.MassMigration(time = 2000*ratio, source = 1, destination = 3, proportion = 0.03), # EU-NEA admix
            msprime.MassMigration(time = 3000*ratio, source = 1, destination = 0, proportion = 1), # EU-AFR split
            msprime.MassMigration(time = 3000*ratio, source = 2, destination = 0, proportion = 1), # BE-AFR split
            msprime.MassMigration(time = 4000*ratio, source = 4, destination = 3, proportion = 1), # Altai-Neanderthal split
            msprime.MassMigration(time = 18000*ratio, source = 3, destination = 0, proportion = 1)] # NEA-AFR split
    
elif model == 'harris_structure':
    # Harris and Nielsen 2006
    
    mut_rate_ori = 2.5e-8
    ratio = mut_rate_ori / mut_rate
    
    pop_conf = [
    msprime.PopulationConfiguration(initial_size = 10000*ratio), # 0 = Africans
    msprime.PopulationConfiguration(initial_size = 20000*ratio), # 1 = Europeans
    msprime.PopulationConfiguration(initial_size = 20000*ratio), # 2 = Europeans 2
    msprime.PopulationConfiguration(initial_size = 1000*ratio), #  3 = Neanderthal
    msprime.PopulationConfiguration(initial_size = 1000*ratio)] #  4 = Altai
    
    g_E = 1/(1100*ratio) * np.log(20000 / 1032)
    
    demo = [
            msprime.PopulationParametersChange(time = 0, growth_rate = g_E, population_id = 1), # EU
            msprime.PopulationParametersChange(time = 0, growth_rate = g_E, population_id = 2), # EU2
            msprime.MassMigration(time = Ta_BE_EU, source = 1, destination = 2, proportion = Aa_BE_EU), # EU-EU2 admix
            msprime.PopulationParametersChange(time = 1100*ratio, initial_size = 1032*ratio, growth_rate = 0, population_id = 1), # EU
            msprime.PopulationParametersChange(time = 2000*ratio, initial_size = 10000*ratio, population_id = 1), # EU
            msprime.PopulationParametersChange(time = 1100*ratio, initial_size = 1032*ratio, growth_rate = 0, population_id = 2), # EU2
            msprime.PopulationParametersChange(time = 2000*ratio, initial_size = 10000*ratio, population_id = 2), # EU2
            msprime.MassMigration(time = 2000*ratio, source = 1, destination = 3, proportion = 0.03), # EU-NEA admix
            msprime.MassMigration(time = 2000*ratio, source = 2, destination = 3, proportion = 0.032), # EU-NEA admix
            msprime.MassMigration(time = 3000*ratio, source = 1, destination = 0, proportion = 1), # EU-AFR split
            msprime.MassMigration(time = 3000*ratio, source = 2, destination = 0, proportion = 1), # EU2-AFR split
            msprime.MassMigration(time = 4000*ratio, source = 4, destination = 3, proportion = 1), # Altai-Neanderthal split
            msprime.MassMigration(time = 18000*ratio, source = 3, destination = 0, proportion = 1)] # NEA-AFR split
   
elif model == 'harris_no_expansion':
    # Harris and Nielsen 2006
    
    mut_rate_ori = 2.5e-8
    ratio = mut_rate_ori / mut_rate
    
    pop_conf = [
    msprime.PopulationConfiguration(initial_size = 10000*ratio), # 0 = Africans
    msprime.PopulationConfiguration(initial_size = 20000*ratio), # 1 = Europeans
    msprime.PopulationConfiguration(initial_size = 20000*ratio), # 2 = Europeans 2
    msprime.PopulationConfiguration(initial_size = 1000*ratio), #  3 = Neanderthal
    msprime.PopulationConfiguration(initial_size = 1000*ratio)] #  4 = Altai
    
    demo = [
            msprime.PopulationParametersChange(time = 0, growth_rate = 0, population_id = 1), # EU
            msprime.PopulationParametersChange(time = 2000*ratio, initial_size = 10000*ratio, population_id = 1), # EU
            msprime.MassMigration(time = 2000*ratio, source = 1, destination = 3, proportion = 0.03), # EU-NEA admix
            msprime.MassMigration(time = 3000*ratio, source = 1, destination = 0, proportion = 1), # EU-AFR split
            msprime.MassMigration(time = 3000*ratio, source = 2, destination = 0, proportion = 1), # BE-AFR split
            msprime.MassMigration(time = 4000*ratio, source = 4, destination = 3, proportion = 1), # Altai-Neanderthal split
            msprime.MassMigration(time = 18000*ratio, source = 3, destination = 0, proportion = 1)] # NEA-AFR split
    

else:
    sys.exit('no model provided')

demo = sorted(demo, key = lambda x: x.time)

dp = msprime.DemographyDebugger(population_configurations = pop_conf,
                            migration_matrix = migr_matrix,
                            demographic_events = demo)
dp.print_history()

IND = open(prefix+'.ind', 'w')
n = 1
for i in range(n_AFR):
    IND.write(str(n)+' U AFR\n')
    n += 1
for i in range(n_EU):
    IND.write(str(n)+' U EU\n')
    n += 1
for i in range(n_NEA):
    IND.write(str(n)+' U NEA\n')
    n += 1
IND.close()

IND = open(prefix+'.1.ind', 'w')
n = 1
for i in range(n_EU):
    IND.write(str(n)+' U EU\n')
    n += 1
IND.close()


print('NOTE.\nExporting only ascertained SNP, i.e. fixed ancestral in AFR and at least one derived in Altai.\n')

for chrom in range(nrep):
    print(chrom+1)
    
    tree = msprime.simulate(samples = samples,
                            population_configurations = pop_conf,
                            migration_matrix = migr_matrix,
                            length = L,
                            recombination_rate = recomb_rate,
                            mutation_rate = mut_rate,
                            demographic_events = demo)
    
    G = tree.genotype_matrix()
    G = np.asmatrix(G)
    
    POS = []
    for site in tree.sites():
        POS.append(site.position)
    D = [x*recomb_rate*100 for x in POS] # in cM
    
    for i in range(len(POS)):
        g = G[i,:].A1
        g_AFR = g[0:n_AFR] + g[n_AFR:(n_AFR*2)]
        g_EU = g[(n_AFR*2):(n_AFR*2+n_EU)] + g[(n_AFR*2+n_EU):(n_AFR*2+n_EU*2)]
        g_NEA = g[(n_AFR*2+n_EU*2):(n_AFR*2+n_EU*2+n_NEA)] + g[(n_AFR*2+n_EU*2+n_NEA):(n_AFR*2+n_EU*2+n_NEA*2)]
        
        # ascertainment check
        asc = (sum(g_AFR)==0) and (sum(g_NEA)>=1)
        
        if asc == True:
            GENO.write(''.join([str(x) for x in g_AFR]) + ''.join([str(x) for x in g_EU]) + ''.join([str(x) for x in g_NEA])+'\n')
            SNP.write('. '+str(chrom+1)+' '+str(D[i])+' '+str(int(POS[i]))+' . .\n')
            GLHOOD.write(''.join([str(0) for x in g_AFR]) + ''.join([str(0) for x in g_EU]) + ''.join([str(0) for x in g_NEA])+'\n')
            GENO1.write(''.join([str(x) for x in g_EU])+'\n')
            SNP1.write('. '+str(chrom+1)+' '+str(D[i])+' '+str(int(POS[i]))+' . .\n')
            GLHOOD1.write(''.join([str(0) for x in g_EU])+'\n')
    
SNP.close()
GENO.close()
GLHOOD.close()

SNP1.close()
GENO1.close()
GLHOOD1.close()

