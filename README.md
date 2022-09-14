Method for estimating the age of admixture using a single diploid sample and genotype likelihoods.

## admixfrog2eigenstrat_v1.py
Converts a \.snp file output by AdmixFrog to 4 eigenstrat-formatted files with extensions .snp / .geno / .ind / .glhood. The .glhood file contains the likelihood of the genotypes.

**Parameters:**

    -p the AdmixFrog .snp file to convert

    -o prefix of the EIGENSTRAT output files

    -s add this option with the name of a file containing a list of SNPs that should be output (all SNPs present in the input .snp that are not listed in this -s file will be discarded), the -s file must have the same format as an EIGENSTRAT .snp file

## ASCovariance_SingleSample_v1c.py
Computes the weighted single sample statistic. Note that by default the script assumes that the **genetic positions provided in the \*.snp file are in centiMorgans**.

**Parameters:**

    -f  prefix of the eigenstrat files to analyze

    -p  name of the target population

    -o  output file prefix

    -minD minimum geneticdistance in centiMorgans

    -maxD maximum genetic distance in centiMorgans

    -stepD  bin size in centiMorgans

    --Morgans add this option if genetic distances are in Morgans (by default assumes centiMorgans)

    --chrom add this option with an integer to restrict the analysis to a specific chromosome
    
Note that if you add the switch --Morgans, the output will still report distances in cM as this is the default unit used by all the scripts provided here. Therefore the switch only refers to the unit of the input \*.snp file.
    
 **Output:**
 
 A single file with suffix *.out*, containing the decay curves calculated for each chromosome and each distance bin, with 7 tab-delimited columns:
 
     chrom  the chromosome ID
     bin.left.bound the left boundary of the distance bin (in cM)
     bin.right.bound    the right boundary of the distance bin (in cM)
     cov    the value of the covariance using the called genotype reported in the *.geno file
     weighted.cov   the value of the covariance weighted by the genotype likelihood reported in the *.glhood file
     n.pairs    the number of SNP pairs used to calculate the covariance in the distance bin
     sum.log10.lhoods   the sum of the log10 of the genotype likelihoods of all SNPs used in the distance bin
    
## global_decay_curve_admixture.R
Combines the per-chromosome decay curves into a single genome-wide decay curve and fits a single term exponential, estimates the error per parameters using a weighted jackknife procedure.
  
**Parameters (positional):**
  
    file containing the per-chromosome decay curves (i.e. the \*.out file output by ASCovariance_SingleSample_v1c.py)

    minimum genetic distance in cM

    maximum genetic distance in cM

    file containing the weight per chromosome for the jackknife (same format as for expfit_9.py)
    
/!\ IMPORTANT. The genetic distances in the \*.out file **must be in centiMorgans** as there is a systematic conversion in the R script afterwards.
    
**Output:**

Three files using the called genotype (\*.hard.\*) and three files using the genotype likelihoods (\*.soft.\*):

**\*.hard** or **\*.soft** This file gives the decay curve (using the called genotypes or the genotype likelihoods, respectively) averaged over all chromosomes; it has 2 columns: the distance bin (in cM) and the covariance value at the bin

**\*.jin** Provides the rate of the exponential decay (i.e. the admixture age, in units of numbers of generations) for each run of the block weighted jackknife (e.g. first line is the estimate when removing the first chromosome)

**\*.fit** Provides the rate of the exponential decay (i.e. the admixture age, in units of nubmers of generations) as calculated from the decay curve averaged over all chromosome (2nd line); as the jackknife-average (3rd line); and its standard error from the weighted jackknife procedure (4th line). We advise reported the jackknife-average (instead of the global average from the line before it).
  
## expfit_v9.py
Performs the fitting of a single exponential or a sum of two exponentials.
  
**Parameters:**
  
    -f  name of the file containing a the decay curve per chromosome
    
    -p  name of the target population
    
    -o  output file prefix
    
    -n  block sizes
    
    -minD minimum geneticdistance in centiMorgans
    
    -maxD maximum genetic distance in centiMorgans
    
    --noBgLDSubstraction  add this option if you do not want to subtract the within-population statistic by a cross-population statistic
    
    --exp2  add this option to fit a sum of two exponentials (by default will fit a single exponential with affine term)
 
## LRT_1exp_2exp_lomax.py
A very generic script to fit a single exponential, a sum of two exponentials and a Lomax distribution to a (X,Y) distribution. Also performs likelihood ratio tests between these 3 models.
  
**Parameter (positional):**
  
    A file containing two columns: 
    
    (1) genetic distance in cM; 
    
    (2) a statistic (e.g. weighted covariance, allele sharing correlation, etc.)
 
