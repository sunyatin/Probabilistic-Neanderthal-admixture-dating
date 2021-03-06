Method for estimating the age of admixture using a single diploid sample and genotype likelihoods.

## admixfrog2eigenstrat_v1.py
Converts a \.snp file output by AdmixFrog to 4 eigenstrat-formatted files with extensions .snp / .geno / .ind / .glhood. The .glhood file contains the likelihood of the genotypes.

**Parameters:**

    -p the AdmixFrog .snp file to convert

    -o prefix of the EIGENSTRAT output files

    -s add this option with the name of a file containing a list of SNPs that should be output (all SNPs present in the input .snp that are not listed in this -s file will be discarded), the -s file must have the same format as an EIGENSTRAT .snp file

## ASCovariance_SingleSample_v1c.py
Computes the weighted single sample statistic.

**Parameters:**

    -f  prefix of the eigenstrat files to analyze

    -p  name of the target population

    -o  output file prefix

    -minD minimum geneticdistance in centiMorgans

    -maxD maximum genetic distance in centiMorgans

    -stepD  bin size in centiMorgans

    --Morgans add this option if genetic distances are in Morgans (by default assumes centiMorgans)

    --chrom add this option with an integer to restrict the analysis to a specific chromosome
    
## global_decay_curve_admixture.R
Combines the per-chromosome decay curves into a single genome-wide decay curve and fits a single term exponential, estimates the error per parameters using a weighted jackknife procedure.
  
**Parameters (positional):**
  
    file containing the per-chromosome decay curves

    minimum genetic distance in cM

    maximum genetic distance in cM

    file containing the weight per chromosome for the jackknife (same format as for expfit_9.py)
  
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
 
