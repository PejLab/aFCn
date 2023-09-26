#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import gzip
import pdb
import pysam
import warnings
import copy

import multiprocessing
from joblib import Parallel, delayed
from tqdm import tqdm
warnings.filterwarnings("ignore")

from calc import *
from parse import * 
from thread_wrapper import *


def main():
    
    ###############################################################################
    #
    # Input args
    #
    ###############################################################################

    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True, help="Genotype VCF")
    parser.add_argument("--expr", required=True, help="Expressions file")
    parser.add_argument("--eqtl", required=True, help="File containing QTL to calculate allelic fold change for")
    parser.add_argument("--output", "-o", required=False, default="afcs_out.csv", help="Output file name")
    parser.add_argument("--nthreads", "-j", required=False, default="1", help="Number of threads to do fitting on")
    parser.add_argument("--conf", action="store_true", required=False, help="Calculate confidence intervals for aFC estimates")

    parser.add_argument("--normalize", action="store_true", required=False, help="This will log-transformm and normalize the expression file")
    parser.add_argument("--logtransform", action="store_true", required=False, help="This will log-transform the expression file")
    parser.add_argument("--redo", action="store_true", required=False, help="Filter out low quality eqtls based on confidence intervals")
    parser.add_argument("--splitexpr", action="store_true", required=False, help='If set, the individual names in the expressions file will be split on "-" characters and the parts of the name on the two side of the first "-" character will be retained.')
    parser.add_argument("--gct", action="store_true", required=False, help="Expressions matrix is in gct format (tab separated with a Description column)")

    args = parser.parse_args()
    args.cov = None 

    ###############################################################################
    #
    # Read in the expressions, haplotype data and eqtls 
    #
    ###############################################################################
    
    eqtl_filename = args.eqtl
    eqtl_dataframe = read_eqtls(eqtl_filename)
    #create a copy of this for storing the results later
    eqtl_dataframe_copy = copy.deepcopy(eqtl_dataframe)
    print("Done reading eqtls")    
 
    expressions_filename = args.expr
    expr_dataframe = read_expressions(expressions_filename, eqtl_dataframe, args)

    print("Done reading expressions")    
    #make the gene names the index
    expr_dataframe.index = expr_dataframe.Name
    #now drop it. cast the rest to a float
    expr_dataframe = expr_dataframe.drop(columns=["Name"]).astype(float) 
   
    vcfname = args.vcf
    vcf_dataframe = read_haplotypes(vcfname, eqtl_dataframe)
    print("Done reading vcf..")    
    
   

    ###############################################################################
    #
    # Do user requested input matrix transforms
    #
    ###############################################################################

    #Do we need to split the individual names by a - ?
    if args.splitexpr != False:
        expr_dataframe.columns = list(expr_dataframe.columns.str.split('-').str[:2].str.join('-'))           
    

    #Did the user request a log transform for the expressions matrix?
    if args.logtransform == True:

        print("Log transforming dataset")
        expr_dataframe = np.log2(expr_dataframe + 1)
 
    #Does the expressions matrix need normalization?
    if args.normalize == True:
        print("Normalizing dataset")
        #Lets pass the dataframe - this should contain the gene names as the
        #index and the rest of the dataset should be the expressions
        expr_dataframe = gene_expression_normalization(expr_dataframe)


    #Drop genes and variants that were not found in expressions or eqtl file
    eqtl_dataframe = drop_novcf_variants(eqtl_dataframe, vcf_dataframe)
    eqtl_dataframe = drop_noexpr_genes(eqtl_dataframe, expr_dataframe)
    

    ###############################################################################
    #
    # Filter out all the non-common individuals in the expressions and vcf matrix 
    #
    ###############################################################################
    
    #cast the individual list into required format inside vcf
    vcf_inds = list(vcf_dataframe.columns[9:])
    #split the individual ID and put it back together into required form if required

    expr_inds = expr_dataframe.columns.tolist()
    common_inds = list(set(vcf_inds).intersection(set(expr_inds)))
    
    
    ###############################################################################
    #
    # Use these individuals to create our beginning matrices
    #
    ###############################################################################

    #only use the common individuals in the expressions matrix    
    expression_matrix = expr_dataframe[common_inds]

    #lets replace NAs in the expression matrix with 0 expression
    expression_matrix = expression_matrix.fillna(0).astype(float)

    haplotype_matrix = vcf_dataframe[common_inds]
    #make the variant names the index
    haplotype_matrix.index = vcf_dataframe.ID

    #eqtl_matrix = eqtl_dataframe.drop(columns=['gene_id'])
    eqtl_matrix = eqtl_dataframe
    
    
    
    ###############################################################################
    #
    # Now do the fitting for each gene, then write the output matrix
    #
    ###############################################################################
    
    #set up the output matrix, based on the eqtl dataframe
    eqtl_matrix = init_output_df(eqtl_matrix, args)

    #this stays an empty string if no covariates were provided
    cov_dataf_full = ""

    if args.cov != None:

        covar_filename = args.cov
        cov_dataf_full = read_covar(covar_filename).astype(float)

    genes2process = list(eqtl_matrix.gene_id_clean.unique())
    genecounter = 0
    #total number of genes to process using threads
    totalgenes = len(genes2process)

    #this list contains our processed eqtls
    processed_list = []

    #for this gene, get expressions, etc
    inputs = {}
    print("Preparing data for parallelization..")
    for gene in genes2process:

        thisgene_expressions = expression_matrix[expression_matrix.index.isin([gene])]
        #variants in this gene
        thisgene_variants = eqtl_matrix[eqtl_matrix.gene_id_clean == gene].variant_id_clean
        #haplotypes for these variants for all inds
        thisgene_haplotypes = haplotype_matrix[haplotype_matrix.index.isin(list(thisgene_variants))]
        #the eqtl entries for this gene
        thisgene_eqtls = eqtl_matrix[(eqtl_matrix.gene_id_clean == gene) & (eqtl_matrix.variant_id_clean.isin(list(thisgene_variants)))]    
        #save in dictionary
        inputs[gene] = [thisgene_expressions, thisgene_eqtls, thisgene_haplotypes, cov_dataf_full, args]

    print("Processing genes on " + str(args.nthreads) + " thread(s)" )

    ingenes = tqdm(genes2process)
    processed_list = Parallel(n_jobs=int(args.nthreads))(delayed(calc_gene)(gene,inputs[gene]) for gene in ingenes)

    #expand the results, since there are several eqtls per gene
    expanded_results = [item for sublist in processed_list for item in sublist]
    #combine processed dataframes
    final_df = pd.DataFrame(expanded_results, columns=eqtl_matrix.columns)

    #check that we have all eqtls in the dataframe
    original_columns = eqtl_dataframe_copy.columns.tolist()
    #separate the gene_id column to be reattached after the merge 
    out_df = pd.merge(eqtl_dataframe_copy, final_df, how='left',on=original_columns )

    #now drop the gene_id_clean column
    out_df = out_df.drop(columns=["gene_id_clean", "variant_id_clean"])
    
    ###############################################################################
    #
    # Now do some post-processing checks and save results
    #
    ###############################################################################
    
    #if more than 25% of the aFCs are missing, print a warning
    if len(out_df.log2_aFC.dropna())/len(out_df.log2_aFC) < 0.75:

        print("WARNING: aFCs could not be calculated for more than 25% of EQTLs")

    #Save results
    print("Writing results to file: " + str(args.output)) 
    out_df.to_csv(args.output, index=None)
    

     


if __name__ == "__main__":
    
    main()
    
    
