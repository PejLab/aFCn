import pandas as pd
import numpy as np
import argparse
import gzip
import pdb
import pysam
import warnings
import sys
warnings.filterwarnings("ignore")

def is_variant_id_format_valid(eqtl_file):

    '''Checks if the variant ids in the eqtl matrix begin with a valid character'''

    #get the first character of all rows
    variant_ids = eqtl_file['variant_id'].tolist()
    firstchar_isdigit = [str(var_id)[0].isalpha() for var_id in variant_ids]

    #check if any of these characters is a digit
    if True in firstchar_isdigit:

        return True

    else:

        return False


def read_eqtls(eqtl_filename):
    
    '''Used to read in the eqtl file into a pandas dataframe'''

    #open the eqtl file
    eqtl_file = pd.read_csv(eqtl_filename, sep='\t')
    
    #drop gene versions
    eqtl_file['gene_id_clean'] = eqtl_file.gene_id.str.split('.').str[0]
    
    if is_variant_id_format_valid(eqtl_file) == False:

        raise Exception('''Variant IDs must not begin with a digit, 
                reformat your vcf and EQTL matrix''')
 

    return eqtl_file


def read_expressions(expressions_filename, eqtl_file, args):

    '''Used to read the relevant expressions into a pandas dataframe'''
    
    chosen_genes = list(eqtl_file.gene_id_clean)
    expr_df = []
    expr_df_cols = []
    #open the expression file, get the expressions for the genes we have eqtls for
    with gzip.open(expressions_filename) as f:

        for line in f:

            line = line.decode()
            line = line.replace('\n','')

            #if the data is in gct format, use tabs
            if args.gct != False: 
                expr_data = line.split('\t')
            else:
                expr_data = line.split(',')

            if expr_data[0] == 'Name':

                expr_df_cols = expr_data


            else:

                #try splitting the gene by version number, in 
                #case our user forgot
                thisgene = expr_data[0].split(".")[0]
                
                if thisgene in chosen_genes:
                    expr_df.append(expr_data)

    expr_dataframe = pd.DataFrame(expr_df, columns=expr_df_cols)

    #if no genes were found
    if expr_dataframe.empty:
        raise Exception("No matching genes found between the EQTL and expressions files")

    #if data is in gct format, drop description column
    if args.gct != False: 

        if "Description" not in expr_df_cols:
            raise Exception("Description column not found in expressions matrix")

        expr_dataframe = expr_dataframe.drop(columns=['Description'])

 
    #now do the gene version correction once and for all if needed
    expr_dataframe['Name'] = expr_dataframe['Name'].str.split(".").str[0]
    
    return expr_dataframe



def get_vcf_header(vcfName):

    ''' This function returns the header of a gzipped or non-gzipped vcf file '''

    #open the gzipped file 
    vcf_stream = gzip.open(vcfName)


    #find the first line that does not begin with "#"
    for line in vcf_stream:

        if isinstance(line, bytes) and not isinstance(line, str):

            line = line.decode()

        if '#CHROM' in line:

            return line.replace('\n','').split('\t')
        
        
        
def read_haplotypes(vcfName, eqtl_file):

    '''Reads in the variants we have eqtls and expressions(genes) for'''
    
    chosen_variants = eqtl_file.variant_id.unique()

    #match these to the vcf file
    vcf_df = []
    vcf_df_cols = []
    #the eqtl variants
    counter = 0
    tabix_haplotypes = pysam.Tabixfile( vcfName,"r")


    #write the header first
    vcf_df_cols = get_vcf_header(vcfName)

    for variant in chosen_variants:

        #get coordinates
        [chrom, pos] = variant.split("_")[:2]

        #query vcf using tabix
        try:

            records = tabix_haplotypes.fetch(chrom, int(pos) - 1, int(pos))

            for record in records:

                vcf_df.append(record.split('\t'))


        except:

            #variant not found
            continue

    vcf_dataframe = pd.DataFrame(vcf_df, columns=vcf_df_cols).drop_duplicates()

    #do a quick check and raise an exception if no eqtls were found
    if vcf_dataframe.empty:

        raise Exception("No mathing EQTLs were found in VCF file - check variant ID formatting")

    
    return vcf_dataframe


def drop_noexpr_genes(eqtl_dataframe, expr_dataframe):
    
    '''Drop variants that belong to genes we have no expressions for from the eqtl matrix'''

    corrected_eqtl_dataframe = eqtl_dataframe[eqtl_dataframe.gene_id_clean.isin(expr_dataframe.index.tolist())]
    
    return corrected_eqtl_dataframe


def drop_novcf_variants(eqtl_dataframe, vcf_dataframe):
    
    '''Drop variants that we couldnt get from the vcf file from the eqtl matrix'''

    corrected_eqtl_dataframe = eqtl_dataframe[eqtl_dataframe.variant_id.isin(vcf_dataframe.ID)]
    
    return corrected_eqtl_dataframe
    


        
  
