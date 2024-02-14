import pandas as pd
import numpy as np
import argparse
import gzip
import pdb
import pysam
import warnings
import sys
warnings.filterwarnings("ignore")


def read_eqtls(eqtl_filename):
    
    '''Used to read in the eqtl file into a pandas dataframe'''

    #open the eqtl file
    eqtl_file = pd.read_csv(eqtl_filename, sep='\t')
    
    # #drop gene versions
    # eqtl_file['gene_id_clean'] = eqtl_file.gene_id.str.split('.').str[0]
    
    # Clean up gene and variant IDs
    gene_map = {gene_id: f'g{i}' for i, gene_id in enumerate(eqtl_file.gene_id.unique())}
    eqtl_file['gene_id_clean'] = eqtl_file.gene_id.map(gene_map)

    variant_map = {variant_id: f'v{i}' for i, variant_id in enumerate(eqtl_file.variant_id.unique())}
    eqtl_file['variant_id_clean'] = eqtl_file.variant_id.map(variant_map)

    # if is_variant_id_format_valid(eqtl_file) == False:

    #     raise Exception('''Variant IDs must not begin with a digit, 
    #             reformat your vcf and EQTL matrix''')
 

    return eqtl_file


def read_expressions(expressions_filename, eqtl_file, args):

    '''Used to read the relevant expressions into a pandas dataframe'''
    
    gene_map = eqtl_file.loc[:, ['gene_id', 'gene_id_clean']].set_index('gene_id')['gene_id_clean'].to_dict()
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
                thisgene = expr_data[0]
                
                if thisgene in gene_map:
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
    expr_dataframe['Name'] = expr_dataframe['Name'].map(gene_map)
    
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
    
    variant_map = eqtl_file.loc[:, ['variant_id', 'variant_id_clean']].set_index('variant_id')['variant_id_clean'].to_dict()

    #match these to the vcf file
    vcf_df = []
    vcf_df_cols = []
    #the eqtl variants
    counter = 0
    tabix_haplotypes = pysam.Tabixfile( vcfName,"r")


    #write the header first
    vcf_df_cols = get_vcf_header(vcfName)

    for variant in variant_map:

        #get coordinates
        chrom = eqtl_file.loc[eqtl_file['variant_id']==variant, ['variant_chr']].iloc[0,0]
        pos = eqtl_file.loc[eqtl_file['variant_id']==variant, ['variant_pos']].iloc[0,0]

        #query vcf using tabix
        try:

            records = tabix_haplotypes.fetch(chrom, int(pos) - 1, int(pos))

            for record in records:
                record_info = record.split('\t')
                if record_info[2] == variant:
                    vcf_df.append(record_info)


        except:

            #variant not found
            continue

    vcf_dataframe = pd.DataFrame(vcf_df, columns=vcf_df_cols).drop_duplicates()

    # Use clean variant ids
    vcf_dataframe['ID'] = vcf_dataframe.ID.map(variant_map)

    #do a quick check and raise an exception if no eqtls were found
    if vcf_dataframe.empty:

        raise Exception("No matching EQTLs were found in VCF file - check variant ID formatting")

    
    return vcf_dataframe


def drop_noexpr_genes(eqtl_dataframe, expr_dataframe):
    
    '''Drop variants that belong to genes we have no expressions for from the eqtl matrix'''

    corrected_eqtl_dataframe = eqtl_dataframe[eqtl_dataframe.gene_id_clean.isin(expr_dataframe.index.tolist())]
    
    return corrected_eqtl_dataframe


def drop_novcf_variants(eqtl_dataframe, vcf_dataframe):
    
    '''Drop variants that we couldnt get from the vcf file from the eqtl matrix'''

    corrected_eqtl_dataframe = eqtl_dataframe[eqtl_dataframe.variant_id_clean.isin(vcf_dataframe.ID)]
    
    return corrected_eqtl_dataframe
    


        
  
