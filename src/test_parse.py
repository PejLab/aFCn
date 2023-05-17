from parse import * 
import pandas as pd
import pdb
import pytest

def test_read_eqtls():
    
    eqtl_test_path = "../tests/dataset/test_eqtls.tsv"
    eqtl_file = read_eqtls(eqtl_test_path)
    assert( len(eqtl_file.columns) == 19 )
    assert( eqtl_file.columns[0] == "gene_id" )

def test_read_eqtls_exception():

    eqtl_test_path = "../tests/dataset/test_wrongformat_eqtls.tsv"
    eqtl_file = pd.read_csv(eqtl_test_path, sep='\t')
   
    #check if function finds malformed variant IDs
    assert(is_variant_id_format_valid(eqtl_file) == False)
    
    
def test_read_expressions():
    
    eqtl_test_path = "../tests/dataset/test_eqtls.tsv"
    eqtl_file = read_eqtls(eqtl_test_path)

    class args:

        gct=False   

    expressions_test_path = '../tests/dataset/test_expressions.csv.gz'
    expressions_file = read_expressions(expressions_test_path, eqtl_file, args)
    
    assert( len(expressions_file.columns.tolist()) == 575 )
    
    
#def test_get_vcf_header():
    
#    vcf_test_path = '../tests/dataset/test.vcf.gz'
#    vcf_file = get_vcf_header(vcf_test_path)
    
#    assert( len(vcf_file) == 847 )

    
#def test_read_haplotypes():
    
#    eqtl_test_path = "../tests/dataset/test_eqtls.tsv"
#    eqtl_file = read_eqtls(eqtl_test_path)
    
#    test_read_haplotypes_path = '../tests/dataset/test.vcf.gz'
#    test_read_haplotypes_file = read_haplotypes(test_read_haplotypes_path, eqtl_file)
    
#    assert( len(test_read_haplotypes_file.columns.tolist()) == 847 )


#def test_read_malformatted_haplotypes():
    
#    eqtl_test_path = "../tests/dataset/test_eqtls.tsv"
#    eqtl_file = read_eqtls(eqtl_test_path)
    
#    test_read_haplotypes_path = '../tests/dataset/test_malformatted.vcf.gz'

    #malformed haplotypes file should raise exception
#    with pytest.raises(Exception):
#        test_read_haplotypes_file = read_haplotypes(test_read_haplotypes_path, eqtl_file)
    
def test_drop_noexpr_genes():
    
    eqtl_test_path = "../tests/dataset/test_eqtls.tsv"
    eqtl_file = read_eqtls(eqtl_test_path)
  
    class args:

        gct=False   

 
    expressions_test_path = '../tests/dataset/test_expressions.csv.gz'
    expressions_file = read_expressions(expressions_test_path, eqtl_file, args)
 
    eqtl_file = drop_noexpr_genes(eqtl_file, expressions_file ) 

    assert ( eqtl_file.columns.tolist()[0] == 'gene_id' ) 

#def test_drop_novcf_variants():
    
#    eqtl_test_path = "../tests/dataset/test_eqtls.tsv"
#    eqtl_file = read_eqtls(eqtl_test_path)
    
#    vcf_test_path = '../tests/dataset/test.vcf.gz'
#    vcf_file = read_haplotypes(vcf_test_path, eqtl_file)
 
#    eqtl_file = drop_novcf_variants(eqtl_file, vcf_file ) 

#    assert ( eqtl_file.columns.tolist()[0] == 'gene_id' ) 
 
#def test_read_covar():
    
#    pass

