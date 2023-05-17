from calc import *

def calc_gene(gene, inputs):

        '''Handles the in-thread aFC calculation'''

        [thisgene_expressions, thisgene_eqtls, thisgene_haplotypes, cov_dataf,args] = inputs

        #check if we have more than one individual for this gene
        if len(thisgene_expressions.columns.tolist()) < 2:

            return thisgene_eqtls.values.tolist()

        #Make numpy arrays out of the haplotype data - h1 is ref, h2 is alt
        [h1,h2] = make_haplotype_arrays(thisgene_haplotypes)

        #Drop values where there are nan values for haplotypes or expression
        [h1, h2, thisgene_expressions] = clean_matrix_nans(h1, h2, thisgene_expressions)

        #Was there a covariate matrix as an input?
        if args.cov != None:

            #Filter out unused invidiuals from our covariates file
            #needed_inds = list(set(h1.columns.tolist()) & set(cov_dataf.columns.tolist()))
            #cov_dataf = cov_dataf[needed_inds]
            #Do the linear fit
            cov_dataf = linear_estimate(h1,h2,thisgene_expressions,cov_dataf)

            #Do the least squares fit
            thisgene_eqtls = least_squares(h1, h2, thisgene_expressions, thisgene_eqtls, cov_dataf, args)

        else:

            #Do the linear fit
            [sa, c0] = nocovar_linear_estimate(h1,h2,thisgene_expressions)

            #Do the least squares fit
            thisgene_eqtls = nocovar_least_squares(h1, h2, thisgene_expressions, thisgene_eqtls, sa, c0, args)

        #return eqtl_matrix[eqtl_matrix.gene_id_clean == gene].values.tolist()
        return thisgene_eqtls.values.tolist()
