import pandas as pd
import numpy as np
import argparse
import gzip
import pdb
import pysam
import warnings
warnings.filterwarnings("ignore")


from lmfit import Parameters
from lmfit import minimize
from sklearn.linear_model import LinearRegression
import lmfit

def fcn2fit(params,h1,h2,variants,expr,covarlist ):
    
        """
        Do least squares optimization with a covariate matrix
        """
    
        #the covariate list
        cov_dataf = covarlist
        used_inds = list(set(covarlist.columns) & set(h1.columns) & set(expr.columns))

        #id_cov is just a list of covariate factors
        id_cov = list(set(covarlist.index).symmetric_difference(set(variants)))

        expr=expr[used_inds].values

        #make sure that we have the correct number of inds

        h1 = h1[used_inds].values
        h2 = h2[used_inds].values

        log2_c0=params['C0'].value
        sa=[]
        wa = []
        for arrayvalue in variants:

            sa.append(params[str(arrayvalue)].value)
           
        #get covariates only from matrix
        k= cov_dataf.transpose()[id_cov].transpose()

        #now make an array of params for each id
        for arrayvalue in id_cov:

            wa.append(params[str(arrayvalue)].value)

        #get expressions for only the common samples
        #calculate model
        dotprod1 = np.dot(np.array(sa) , np.array(h1))
        dotprod2 = np.dot(np.array(sa) , np.array(h2))

        model = np.log2(2**(dotprod1) + 2**(dotprod2)) + log2_c0 + np.dot(k[k.columns[:-1]].transpose().values , wa) #also add the covar matrix

        #expressions were log transformed outside function
        error = (model - expr)

        #return error
        return error
    
    
def nocovar_fcn2fit(params,h1,h2,variants,expr):
    
        """
        Do least squares optimization without a covariate matrix
        """
        
        log2_c0=params['C0'].value
        sa=[]

        for arrayvalue in variants:

            sa.append(params[str(arrayvalue)].value)

        #get expressions for only the common samples
        expr=expr[h1.columns].values
        
        #lets filter h2 and h2 for the variants we have
        h1 = h1.loc[variants].values
        h2 = h2.loc[variants].values
        
        #calculate model
        dotprod1 = np.dot(np.array(sa).reshape((1, len(sa))),h1)
        dotprod2 = np.dot(np.array(sa).reshape((1, len(sa))),h2)

        model = np.log2(2**(dotprod1) + 2**(dotprod2)) + log2_c0 
        #take the log2 of the expr (np.log2(expr)) if the input expressions were not log transformed!!!
        error = model - expr 

        return error
    
    
def gene_expression_normalization(counts):
        
        """
        If the option is passed, normalize gene expression matrix
        """
        original_counts = counts.copy()
        original_counts = original_counts + 1
        counts = counts + 1
        counts = counts[np.alltrue(counts > 100, axis=1)]
        logcounts = np.log2(counts)
        log_original_counts = np.log2(original_counts)
        # median of rows
        loggeommeans = np.median(logcounts, axis=1).reshape(len(logcounts), 1)
        loggeommeans_expand = np.repeat(loggeommeans, logcounts.shape[1] , axis=1)
        # median of columns
        sf = np.median(logcounts - loggeommeans_expand, axis=0).reshape(1,logcounts.shape[1]);
        sf_expand = np.repeat(sf,log_original_counts.shape[0], axis = 0);
        normalized_matrix = log_original_counts - sf_expand;
        return (normalized_matrix);




def init_output_df(eqtl_dataf, args):
    
    '''Adds new output columns to the original eqtl file to save outputs in'''

    #add a new column to eqtl_dataf to store afc values       
    eqtl_dataf['log2_aFC'] = np.full(len(eqtl_dataf.index.values),np.nan) 
    eqtl_dataf['log2_aFC_error'] = np.full(len(eqtl_dataf.index.values),np.nan) 
    eqtl_dataf['log2_aFC_c0'] = np.full(len(eqtl_dataf.index.values),np.nan) 

    if args.conf == True:

        #confidence intervals newly added
        eqtl_dataf['log2_aFC_min_95_interv'] = np.full(len(eqtl_dataf.index.values),np.nan) 
        eqtl_dataf['log2_aFC_plus_95_interv'] = np.full(len(eqtl_dataf.index.values),np.nan) 

        eqtl_dataf['log2_aFC_c0_min_95_interv'] = np.full(len(eqtl_dataf.index.values),np.nan) 
        eqtl_dataf['log2_aFC_c0_plus_95_interv'] = np.full(len(eqtl_dataf.index.values),np.nan) 
    
    return eqtl_dataf


def nocovar_linear_estimate(h1,h2,expr):

        """
        Make an initial guess at the aFCs, using a linear fit, without a covariate matrix
        """

        
        #take only the samples that are in the haplotype matrix - this is the order the 
        #columns are organized in
        expr=expr[h1.columns].values
        added_cov = h1.astype(float) + h2.astype(float)
        X = (added_cov) * 0.5

        #do the linear regression
        reg = LinearRegression().fit(X.transpose(),expr.transpose())
        sa = reg.coef_.transpose() #each coeff is for a different variant
        C0 = reg.intercept_ #the error

        sad = {}
        varc = 0
        variant_list = list(h1.index.tolist())
        for _ in variant_list:

            sad.update( {variant_list[varc] : sa[varc][0]} )
            #sad[var] = sa[varc]
            varc += 1
            
        #return the aFCs as well as the C0
        return [sad,C0]
    
def linear_estimate(h1,h2,expr,cov_dataf):

        """
        Make an initial guess at the aFCs, using a linear fit, with a covariate matrix
        """
        #take only the samples that are in the haplotype matrix - this is the order the 
        #columns are organized in

        used_inds = list(set(h1.columns.tolist()) & set(cov_dataf.columns.tolist()) & set(expr.columns.tolist()))
        #expr=expr[h1.columns].values
        expr=expr[used_inds].values
        #cov_dataf = cov_dataf[h1.columns]
        cov_dataf = cov_dataf[used_inds].astype(float)
        h1 = h1[used_inds]
        h2 = h2[used_inds]
        #cov_dataf = cov_dataf.append(h1 + h2)
        added_hs = h1 + h2
        cov_dataf = pd.concat([cov_dataf, added_hs], axis=0)

        added_cov = cov_dataf
        X = (added_cov) * 0.5



        #print(added_cov[added_cov.isna()])
        #get the intercepts,etc out of the linear fit
        reg = LinearRegression().fit(X.transpose(), expr.transpose())
        sa = reg.coef_ #the last one is the coef for the variant
        C0 = reg.intercept_ #the error
            
        #now store the results in the same matrix, to be passed to the least square optimizer
        cov_dataf = cov_dataf.transpose().append(pd.DataFrame(sa.transpose(), columns=['covar_coeff'], index=cov_dataf.index).transpose()).transpose()

        cov_dataf = cov_dataf.append(pd.DataFrame(np.full((cov_dataf.shape[1],1), C0), columns=['C0'], index=cov_dataf.columns).transpose()) #append C0s


        return cov_dataf

    
    
def nocovar_least_squares(h1, h2, thisgene_expressions, eqtl_dataf, sa, c0, args):
    
    '''
    Do least squares fit based on initial values from a linear fit for expression and haplotype
    data in a gene, without a covariate matrix
    '''
    
    #get the name of the current gene being processed
    gene = thisgene_expressions.index.tolist()[0]

    params = Parameters()

    variant_list = h1.index.tolist()
    for variant in variant_list:

        params.add( str(variant), value=sa[variant] , min=-6.64, max=6.64)
        
    params.add('C0', value=c0[0])

    afc_minimizer = lmfit.Minimizer(nocovar_fcn2fit, params, fcn_args=(h1, h2, variant_list, thisgene_expressions), xtol=0.01, ftol=0.01)
    result = afc_minimizer.minimize( method='leastsq')


    # lets also use this to disable the option to calculate conf intervals if the user doesnt want them 
    ci_unsuccessful = 1
    if args.conf == True:
	    
            # this skips appending the confidence intervals if they couldnt be calculated
            try:

                ci = lmfit.conf_interval(afc_minimizer, result)
                ci_unsuccessful = 0

            except: 

                ci_unsuccessful = 1
		


    zero_list = {} #array of eqtls that will have an aFC of 0, will be empty
                   #if we dont need to redo the fit

    #lets check if we have to redo the fit
    if args.redo == True:

        if ci_unsuccessful == 0:

            refit_params = Parameters()

            #filter out eqtls based on their confidence intervals
            for param in result.params:

                #check if the variants aFC should be set to 0
                if param != 'C0':

                    minus_95 = ci[param][1][1]
                    plus_95 = ci[param][5][1]

                    #if the ci contains 0, add eqtl to zero_list
                    if minus_95 <= 0 <= plus_95:

                        zero_list[param] = [params[param].value, params[param].stderr, minus_95, plus_95]
                        #print(str(minus_95) + " and " + str(plus_95) + " have 0 between them")
                    
                #remove the parameter from the fit we are about to do
                if param not in zero_list:
                    
                    if param != 'C0':

                        refit_params.add( str(param), value=params[param].value , min=-6.64, max=6.64)

                    #add C0 to parameters as it is
                    else:

                        refit_params.add( str(param), value=params[param].value )
                         

            #redo the fit with variants that are not low quality
            remaining_variants = list(set(variant_list).symmetric_difference(set(list(zero_list.keys()))))
            #print("Remaining variants for the fit:")
            #print(remaining_variants)
            #print("Out of original variants:")
            #print(variant_list)
            #print("Excluded variants:")
            #print(zero_list)
            if len(remaining_variants) > 0:
                
                afc_minimizer = lmfit.Minimizer(nocovar_fcn2fit, refit_params, fcn_args=(h1, h2, remaining_variants, thisgene_expressions), xtol=0.01, ftol=0.01)
                #rewrite the previous fits
                result = afc_minimizer.minimize( method='leastsq')
                ci = lmfit.conf_interval(afc_minimizer, result)

        else:       

            raise Exception("Confidence intervals were not set for calculation or an error occurred while calculating them")



    #save parameters and results in results df
    for parameter in result.params:
       
        #if parameter is an s value, variant_order excludes variants with nan values
        if parameter in variant_list:
            
       
            eqtl_dataf.log2_aFC[(parameter == eqtl_dataf.variant_id_clean) & (eqtl_dataf.gene_id_clean == gene)] = result.params[parameter].value
            eqtl_dataf.log2_aFC_error[(parameter == eqtl_dataf.variant_id_clean) & (eqtl_dataf.gene_id_clean == gene)] = result.params[parameter].stderr
       
            if ci_unsuccessful == 0:
       
                #get +-95% intervals
                minus_95 = ci[parameter][1][1]
                eqtl_dataf.log2_aFC_min_95_interv[(parameter == eqtl_dataf.variant_id_clean) & (eqtl_dataf.gene_id_clean == gene)] = minus_95
       
                plus_95 = ci[parameter][5][1]
                eqtl_dataf.log2_aFC_plus_95_interv[(parameter == eqtl_dataf.variant_id_clean) & (eqtl_dataf.gene_id_clean == gene)] = plus_95
       
       
        else: #its C0
            eqtl_dataf.log2_aFC_c0.loc[(eqtl_dataf.gene_id_clean == gene)] = result.params[parameter].value
       
            if ci_unsuccessful == 0:
       
                #get +-95% intervals
                c_minus_95 = ci[parameter][1][1]
                eqtl_dataf.log2_aFC_c0_min_95_interv[eqtl_dataf.gene_id_clean == gene] = c_minus_95

                c_plus_95 = ci[parameter][5][1]
                eqtl_dataf.log2_aFC_c0_plus_95_interv[eqtl_dataf.gene_id_clean == gene] = c_plus_95

       
    

    #now save just the eqtls that have a conf interval spanning 0
    for lowq_eqtl in zero_list.keys():
 
        # No need for confidence intervals
        eqtl_dataf.log2_aFC[( lowq_eqtl == eqtl_dataf.variant_id_clean) & (eqtl_dataf.gene_id_clean == gene)] = 0

    return eqtl_dataf



def least_squares(h1, h2, thisgene_expressions, eqtl_dataf, covar, args):
    
    '''
    Do least squares fit based on initial values from a linear fit for expression and haplotype
    data in a gene, with a covariate matrix
    '''    
    
    #get the name of the current gene being processed
    gene = thisgene_expressions.index.tolist()[0]
    
    params = Parameters()
    
    variant_list = h1.index.tolist()
    for parameter in covar['covar_coeff'].index:
        
        # lets restrict the value for C0 and the variant covariants
        if parameter in variant_list:
                    
            params.add( str(parameter), value=covar['covar_coeff'][parameter], min=-6.64, max=6.64)
            
        elif parameter == 'C0':
            
            params.add( str(parameter), value=covar['covar_coeff'][parameter])

        # or if these are just covariate parameters
        else:
        
            params.add( str(parameter), value=covar['covar_coeff'][parameter])
    
    afc_minimizer = lmfit.Minimizer(fcn2fit, params, fcn_args=(h1, h2, variant_list, thisgene_expressions, covar), xtol=0.01, ftol=0.01)
    result = afc_minimizer.minimize( method='leastsq')

    # lets also use this to disable the option to calculate conf intervals if the user doesnt want them
    ci_unsuccessful = 1
    if args.conf == True:

            # this skips appending the confidence intervals if they couldnt be calculated
            try:

                ci = lmfit.conf_interval(afc_minimizer, result)
                ci_unsuccessful = 0

            except: 

                ci_unsuccessful = 1
	   

    #save parameters and results in results df
    for parameter in result.params:

        #if parameter is an s value, variant_order excludes variants with nan values
        if parameter in variant_list:

            eqtl_dataf.log2_aFC[(parameter == eqtl_dataf.variant_id_clean) & (eqtl_dataf.gene_id_clean == gene)] = result.params[parameter].value
            eqtl_dataf.log2_aFC_error[(parameter == eqtl_dataf.variant_id_clean) & (eqtl_dataf.gene_id_clean == gene)] = result.params[parameter].stderr

            if ci_unsuccessful == 0:

                #get +-95% intervals
                minus_95 = ci[parameter][1][1]
                eqtl_dataf.log2_aFC_min_95_interv[(parameter == eqtl_dataf.variant_id_clean) & (eqtl_dataf.gene_id_clean == gene)] = minus_95

                plus_95 = ci[parameter][5][1]
                eqtl_dataf.log2_aFC_plus_95_interv[(parameter == eqtl_dataf.variant_id_clean) & (eqtl_dataf.gene_id_clean == gene)] = plus_95


        else: #its C0 or one of the peer factors
    

            if parameter == 'C0':
                eqtl_dataf.log2_aFC_c0[eqtl_dataf.gene_id_clean == gene] = result.params[parameter]

                if ci_unsuccessful == 0:

                    #get +-95% intervals
                    c_plus_95 = ci[parameter][1][1]
                    eqtl_dataf.log2_aFC_c0_plus_95_interv[eqtl_dataf.gene_id_clean == gene] = c_plus_95
                    c_minus_95 = ci[parameter][5][1]
                    eqtl_dataf.log2_aFC_c0_min_95_interv[eqtl_dataf.gene_id_clean == gene] = c_minus_95

                
    return eqtl_dataf



def make_haplotype_arrays(thisgene_haplotypes):
    
    '''
    Takes the haplotype df that only contains the variants we are interested in for a certain gene
    and make two numpy arrays out of them, once for the ref, the other for the alt haplotype. 
    The left haplotype value is always ref, the right one is always alt.
    '''

    h1 = thisgene_haplotypes.apply(lambda x: x.str.split(':').str[0].str.split('|').str[0]).replace('./.',np.nan).astype(float)
    h2 = thisgene_haplotypes.apply(lambda x: x.str.split(':').str[0].str.split('|').str[1]).replace('./.',np.nan).astype(float)


    return [h1,h2]
    
    
def get_nan_columns(df):
    
    '''Returns a list of columns in the input dataframe that have nan values'''
    
    nan_values = df.isna()
    nan_columns = nan_values.any()
    columns_with_nan = df.columns[nan_columns].tolist()

    return columns_with_nan
    
    
def clean_matrix_nans(h1, h2, expression_matrix):
    
    '''
    This take in the expression and haplotype matrices and drop individuals that have nans anywhere
    and returns all 3 input dataframes without the offending individuals
    '''
    h1_nan_cols = get_nan_columns(h1)
    h2_nan_cols = get_nan_columns(h2)
    expr_nan_cols = get_nan_columns(expression_matrix)
    
    #get a list of columns that need to be dropped (union)
    offending_inds = list(set(h1_nan_cols) | set(h2_nan_cols) | set(expr_nan_cols))
        
    h1 = h1.drop(columns=offending_inds)
    h2 = h2.drop(columns=offending_inds)
    expression_matrix = expression_matrix.drop(columns=offending_inds)
    
    return [h1, h2, expression_matrix]




    
    
    
