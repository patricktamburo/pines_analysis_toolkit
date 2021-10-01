import pdb 
import numpy as np 
from scipy.stats import pearsonr 
from sklearn import linear_model
import pandas as pd 
from pines_analysis_toolkit.analysis.block_splitter import block_splitter
import matplotlib.pyplot as plt 

def regression(flux, regressors, corr_significance=1e-2, verbose=False):
    """Does a regression of target flux against regressors, if the correlation significance is less than the significance threshold. 

    :param flux: array of flux values
    :type flux: numpy array
    :param regressors: dictionary containing labeled regressors, each a numpy array of the same length as the flux values
    :type regressors: dict
    :param corr_significance: the p-value that a correlation between the flux and an individual regressor must have to be included in the regression, defaults to 1e-2
    :type corr_significance: float, optional
    :param verbose: whether or not to print information about the regressors used, defaults to False
    :type verbose: bool, optional
    :return: regressed flux
    :rtype: numpy array
    """

    keys = np.array(list(regressors.keys()))
    good_locs = np.where(~np.isnan(flux))[0] #Only perform the regression on non-NaN values. 

    sigs = []
    if verbose:
        print('-------------------------')
        print('{:<11s} | {:>11s}'.format('Regressor', 'Signficance'))
        print('-------------------------')
    for i in range(len(regressors)):
        corr, sig = pearsonr(flux[good_locs], regressors[keys[i]][good_locs])
        sigs.append(sig)
        if verbose:
            print('{:<11s} | {:>.2e}'.format(keys[i], sig))

    use_inds = np.where(np.array(sigs) <= corr_significance)
    use_keys = keys[use_inds]

    if verbose:
        print('Using the following regressors: ')
        for i in range(len(use_keys)):
            print('{:<11s}'.format(use_keys[i]))
        print('')

    #Now set up the linear regression.
    regr = linear_model.LinearRegression()
    
    #Set up the regression dict using only the regressors with correlation significances less than corr_significance
    regress_dict = {}
    for i in range(len(use_keys)):
        regress_dict[use_keys[i]] = regressors[use_keys[i]]
    
    #Finally, add target flux
    regress_dict['flux'] = flux

    #Get list of keys
    keylist = list()
    for i in regress_dict.keys():
        keylist.append(i)

    #Create data frame of regressors.
    df = pd.DataFrame(regress_dict,columns=keylist)
    x = df[keylist[0:len(keylist)-1]]
    y = df['flux']

    if np.shape(x)[1] >0:
        regr.fit(x[~np.isnan(y)],y[~np.isnan(y)])

        #Now, define the model.
        linear_regression_model = regr.intercept_

        for i in range(len(use_keys)):
            linear_regression_model += regr.coef_[i]*regress_dict[use_keys[i]]

        #Calculate the chi^2 of the fit. 
        chi_2 = np.nansum((flux - linear_regression_model)**2/(linear_regression_model))

        #Divide out the fit. 
        corrected_flux = flux/linear_regression_model   

    else:
        #print('No regressors used.')
        corrected_flux = flux
    
    return corrected_flux
