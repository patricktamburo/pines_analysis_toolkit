import pdb 
import numpy as np 
from scipy.stats import pearsonr 
from sklearn import linear_model
import pandas as pd 
from pines_analysis_toolkit.analysis.block_splitter import block_splitter
import matplotlib.pyplot as plt 

def regression(flux, regressors, corr_significance=1e-2, verbose=False):
    '''Authors:
        Patrick Tamburo, Boston University, May 2021
        Purpose: 
            Does a regression of target flux against regressors, if the correlation significance is less than the significance threshold. 
        Inputs:
            flux (np array): Array of flux values.
            regressors (dict): Dictionary of regressors. The length of the regressors must match that of flux (i.e., they have to be taken at the same time as the flux observations).
            corr_significance (float): The correlation significance threhsold. If the significance of a parameter's correlation is *less* than this value, it is used in the regression model. 
            verbose (bool): Whether or not to print info about regressor correlation.
        Outputs:
            corrected_flux (np array): Array of corrected flux values. 
    '''  
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

def leave_one_out_regression(times, flux, flux_err, regressors, corr_significance=1e-2, verbose=False):
    '''Authors:
        Patrick Tamburo, Boston University, May 2021
        Purpose: 
            Does a regression of target flux against regressors, if the correlation significance is less than the significance threshold. 
            This version a one regression for each block of data on an observing night, leaving one block out from the regression each time. 
            The model that results in the lowest chi-squared is retained as the best global model. 
        Inputs:
            flux (np array): Array of flux values.
            regressors (dict): Dictionary of regressors. The length of the regressors must match that of flux (i.e., they have to be taken at the same time as the flux observations).
            corr_significance (float): The correlation significance threhsold. If the significance of a parameter's correlation is *less* than this value, it is used in the regression model. 
            verbose (bool): Whether or not to print info about regressor correlation.
        Outputs:
            corrected_flux (np array): Array of corrected flux values. 
    '''  
    keys = np.array(list(regressors.keys()))
    good_locs = np.where(~np.isnan(flux))[0] #Only perform the regression on non-NaN values. 

    block_inds = block_splitter(times)
    chi_squares = np.zeros(len(block_inds))

    plt.ion()
    fig, ax = plt.subplots(nrows=len(block_inds),ncols=1,sharex=True,sharey=True,figsize=(6,10))

    for b in range(len(block_inds)):

        loo_inds = []#"Leave one out" indices
        for k in range(len(block_inds)):
            if k != b: 
                for j in range(len(block_inds[k])):
                    if block_inds[k][j] in good_locs:
                        loo_inds.append(block_inds[k][j])

        sigs = []
        if verbose:
            print('-------------------------')
            print('{:<11s} | {:>11s}'.format('Regressor', 'Signficance'))
            print('-------------------------')
        for i in range(len(regressors)):
            corr, sig = pearsonr(flux[loo_inds], regressors[keys[i]][loo_inds])
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
            regress_dict[use_keys[i]] = regressors[use_keys[i]][loo_inds]
        
        #Finally, add target flux
        regress_dict['flux'] = flux[loo_inds]

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
            chi_squares[b] = np.sum((y-linear_regression_model)**2/(flux_err[loo_inds]**2))
            print(len(y))
            ax[b].plot(times[loo_inds], y, 'ko')
            ax[b].plot(times[loo_inds], linear_regression_model, 'ro')

            #Divide out the fit. 
            #corrected_flux = flux/linear_regression_model   

        else:
            #print('No regressors used.')
            corrected_flux = flux
    pdb.set_trace()

    return corrected_flux