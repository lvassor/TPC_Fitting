#!/usr/bin/env python3
"""Code to input and wrangle data, fit models and output model results"""

__appname__ = 'mini_project_wrangling_fitting.py'
__author__ = 'Luke Vassor (ljv3618@ic.ac.uk)'
__version__ = '0.0.1'
__license__ = "M.Sc. CMEE, Imperial College London"
# Date: Nov 2018

# Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import lmfit
from lmfit import Minimizer, Parameters, minimize
from scipy import stats
from datetime import datetime
import time
from collections import OrderedDict

# Data wrangling
data_path = "../Data/BioTraits.csv"
data_frame = pd.read_csv(data_path, encoding = "ISO-8859-1", low_memory = False)

def split_datasets(data_frame, n_points):
    """ Returns a dictionary with each key
    as a TPC ID and its corresponding trait values"""

    print("\033[1;32;40m Obtaining unique curves... \033[0m")
    key_values = data_frame["FinalID"].unique().tolist()
    my_dict = dict.fromkeys(key_values)
    for ID in data_frame["FinalID"].unique():
        my_dict[ID] = [np.array([data_frame["ConTemp_K"][k] for k in data_frame[data_frame["FinalID"]==ID].index]), np.array([data_frame["OriginalTraitValue"][j] for j in data_frame[data_frame["FinalID"]==ID].index]), np.array([data_frame["StandardisedTraitName"][l] for l in data_frame[data_frame["FinalID"]==ID].index])] # returns dictionary with temperature array and trait value array for every curve with curves as keys
        if my_dict[ID][1][np.argmin(my_dict[ID][1])] < 0:
            my_dict[ID][1] = (my_dict[ID][1] + 1) - my_dict[ID][1][np.argmin(my_dict[ID][1])]

        
    new_dict = {k: v for k, v in my_dict.items() if len(np.unique(v[0])) > n_points} # remove any keys which have less than n unique temperature measurements
    print("\033[1;32;40m Unique curves obtained! \033[0m")
    return new_dict

def wrangle_and_split(data_frame, n_points):
    """ Takes dataset, wrangles it, isolates curves which meet
    criteria n_points and returns wrangled data frame"""
    print("\033[1;32;40m Beginning data wrangling... \033[0m")
    subset_1 = data_frame[['FinalID', 'StandardisedTraitName', 'OriginalTraitValue', 'ConTemp']]         # Strips out other columns
    subset_2 = subset_1.dropna().drop_duplicates().sort_values(["FinalID"], ascending = [True]) # removes NAs and rows that are temperature duplicates, then sorts by ID
    subset_2["ConTemp_K"] = (subset_2["ConTemp"] + 273.15)
    my_dict = split_datasets(subset_2, n_points)                                 # Not an ideal solution, easier to convert to dictionary and scale negatives/select for n_points from there, then convert back to a pd and save but requires looping (slow)
    wrangled_data = pd.DataFrame(columns = ["FinalID", "ConTemp_K", "OriginalTraitValue", "StandardisedTraitName"])
    for key, value in my_dict.items():
        temp = pd.DataFrame(value[0])
        trait = pd.DataFrame(value[1])
        name = pd.DataFrame(value[2])
        curve = pd.merge(temp, trait, left_index=True, right_index=True)
        curve = pd.merge(curve, name, left_index=True, right_index=True)
        curve.columns = ["ConTemp_K", "OriginalTraitValue", "StandardisedTraitName"]
        curve["FinalID"] = key
        wrangled_data = pd.concat([wrangled_data, curve], axis = 0)  # ideally want to find better way of wrangling for unique temps inside original frame, as now have a dictionary and two pd frames in memory
    wrangled_data.sort_values(["FinalID"], ascending = [True]).set_index("FinalID").to_csv('../Results/wrangled_data.csv')
    print("\033[1;32;40m Data wrangling complete! \033[0m")
    return my_dict

my_dict = wrangle_and_split(data_frame, 6) # obtain our dictionary dataset

model_names = ["schoolfield","briere","cubic"] ## cant get eear to work, future work. As a package project, a large range of model functions could be scripted
models = [i.lower() for i in model_names] #convert model names to lower case for robustness
fitted_models = [] #Use to store model objects
not_fitted = [] #List to detail any models which fail to fit

# Define model functions
def cubic(params, x, data):
    """ Phenomenological cubic polynomial model"""
    B_0  = params['B_0']
    B_1  = params['B_1']
    B_2  = params['B_2']
    B_3  = params['B_3']
    T = x
    model = B_0 + B_1*T + B_2*T**2 + B_3*T**3

    return model - data # residual array

def briere(params, x, data):
    """ Phenomenological model"""
    B_0 = params['B_0']
    T = x
    T_0 = params['T_0']
    T_m = params['T_m']

    model = B_0*T*(T-T_0)*((T_m-T)**0.5)

    return model - data # residual array

def schoolfield(params, x, data):
    """Mechanistic model"""
    B_0 = params['B_0']
    E = params['E']
    k = params['k']
    T_h = params['T_h']
    E_h = params['E_h']
    T = x

    model = (B_0*np.exp((-E/k)*((1/T)-(1/283.15))))/(1+np.exp((E_h/k)*(1/T_h - 1/T)))

    return model - data # residual array

def log_eaar(params, x, data): # could not achieve this, future work
    """Mechanistic model"""
    A_0 = params['A_0']
    E_b = params['E_b']
    k = params['k']
    E_dh = params['E_dh']
    T_m = params['T_m']
    E_dCp = params['E_dCp']
    T = x

    ln_model = np.log(A_0) - (E_b - E_dh*(1 - T/T_m) - E_dCp*(T - T_m - T*np.log(T/T_m)))/(k*T)
    return ln_model - np.log(data) # residual array

def eaar(params, x, data):  # could not achieve this, future work
    """Mechanistic model"""
    A_0 = params['A_0']
    E_b = params['E_b']
    k = params['k']
    E_dh = params['E_dh']
    T_m = params['T_m']
    E_dCp = params['E_dCp']
    T = np.array(x, dtype=float)

    model = A_0*np.exp(-(E_b - (E_dh*(1 - T/T_m) + E_dCp*(T - T_m - T*np.log(T/T_m))))/k*T)

    return model - data # residual array

# Estimate starting values for parameters
def estimate_parameters(model, key): # need to use regression for schoolfield
    """Estimates starting values for model parameters"""
    T_ref = 283.15
    k = 8.617e-5
    x = my_dict[key][0]                         # temp values
    y = my_dict[key][1]                         # trait values
    y_left_half = y[:np.argmax(y)+1]            # take y data to the left of the peak
    x_left_half = x[:np.argmax(y)+1]            # take x data to the left of the peak
    y_right_half = y[(np.argmax(y)+1):]         # take y data to the right of the peak
    x_right_half = x[(np.argmax(y)+1):]         # take x data to the right of the peak
    T_pk = x[np.argmax(y)]

    # Calculate fixed values
    params = Parameters()
    params.add("k", value = k, vary=False)
    params.add("T_ref", value = T_ref, vary = False)
    params.add("T_pk", value = T_pk, vary = False) # optumum temperature
    params.add("max_response", value = y[np.argmax(y)], vary = False)
    params.add("max_temp", value = max(x), vary = False)
    params.add("n_left", value = len(y_left_half), vary = False)
    params.add("n_right", value = len(y_right_half), vary = False)  

    if model == "cubic":
        # params = Parameters()
        params.add("B_0", value = 0.)
        params.add("B_1", value = 0.)
        params.add("B_2", value = 0.)
        params.add("B_3", value = 0.)
    elif model == "briere":
        # params = Parameters()
        params.add("B_0", value = 0.01) # normalising constant
        params.add("T_0", value = min(x), min = 0, max = T_pk) # lowest temp cant be lower than zero or greater than peak
        params.add("T_m", value = max(x), min = x[4], max = 500) # setting highest temp to peak is too restrictive, so it to intermediate temp, allows model to adapt

    elif model == "eaar":
        params.add('A_0', value = np.exp(20))
        params.add('E_b', value = 0.65)
        params.add('E_dh', value = 2.5)
        params.add('T_m', value = 310)
        params.add('E_dCp', value = 0.1)

    elif model == "schoolfield":
        ## Transforming the data for a schoolfield LEFT OF PEAK
        
        B_0_index = np.abs(x - T_ref).argmin()
        
        y_left_trans = np.log(y_left_half)              # log-transform y data
        x_left_trans = 1/(k*(x_left_half))              # transform x data 

        ## Transforming the data for a schoolfield RIGHT OF PEAK
        y_right_trans = np.log(y_right_half)            # log-transform y data
        x_right_trans = 1/(k*(x_right_half))            # transform x data

        ## Don't require estimation
        params.add("B_0", value = y[B_0_index])
        params.add("T_h", value = max(x), min = T_pk, max = 500) # cant be lower than peak

        # Calculate values from left half of data
        if len(y_left_half) > 2: # only use linregress on 3 points minimum
            try:
                slope_left_trans, intercept_left_trans, r_value_left_trans, p_value_left_trans, std_err_left_trans = stats.linregress(x_left_trans,y_left_trans) # use linear regression to estimate E
                if slope_left_trans > 0.2 and slope_left_trans < 1.2:
                    params.add("E", value = slope_left_trans, min = 0.2, max = 1.2) # Dell et al. 2011
                    params.add("E_h", value = 1.15, min = 0.76, max = 10*slope_left_trans) ## New
                else:
                    params.add("E", value = 0.66, min = 0.2, max = 1.2) # Dell et al. 2011
                    params.add("E_h", value = 1.15, min = 0.76, max = 6.6) ## use 2 * E for E_h
            except Exception:
                params.add("E", value = 0.66, min = 0.2, max = 1.2)
                params.add("E_h", value = 1.15, min = 0.76, max = 6.6)
        else:
            params.add("E", value = 0.66, min = 0.2, max = 1.2)
            params.add("E_h", value = 1.15, min = 0.76, max = 6.6)

        # safeguard for nan values being produced
        if np.isnan(params.valuesdict()["E"]):
            params.add("E", value = 0.66, min = 0.2, max = 1.2)
        if np.isnan(params.valuesdict()["E_h"]):
            params.add("E_h", value = 1.15, min = 0.76, max = 6.6)
        if np.isnan(params.valuesdict()["T_h"]):
            params.add("T_h", value = max(x), min = T_pk, max = 500)
    return params

def randomise_estimates(model, key, params):
    "Randomises estimates for model parameter values using random gaussian fluctuation"
    if model == "schoolfield":
        params.add("E", value = np.random.normal(loc=0.66, scale = 1), min = 0.001, max = 10) # if model cannot converge, taking random estimates with default E as mean and try for any growth response using relaxed bounds from Dell et al.
        params.add("T_h", value = np.random.normal(loc=params.valuesdict()["max_temp"], scale = 1), min = params.valuesdict()["T_pk"])
        params.add("E_h", value = np.random.normal(loc = 1.15, scale = 1), min = 0.001, max = 10)
    if model == "briere":
        B_0 = params.add("B_0", value = np.random.normal(loc=params.valuesdict()["B_0"]))
        T_0 = params.add("T_0", value = np.random.normal(loc=params.valuesdict()["T_0"]))
        T_m = params.add("T_m", value = np.random.normal(loc=params.valuesdict()["T_m"]))

    return params


## Attempt to fit SS to all

def fit_schoolfield(dictionary, key, columns):
    """Fits schoolfield model to input data from dictionary"""
    successful_fits = 0
    unsuccessful_fits = 0
    attempts = 0
    x = my_dict[key][0]
    y = my_dict[key][1]
    params = estimate_parameters(model = "schoolfield", key = key)
    try:         
        school_min = minimize(schoolfield, params, args = (my_dict[key][0], my_dict[key][1])) # try to minimize based on estimates
    except ValueError:
        timeout = time.time() + 5 # limit attempts to 5 seconds
        while attempts <= 30 and time.time() < timeout: 
            params = randomise_estimates(model="schoolfield", key=key, params = params)
            try:                            
                school_min = minimize(schoolfield, params, args = (my_dict[key][0], my_dict[key][1]))
            except Exception:
                if attempts == 30:
                    unsuccessful_fits = unsuccessful_fits + 1
                    school_min = False
                    break
                attempts = attempts + 1
                continue
    sch_results = OrderedDict.fromkeys(columns)
    if school_min:
        rss_sch = np.sum(school_min.residual**2)
        tss_sch = np.sum((y - np.mean(y))**2)
        rsq_sch = 1 - (rss_sch/tss_sch)
        # attempts_rsq = 0                      ### LEFT THIS IN FOR OPTIONAL USE, INCREASES TIME BUT BETTER FITS
        # rsq_timer = time.time() + 1
        # while rsq_sch < 0.5 or time.time() < rsq_timer:
        #     # print("inside while loop")
        #     try: 
        #         # "trying to improve score"
        #         params = randomise_estimates(model="schoolfield", key=key, params = params)
        #         school_min_try = minimize(schoolfield, params, args = (my_dict[key][0], my_dict[key][1]))
        #         rss_sch = np.sum(school_min.residual**2)
        #         tss_sch = np.sum((y - np.mean(y))**2)
        #         if 1 - (rss_sch/tss_sch) > rsq_sch:
        #             school_min = school_min_try
        #         rsq_sch = 1 - (rss_sch/tss_sch)
        #     except Exception:
        #         continue
        #     attempts_rsq = attempts_rsq + 1
        #     # print(attempts_rsq)
        #     if attempts_rsq == 30:
        #         break
        # # print("outside while loop")
        # successful_fits = successful_fits + 1

        school_pars = school_min.params.valuesdict() # obtain dictionary of values
        sch_results["sch_final_B_0"] = school_pars["B_0"]
        sch_results["sch_final_E"]= school_pars["E"]
        sch_results["sch_final_E_h"]= school_pars["E_h"]
        sch_results["sch_final_T_h"]= school_pars["T_h"]
        sch_results["sch_R_Squared"]= rsq_sch
        sch_results["sch_BIC"]= school_min.bic
        sch_results["sch_AICc"]= school_min.aic + ((2*(4)**2 + 2*(4))/(school_min.ndata - 4 - 1))
    else:  # if model doesnt converge values are nan
        sch_results["sch_final_E"]= np.nan
        sch_results["sch_final_B_0"]= np.nan
        sch_results["sch_final_T_h"]= np.nan
        sch_results["sch_final_E_h"]= np.nan
        sch_results["sch_R_Squared"]= np.nan
        sch_results["sch_AICc"]= np.nan
        sch_results["sch_BIC"]= np.nan

    # these are start values dont require model to converge    
    sch_results["FinalID"] = key
    sch_results["sch_start_B_0"]= params.valuesdict()["B_0"]     # starting values
    sch_results["sch_start_E"]= params.valuesdict()["E"]         # starting values
    sch_results["sch_start_E_h"]= params.valuesdict()["E_h"]     # starting values
    sch_results["sch_start_T_h"]= params.valuesdict()["T_h"]     # starting values
    sch_results["Est_Tpk"]= params.valuesdict()["T_pk"]
    sch_results["Max_response"]= params.valuesdict()["max_response"]
    sch_results["Number_Of_Data_Points"]= len(x)
    sch_results["Points_Before_Peak"]= params.valuesdict()["n_left"]
    sch_results["Points_After_Peak"]= params.valuesdict()["n_right"]
    
    return sch_results

def fit_briere(dictionary, key, columns):
    """ Fits briere model to input data from dictionary"""
    successful_fits = 0
    unsuccessful_fits = 0
    attempts = 0
    x = my_dict[key][0]
    y = my_dict[key][1]
    params = estimate_parameters(model = "briere", key = key) 
    try:                    
        briere_min = minimize(briere, params, args = (my_dict[key][0], my_dict[key][1]))
    except ValueError:
        timeout = time.time() + 5
        while attempts <= 30 and time.time() < timeout:
            params = randomise_estimates(model = "briere", key=key, params = params)
            try:                            
                briere_min = minimize(briere, params, args = (my_dict[key][0], my_dict[key][1]))
            except Exception:
                if attempts == 30:
                    unsuccessful_fits = unsuccessful_fits + 1
                    briere_min = None
                    break
                attempts = attempts + 1
                continue
    bri_results = OrderedDict.fromkeys(columns) 
    if briere_min:
        briere_pars = briere_min.params.valuesdict()
        successful_fits = successful_fits + 1
        rss_bri = np.sum(briere_min.residual**2)
        tss_bri = np.sum((y - np.mean(y))**2)
        rsq_bri = 1 - (rss_bri/tss_bri)
        # attempts_rsq = 0
        # rsq_timer = time.time() + 1
        # while rsq_bri < 0.5 or time.time() < rsq_timer:
        #     # print("inside while loop")
        #     try: 
        #         params = randomise_estimates(model="briere", key=key, params = params)
        #         briere_min_try = minimize(briere, params, args = (my_dict[key][0], my_dict[key][1]))
        #         rss_bri = np.sum(briere_min.residual**2)
        #         tss_bri = np.sum((y - np.mean(y))**2)
        #         if 1 - (rss_bri/tss_bri) > rsq_bri:
        #             briere_min = briere_min_try
        #         rsq_bri = 1 - (rss_bri/tss_bri)
        #     except Exception:
        #         continue
        #     attempts_rsq = attempts_rsq + 1
        #     if attempts_rsq == 30:
        #         break
        
        bri_results["bri_final_B_0"] = briere_pars["B_0"]
        bri_results["bri_final_T_0"] = briere_pars["T_0"]
        bri_results["bri_final_T_m"] = briere_pars["T_m"]
        bri_results["bri_R_Squared"] = rsq_bri
        bri_results["bri_AICc"] = briere_min.aic + ((2*(3)**2 + 2*(3))/(briere_min.ndata - 3 - 1))
        bri_results["bri_BIC"] = briere_min.bic
        
    else:
        bri_results["bri_final_B_0"] = np.nan
        bri_results["bri_final_T_0"] = np.nan
        bri_results["bri_final_T_m"] = np.nan
        bri_results["bri_R_Squared"] = np.nan
        bri_results["bri_AICc"] = np.nan
        bri_results["bri_BIC"] = np.nan
    bri_results["FinalID"] = key 
    bri_results["bri_start_B_0"] = params.valuesdict()["B_0"]
    bri_results["bri_start_T_0"] = params.valuesdict()["T_0"]
    bri_results["bri_start_T_m"] = params.valuesdict()["T_m"]
    bri_results["Est_Tpk"] = params.valuesdict()["T_pk"]
    bri_results["Max_response"] = params.valuesdict()["max_response"]
    bri_results["Number_Of_Data_Points"] = len(x)
    bri_results["Points_After_Peak"] = params.valuesdict()["n_right"]
    bri_results["Points_Before_Peak"] = params.valuesdict()["n_left"]
    return bri_results


def fit_eaar(dictionary, key, columns):
    """Fits eaar model to input data from dictionary""" # does not work so won't be called -requires further work
    eaar_results = []
    successful_fits = 0
    unsuccessful_fits = 0
    attempts = 0
    x = my_dict[key][0]
    y = my_dict[key][1]
    # print(key)
    params = estimate_parameters(model = "eaar", key = key)
    print(params)
    try:                    
        log_eaar_min = minimize(log_eaar, params, args = (my_dict[key][0], my_dict[key][1]))
        params.add("A_0", value = log_eaar_min.params.valuesdict()["A_0"])
        params.add("E_b", value = log_eaar_min.params.valuesdict()["E_b"])
        params.add("E_dh", value = log_eaar_min.params.valuesdict()["E_dh"])
        params.add("T_m", value = log_eaar_min.params.valuesdict()["T_m"])
        params.add("E_dCp", value = log_eaar_min.params.valuesdict()["E_dCp"])
        eaar_min = minimize(eaar, params, args = (my_dict[key][0], my_dict[key][1]))   
        print("estimated parameters")
        # print("tried briere_min")
    except ValueError:
        # timeout = time.time() + 10
        print("could not estimate parameters")
        # eaar_min = None
        
        # while attempts <= 10 and time.time() < timeout:
        #     params = randomise_estimates(model = "briere", key=key, params = params)
        #     try:                            
        #         briere_min = minimize(briere, params, args = (my_dict[key][0], my_dict[key][1]))
        #         # print("worked on first attempt")
        #     except Exception:
        #         if attempts == 10:
        #             # print("unsuccessful after 10 attempts, moving to next key...")
        #             unsuccessful_fits = unsuccessful_fits + 1
        #             briere_min = None
        #             break
        #         # print(attempts)
        #         # print("here it is")
        #         attempts = attempts + 1
        #         continue
    # print("worked without randomising")
    # print(params.valuesdict()["n_right"])
    # print(params.valuesdict()["n_left"])
    eaar_results = OrderedDict.fromkeys(columns) 
    eaar_results["FinalID"] = key
    eaar_results["eaar_start_A_0"] = params.valuesdict()['A_0']
    eaar_results["eaar_start_E_b"] = params.valuesdict()['E_b']
    eaar_results["eaar_start_E_dh"] = params.valuesdict()['E_dh']
    eaar_results["eaar_start_T_m"] = params.valuesdict()['T_m']
    eaar_results["eaar_start_E_dCp"] = params.valuesdict()['E_dCp']    
    eaar_results["Est_Tpk"] = params.valuesdict()["T_pk"]
    eaar_results["Max_response"] = params.valuesdict()["max_response"]
    eaar_results["Number_Of_Data_Points"] = len(x)
    eaar_results["Points_After_Peak"] = params.valuesdict()["n_right"]
    eaar_results["Points_Before_Peak"] = params.valuesdict()["n_left"]
    if eaar_min:
        eaar_pars = eaar_min.params.valuesdict()
        successful_fits = successful_fits + 1
        rss_eaar = np.sum(eaar_min.residual**2)
        tss_eaar = np.sum((y - np.mean(y))**2)
        rsq_eaar = 1 - (rss_eaar/tss_eaar)
        eaar_results["eaar_final_A_0"] = eaar_pars["A_0"]
        eaar_results["eaar_final_E_b"] = eaar_pars["E_b"]
        eaar_results["eaar_final_E_dh"] = eaar_pars["E_dh"]
        eaar_results["eaar_final_T_m"] = eaar_pars["T_m"]
        eaar_results["eaar_final_E_dCp"] = eaar_pars["E_dCp"]         
        eaar_results["eaar_R_Squared"] = rsq_eaar
        eaar_results["eaar_AICc"] = eaar_min.aic + ((2*(3)**2 + 2*(3))/(eaar_min.ndata - 3 - 1))
        eaar_results["eaar_BIC"] = eaar_min.bic
        
    else:
        eaar_results["eaar_final_A_0"] = np.nan
        eaar_results["eaar_final_E_b"] = np.nan
        eaar_results["eaar_final_E_dh"] = np.nan
        eaar_results["eaar_final_T_m"] = np.nan
        eaar_results["eaar_final_E_dCp"] = np.nan
        eaar_results["eaar_R_Squared"] = np.nan
        eaar_results["eaar_AICc"] = np.nan
        eaar_results["eaar_BIC"] = np.nan

    return eaar_results


def fit_cubic(dictionary, key, columns):
    """Fits cubic polynomial model to input data from dictionary"""
    cubic_results = []
    successful_fits = 0
    unsuccessful_fits = 0
    attempts = 0
    x = my_dict[key][0]
    y = my_dict[key][1]
    params = estimate_parameters(model = "cubic", key = key)
    try:                    
        cubic_min = minimize(cubic, params, args = (my_dict[key][0], my_dict[key][1]))

    except ValueError:
        timeout = time.time() + 5
        while attempts <= 30 and time.time() < timeout:
            params = randomise_estimates(model = "cubic", key=key, params = params)
            try:                            
                cubic_min = minimize(cubic, params, args = (my_dict[key][0], my_dict[key][1]))
            except Exception:
                if attempts == 30:
                    unsuccessful_fits = unsuccessful_fits + 1
                    cubic_min = None
                    break
                attempts = attempts + 1
                continue
    cubic_results = OrderedDict.fromkeys(columns) 
        
    if cubic_min:
        cubic_pars = cubic_min.params.valuesdict()
        successful_fits = successful_fits + 1
        rss_cub = np.sum(cubic_min.residual**2)
        tss_cub = np.sum((y - np.mean(y))**2)
        rsq_cub = 1 - (rss_cub/tss_cub)
        # attempts_rsq = 0
        # rsq_timer = time.time() + 1
        # while rsq_cub < 0.5 and time.time() < rsq_timer:
        #     # print("inside while loop")
        #     try: 
        #         params = randomise_estimates(model="cubic", key=key, params = params)
        #         cubic_min_try = minimize(cubic, params, args = (my_dict[key][0], my_dict[key][1]))
        #         rss_cub = np.sum(cubic_min.residual**2)
        #         tss_cub = np.sum((y - np.mean(y))**2)
        #         if 1 - (rss_cub/tss_cub) > rsq_cub:
        #             cubic_min = cubic_min_try
        #         rsq_cub = 1 - (rss_cub/tss_cub)
        #     except Exception:
        #         continue
        #     attempts_rsq = attempts_rsq + 1
        #     if attempts_rsq == 30:
        #         break
        cubic_results["cub_final_B_0"] = cubic_pars["B_0"]
        cubic_results["cub_final_B_1"] = cubic_pars["B_1"]
        cubic_results["cub_final_B_2"] = cubic_pars["B_2"]
        cubic_results["cub_final_B_3"] = cubic_pars["B_3"]
        cubic_results["cub_R_Squared"] = rsq_cub
        cubic_results["cub_AICc"] = cubic_min.aic + ((2*(4)**2 + 2*(4))/(cubic_min.ndata - 4 - 1))
        cubic_results["cub_BIC"] = cubic_min.bic
        
    else:
        cubic_results["cub_final_B_0"] = np.nan
        cubic_results["cub_final_B_1"] = np.nan
        cubic_results["cub_final_B_2"] = np.nan
        cubic_results["cub_final_B_3"] = np.nan
        cubic_results["cub_R_Squared"] = np.nan
        cubic_results["cub_AICc"] = np.nan
        cubic_results["cub_BIC"] = np.nan
    cubic_results["FinalID"] = key 
    cubic_results["cub_start_B_0"] = params.valuesdict()["B_0"]
    cubic_results["cub_start_B_1"] = params.valuesdict()["B_1"]
    cubic_results["cub_start_B_2"] = params.valuesdict()["B_2"]
    cubic_results["cub_start_B_3"] = params.valuesdict()["B_3"]
    cubic_results["Est_Tpk"] = params.valuesdict()["T_pk"]
    cubic_results["Max_response"] = params.valuesdict()["max_response"]
    cubic_results["Number_Of_Data_Points"] = len(x)
    cubic_results["Points_After_Peak"] = params.valuesdict()["n_right"]
    cubic_results["Points_Before_Peak"] = params.valuesdict()["n_left"]
    return cubic_results


def calculate_delta_AIC(model, results_dict):
    """Calculates delta AIC values depending on which models have converged or not""" # this is required because can only calculate delta AIC from models which have converged, so have to check which are present
    if model == "schoolfield":
        if np.isnan(results_dict["sch_AICc"]) == True:
            delta_AICc = np.nan
        elif np.isnan(results_dict["bri_AICc"]) == False and np.isnan(results_dict["cub_AICc"]) == False: # when sch AIC is not an nan and neither are the others
            AIC_min = min(results_dict["sch_AICc"], results_dict["bri_AICc"], results_dict["cub_AICc"])
            delta_AICc = results_dict["sch_AICc"] - AIC_min
        elif np.isnan(results_dict["bri_AICc"]) == True and np.isnan(results_dict["cub_AICc"]) == False:
            AIC_min = min(results_dict["sch_AICc"], results_dict["cub_AICc"]) # if briere has no AICc, only use school and cubic
            delta_AICc = results_dict["sch_AICc"] - AIC_min
        elif np.isnan(results_dict["bri_AICc"]) == False and np.isnan(results_dict["cub_AICc"]) == True:
            AIC_min = min(results_dict["sch_AICc"], results_dict["bri_AICc"]) # if cubic has no AICc, only use school and briere
            delta_AICc = results_dict["sch_AICc"] - AIC_min
    if model == "briere":
        if np.isnan(results_dict["bri_AICc"]):
            delta_AICc = np.nan
        elif np.isnan(results_dict["sch_AICc"]) == False and np.isnan(results_dict["cub_AICc"]) == False: # when bri AIC is not an nan and neither are the others
            AIC_min = min(results_dict["sch_AICc"], results_dict["bri_AICc"], results_dict["cub_AICc"])
            delta_AICc = results_dict["bri_AICc"] - AIC_min
        elif np.isnan(results_dict["sch_AICc"]) == True and np.isnan(results_dict["cub_AICc"]) == False:
            AIC_min = min(results_dict["bri_AICc"], results_dict["cub_AICc"]) # if sch has no AICc, only use briere and cubic
            delta_AICc = results_dict["bri_AICc"] - AIC_min
        elif np.isnan(results_dict["sch_AICc"]) == False and np.isnan(results_dict["cub_AICc"]) == True:
            AIC_min = min(results_dict["sch_AICc"], results_dict["bri_AICc"]) # if cubic has no AICc, only use school and briere
            delta_AICc = results_dict["bri_AICc"] - AIC_min
    if model == "cubic":
        if np.isnan(results_dict["cub_AICc"]):
            delta_AICc = np.nan
        elif np.isnan(results_dict["sch_AICc"]) == False and np.isnan(results_dict["bri_AICc"]) == False: # when cub AIC is not an nan and neither are the others
            AIC_min = min(results_dict["sch_AICc"], results_dict["bri_AICc"], results_dict["cub_AICc"])
            delta_AICc = results_dict["cub_AICc"] - AIC_min
        elif np.isnan(results_dict["sch_AICc"]) == True and np.isnan(results_dict["bri_AICc"]) == False:
            AIC_min = min(results_dict["bri_AICc"], results_dict["cub_AICc"]) # if sch has no AICc, only use briere and cubic
            delta_AICc = results_dict["cub_AICc"] - AIC_min
        elif np.isnan(results_dict["sch_AICc"]) == False and np.isnan(results_dict["bri_AICc"]) == True:
            AIC_min = min(results_dict["sch_AICc"], results_dict["cub_AICc"]) # if cubic has no AICc, only use school and briere
            delta_AICc = results_dict["cub_AICc"] - AIC_min
    return delta_AICc

def mergeDictsOverwriteEmpty(d1, d2, d3, main_cols): # alter depending on size of model list
    """ merges results from all dictionaries into global dictionary
    which is then concatenated to pd frame"""
    all_results = OrderedDict.fromkeys(main_cols)
    for key in d1:
        if d1[key] != None:
            all_results[key] = d1[key]
    for key in d2:
        if d2[key] != None:
            all_results[key] = d2[key]
    for key in d3:
        if d3[key] != None:
            all_results[key] = d3[key]
    # for key in d4:
    #     if d4[key] != None:
    #         all_results[key] = d4[key]
    return all_results


def fit_models(models, dictionary):
    """Fits models and output results to pandas data frame"""
    print("\033[1;32;40m Beginning model fitting... \033[0m")
    main_cols = ["FinalID", "sch_start_B_0", "sch_start_E", "sch_start_E_h", "sch_start_T_h", 
                "bri_start_B_0", "bri_start_T_0", "bri_start_T_m", 
                "cub_start_B_0", "cub_start_B_1", "cub_start_B_2", "cub_start_B_3", 
                # "eaar_start_A_0", "eaar_start_E_b", "eaar_start_E_dh", "eaar_start_T_m", "eaar_start_E_dCp",
                "sch_final_B_0", "sch_final_E", "sch_final_E_h", "sch_final_T_h",
                "bri_final_B_0", "bri_final_T_0", "bri_final_T_m", 
                "cub_final_B_0", "cub_final_B_1", "cub_final_B_2", "cub_final_B_3",
                # "eaar_final_A_0", "eaar_final_E_b", "eaar_final_E_dh", "eaar_final_T_m", "eaar_final_E_dCp",
                "Est_Tpk", "Max_response", 
                "sch_R_Squared", "sch_AICc", "sch_BIC",
                "bri_R_Squared", "bri_AICc", "bri_BIC",
                "cub_R_Squared", "cub_AICc", "cub_BIC",
                # "eaar_R_Squared", "eaar_AICc", "eaar_BIC",
                "sch_delta_AICc", "bri_delta_AICc", "cub_delta_AICc",
                "sch_Wi_AICc", "bri_Wi_AICc", "cub_Wi_AICc",
                "Number_Of_Data_Points", "Points_Before_Peak", "Points_After_Peak"] #, "Number_Of_Variables"]
    results = pd.DataFrame(columns = main_cols) # would be good to find a better way to make the index the key as you move through the for loop, could use set_index and make columns of keys
    for key in my_dict.keys(): # loops through curves
        # print(key)
        # print(key)
        for model in models: # loops through models within each curve
            
            if model == 'schoolfield':
                sch_results = fit_schoolfield(dictionary, key, main_cols) # could make a global dictionary which I ammend with each model function, but this would need to be returned and passed onto next function, which defeats point of iterating through what ever models the user has chosen
            if model == 'briere':
                bri_results = fit_briere(dictionary, key, main_cols)
            if model == 'cubic':
                cub_results = fit_cubic(dictionary, key, main_cols)
            if model == 'eaar':
                eaar_results = fit_eaar(dictionary, key, main_cols)
        
        # Merge dictionary results from each model
        all_results = mergeDictsOverwriteEmpty(sch_results, bri_results, cub_results, main_cols) #adjust for more models
        # Calculate delta AICc
        all_results["sch_delta_AICc"] = calculate_delta_AIC(model = "schoolfield", results_dict = all_results)
        all_results["bri_delta_AICc"] = calculate_delta_AIC(model = "briere", results_dict = all_results)
        all_results["cub_delta_AICc"] = calculate_delta_AIC(model = "cubic", results_dict = all_results)

        # Calculate Aikake weights
        all_results["sch_Wi_AICc"] = np.exp(-0.5*all_results["sch_delta_AICc"])/(np.exp(-0.5*all_results["sch_delta_AICc"])+np.exp(-0.5*all_results["bri_delta_AICc"])+np.exp(-0.5*all_results["cub_delta_AICc"]))
        all_results["bri_Wi_AICc"] = np.exp(-0.5*all_results["bri_delta_AICc"])/(np.exp(-0.5*all_results["sch_delta_AICc"])+np.exp(-0.5*all_results["bri_delta_AICc"])+np.exp(-0.5*all_results["cub_delta_AICc"]))
        all_results["cub_Wi_AICc"] = np.exp(-0.5*all_results["cub_delta_AICc"])/(np.exp(-0.5*all_results["sch_delta_AICc"])+np.exp(-0.5*all_results["bri_delta_AICc"])+np.exp(-0.5*all_results["cub_delta_AICc"]))


        all_results_df = pd.DataFrame([all_results], columns = all_results.keys())
        results = pd.concat([results, all_results_df], axis = 0)
    results = results.reindex(columns=main_cols)
    results = results.sort_values(["FinalID"], ascending = [True]).set_index("FinalID")
    print("\033[1;32;40m Model fitting complete! See '../Results' directory. \033[0m")
    return results

starttime = datetime.now() 
results = fit_models(models, my_dict)
results.to_csv('../Results/model_fit_results.csv')
print("\033[1;32;40m Completed in: {}\033[0m".format(datetime.now() - starttime))