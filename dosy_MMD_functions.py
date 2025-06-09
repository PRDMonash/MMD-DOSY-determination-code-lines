# Author: Igor W.F. Silva -- PRD group -- polymatter.net
# v1.0 - 2024 Nov. 11
# Check our paper in DOI: xxx
# When using this library, please remind to cite us "xxx"
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson

def integrals(coef_d, prob_d, K, V):
    delta_coef_d = np.diff(coef_d) # Normalization of P(D) by area
    area_prob_d_total = simpson(y=prob_d, x=coef_d, dx=delta_coef_d) # area calculated by the Simpson's rule (numerical calculation of areas and volumes)
    # possible to use *from numpy import traz // a =  trapz(y, dx)* as well where the calculation is done by the trapezoidal rule
    prob_d_area_norm = (prob_d / coef_d)/area_prob_d_total #P(D) normalized by the area of P(D) curve + adjusting for the dlogD = dD/D
    
    
    m_i = (coef_d / K) ** -(1 / V) # transforming Coefficient of Diffusion to Molar Mass in i
    delta_m = np.diff(np.linspace(m_i.min(), m_i.max(), len(m_i))).mean()
    w_m = - np.gradient(coef_d, m_i) * prob_d_area_norm #V * K * (m_i ** (-V-1)) * prob_d # weight fraction distribution
    m_n = simpson(y=w_m, x=m_i, dx=delta_m) / simpson(y=(w_m / m_i), x=m_i, dx=delta_m) #sum(w_m) / sum(w_m / m_i) # number average molar mass (Mn) a.k.a. mean
    m_w = simpson(y=(w_m * m_i), x=m_i, dx=delta_m) / simpson(y=w_m, x=m_i, dx=delta_m) #sum(w_m * m_i) / sum(w_m) # mass weight average molar mass (Mw)
    pdi = m_w / m_n # dispersity
    cv = np.sqrt(pdi - 1) # coefficient of variation = standard deviation / mean
        
    # Create a dictionary to store the results
    results = {
        "Mn = ": m_n,
        "Mw = ": m_w,
        "Mw / Mn = ": pdi,
        "Coef. Variation = ": cv
        }
     
    return results

def summation(coef_d, prob_d, K, V):
    m_i = (coef_d / K) ** (-1 / V) # transforming Coefficient of Diffusion to Molar Mass in i
    m_n = sum(prob_d) / sum(prob_d / m_i) # number average molar mass (Mn) a.k.a. mean
    m_w = sum(prob_d * m_i) / sum(prob_d) # mass weight average molar mass (Mw)
    pdi = m_w / m_n # dispersity
    cv = np.sqrt(pdi - 1) # coefficient of variation = standard deviation / mean
        
    # Create a dictionary to store the results
    results_sum = {
        "Mn = ": m_n,
        "Mw = ": m_w,
        "Mw / Mn = ": pdi,
        "Coef. Variation = ": cv
        }
    
    return results_sum


def chromatograms(coef_d, prob_d, error, filename, K, V):
    delta_coef_d = np.diff(coef_d) # Normalization of P(D) by area
    area_prob_d_total = simpson(y=prob_d, x=coef_d, dx=delta_coef_d) # area calculated by the Simpson's rule (numerical calculation of areas and volumes)
    # possible to use *from numpy import traz // a =  trapz(y, dx)* as well where the calculation is done by the trapezoidal rule
    prob_d_area_norm = (prob_d / coef_d)/area_prob_d_total #P(D) normalized by the area of P(D) curve + adjusting for the dlogD = dD/D
    
    
    m_i = (coef_d / K) ** -(1 / V) # transforming Coefficient of Diffusion to Molar Mass in i
    w_m = - np.gradient(coef_d, m_i) * prob_d_area_norm #V * K * (m_i ** (-V-1)) * prob_d # weight fraction distribution
    n_m = w_m / m_i
    w_logm = w_m * m_i
    
    chromatogram_curve = {
        "D": coef_d,
        "P(D)": prob_d,
        "error": error,
        "M": m_i,
        "w(M)": w_m,
        "n(M)": n_m,
        "w(dlogM)": w_logm
    }

    #plot for logM versus w(M) and n(M)
    plt.figure(figsize=(12, 4))
    #plt.subplot(row, column, order)
    plt.subplot(131)
    plt.plot(coef_d, prob_d)
    plt.xscale('log')
    plt.xlabel('D / m2s-1'); plt.ylabel('P(D)')
    plt.subplot(132)
    plt.plot(m_i, w_m)
    plt.xscale('log')
    plt.xlabel('M / g mol-1'); plt.ylabel('w(M)')
    plt.subplot(133)
    plt.plot(m_i, w_m/m_i)
    plt.xscale('log')
    plt.xlabel('M / g mol-1'); plt.ylabel('n(M)')
    plt.suptitle(f'Plotting data to P(D) vs D; w(M) vs M; n(M) vs M for {filename}')

    df_chromatogram = pd.DataFrame(chromatogram_curve) # Create a DataFrame for the chromatogram curves
    df_chromatogram.to_csv(filename, index=False)
    print(f"Chromatogram saved to {filename}")

    return plt.show