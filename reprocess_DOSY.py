import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.integrate import simpson
import dosy_MMD_functions as m2d3

folder = r'G:\My Drive\Research\PhD\DOSY project\Dispersity\Code lines/'
chromatogram_folder = "processed"
os.chdir(folder)


#Define calibration curve parameters v and K logD = logK -v*logM
LOGK = -7.7237 #where LOGK is the logarithm of K on the calibration curve
V = 0.52830; #where V is v from the calibration curve
K = 10 ** LOGK

print('')
print('calibration curve parameter input:')
print('K =', K, '        ','log K =', LOGK)
print('v =', V)


# Create the output directory if it doesn't exist
if not os.path.exists(chromatogram_folder):
    os.makedirs(chromatogram_folder)
    print(f"Created folder: {chromatogram_folder}")

new_data = []


for filename in os.listdir(folder):
    if filename == 'report.csv':
        continue
    if filename.endswith(".csv"):
        full_path = os.path.join(folder, filename)
        
        df = pd.read_csv(filename, usecols=[0, 1, 2])  # guarantee it reads the 3 columns
        d = np.array(df.iloc[:,0]); # matrix of values of D
        p_d = np.array(df.iloc[:,1]); # matrix of the probability of D, a.k.a. intensity, gradient...
        error_pd = np.array(df.iloc[:,2]) # matrix of the error of the probability of D
                
        
        print('')
        print('Filename =', filename)
        print('')
        print('Molar mass distribution calculated from the INTEGRALS of the curve distributions')
        results = m2d3.integrals(d, p_d, K, V)
        for key, value in results.items():
            print(f"{key}: {value}")
        
        print('')
        print('Molar mass distribution calculated from the SUMMATION of intensities')
        results_sum = m2d3.summation(d, p_d, K, V)
        for key, value in results_sum.items():
            print(f"{key}: {value}")
        
        print('')
        
        
        # Append the filename and PDI error to the list
        new_data.append({
            'filename': filename,
            'Mn': results["Mn = "],
            'Mw': results["Mw = "],
            'PDI': results["Mw / Mn = "],
            'CV': results["Coef. Variation = "],
            'Mn (sum)': results_sum["Mn = "],
            'Mw (sum)': results_sum["Mw = "],
            'PDI (sum)': results_sum["Mw / Mn = "],
            'CV (sum)': results_sum["Coef. Variation = "],
              })
        
        processed_filename = filename.replace(".csv", "_processed.csv")
        output_path = os.path.join(chromatogram_folder, processed_filename)
        m2d3.chromatograms(d, p_d, error_pd, output_path, K, V)
       


# Create a DataFrame from the list
df_new_data = pd.DataFrame(new_data)

# Save the DataFrame to a CSV file
output_file = 'report.csv'
df_new_data.to_csv(output_file, index=False)

print(f"Values saved to {output_file}")

plt.show()