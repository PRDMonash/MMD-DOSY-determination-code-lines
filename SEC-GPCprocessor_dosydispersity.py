# This code line is designed for a PSS report
# It might need adaptations for anothers SEC softwares
# Author: Igor W. F. Silva -- PRD group -- polymatter.net
# v0.1 - 2025 Feb. 05
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

import pandas as pd
import numpy as np
from scipy.integrate import simpson
import matplotlib.pyplot as plt
import os

folder = r'G:\My Drive\Research\PhD\DOSY project\Dispersity\gpc/'
os.chdir(folder)

excel_file = 'SEC-GPCprocessor_dosydispersity.xlsx'
file_list_df = pd.read_excel(excel_file)


results = []

# Iterate through each file name in the Excel sheet
for index, row in file_list_df.iterrows():
    file_name = row['File Name']
    if pd.notna(file_name):
        file_name = str(file_name)
        file_path = os.path.join(os.getcwd(), file_name)

        # Process each file
        skip_rows = 0
        with open(file_path, 'r') as file:
            for line in file:
                if 'ELUstart :' in line:
                    break
                skip_rows += 1

        data_lines = []
        with open(file_path, 'r') as file:
            capture = False
            for line in file:
                if 'ELUstart :' in line:
                    capture = True
                    continue  # Skip the 'RAWstart :' line
                if 'ELUstop :' in line:
                    break  # Stop capturing when 'RAWstop :' is found
                if capture:
                    segments = line.strip().split(':')
                    if len(segments) > 0:
                        # Append the whole line instead of just the last segment
                        data_lines.append(line.strip())
                    else:
                        data_lines.append(line.strip())  # If no colon, just append the line

    # Create a DataFrame from the captured lines
    from io import StringIO

    data_str = '\n'.join(data_lines)
    data = StringIO(data_str)

    # Read the data using regular expression for delimiter
    df = pd.read_csv(data, sep=r'\s+', engine='python', skiprows=0, header=None)  # Skip the header row

    # Ensure the 'Time' and 'I1' columns are numerical
    df.iloc[:, 0] = pd.to_numeric(df.iloc[:, 0], errors='coerce')  # Assuming 'Time' is the first column
    df.iloc[:, 1] = pd.to_numeric(df.iloc[:, 1], errors='coerce')  # Assuming 'I1' is the second column
    df.iloc[:, 2] = pd.to_numeric(df.iloc[:, 2], errors='coerce')

    # Drop rows with NaN values in 'Time' or 'I1' columns
    df.dropna(subset=df.columns[[0, 1]], inplace=True)

    # Open the file again to extract calibration coefficients
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Initialize variables to store calibration coefficients
    calibration_coefficients = {}

    in_coeff_section = False

    # Iterate through the lines to find the Calibration Coefficients section
    for line in lines:
        if 'Calibration Coefficients:' in line:
            in_coeff_section = True
            continue  # Skip the header line

        if in_coeff_section:
            if not line.strip():  # Stop if a blank line is encountered
                continue

            if 'ELUstart :' in line or 'ELUstop :' in line:
                # Stop capturing if we encounter 'RAWstart :' or 'RAWstop :' which indicate data sections
                break

            segments = line.strip().split(':')
            if len(segments) > 1:
                key = segments[0].strip()
                value = segments[1].strip()
                calibration_coefficients[key] = value

    # Example: Accessing individual coefficients
    vol_min = float(calibration_coefficients.get('Vol. min', 0.0))  # Default value of 0.0 if key not found
    vol_max = float(calibration_coefficients.get('Vol. max', 0.0))  # Default value of 0.0 if key not found
    r_value = float(calibration_coefficients.get('R', 0.0))  # Default value of 0.0 if key not found
    fit = calibration_coefficients.get('Fit', 'Polynom 5')  # Default value of 'Polynom 5' if key not found

    # Extracting the polynomial coefficients 
    # For a calibration curve: Log M = A + B V + C V^2 + D V^3 + E V^4 + F V^5 + G V^6 + H V^7
    A = float(calibration_coefficients.get('Const.', 0.0))
    B = float(calibration_coefficients.get('Coef.2', 0.0))
    C = float(calibration_coefficients.get('Coef.3', 0.0))
    D = float(calibration_coefficients.get('Coef.4', 0.0))
    E = float(calibration_coefficients.get('Coef.5', 0.0))
    F = float(calibration_coefficients.get('Coef.6', 0.0))
    G = float(calibration_coefficients.get('Coef.7', 0.0))
    H = float(calibration_coefficients.get('Coef.8', 0.0))


    #Vmin and Vmax given by the PSS software   18.446150
    Vmin = df.iloc[0,0] #row.get('Vmin')
    Vmax = df.iloc[-1,0] #row.get('Vmax')
    ve_raw = np.clip(df.iloc[:, 0], Vmin, Vmax)
    hv_raw = df.iloc[:, 2] - np.min(df.iloc[:, 2])
    m_raw = df.iloc[:, 1]

    # Filter data where Vmin < ve < Vmax
    filtered_indices = np.where((ve_raw > Vmin) & (ve_raw < Vmax))[0]
    ve = np.array(ve_raw.iloc[filtered_indices])
    hv = hv_raw[filtered_indices]
    m = np.array(m_raw[filtered_indices])

    
    hv_max = np.max(hv) #heighest value of h(V)
    hv_height_norm = hv/hv_max #h(V) normalized by the heighest value of h(V) - range 0 to 1.0
  

    # Normalization of h(V) by area
    delta_ve = ve[1]-ve[0]
    area_hv_total = simpson(y=hv,x=ve,dx=delta_ve) # area calculated by the Simpson's rule (numerical calculation of areas and volumes)
    # possible to use *from numpy import traz // a =  trapz(y, dx)* as well where the calculation is done by the trapezoidal rule
    hv_area_norm = hv/area_hv_total #h(V) normalized by the area of h(V) curve


    logM = A + B * ve + C * (ve ** 2) + D * (ve ** 3) + E * (ve ** 4) + F * (ve ** 5) + G * (ve ** 6) + H * (ve ** 7) # applying calibration curve to the values
    # m = 10 ** logM #matrix of Molar Mass for every given elution volume (Ve)
    slope = B + 2 * C * ve + 3 * D * (ve ** 2) + 4 * E * (ve ** 3) + 5 * F * (ve ** 4) + 6 * G * (ve ** 5) + 7 * H * (ve ** 6) #slope of the calibration curve to every dV
    

    xm = (-1 / slope) * hv_area_norm #x(M)
    wm = np.log10(np.exp(1))*(xm/m) #w(M), logarithm of base 10 of 'e'
    wm_m = wm*m #w(M)M
    nm = wm/m #n(M)

    delta_m = m[1]- m[0]
    area_wm_total = simpson(y=wm,x=m,dx=delta_m)
    area_wm_m_total = simpson(y=wm_m,x=m,dx=delta_m)
    area_nm_total = simpson(y=nm,x=m,dx=delta_m)
    
    # Ensure m, xm, wm, and nm have the same length
    min_length = min(len(m), len(xm), len(wm), len(nm))
    m = m[:min_length]
    xm = xm[:min_length]
    wm = wm[:min_length]
    nm = nm[:min_length]


    if area_wm_total != 0 and area_nm_total != 0:
        total_mw = area_wm_m_total / area_wm_total # Mw
        total_mn = area_wm_total / area_nm_total # Mn
        total_pdi = total_mw / total_mn # PDI
        total_sd = np.sqrt(total_pdi - 1) * total_mn

        print('')
        print('The',file_name,'chromatogram presents the following values:')
        print('')
        print('Mw =', total_mw)
        print('Mn =', total_mn)
        print('Mw/Mn =', total_pdi)
        print('')
    else:
        print("Error: Division by zero or NaN encountered.")

    results.append({
        'filename': file_name,
        'Mw': total_mw,
        'Mn': total_mn,
        'PDI': total_pdi,
        'SD': total_sd,
    })

    generate_data = row.get('Process table? (y/n)', 'n')

    if isinstance(generate_data, str) and generate_data.lower() == 'y':
        # New file for the processed data
        data = {
            'Ve': ve,
            'h(V)': hv,
            'h(V) height norm': hv_height_norm,
            'h(V) area norm': hv_area_norm,
            'log M': logM,
            'M': m,
            'x(M)': xm,
            'w(M)': wm,
            'w(M)M': wm_m,
            'n(M)': nm
        }

        df_new = pd.DataFrame(data)

        file_name_wout_extension = file_name.replace('.TXT', '')
        
        # Write the DataFrame to a new Excel file
        output_file = file_name_wout_extension + "_processed.txt"
        df_new.to_csv(output_file, index=False)

        print(f"Calculated values saved to {output_file}")

    plot_data = row.get('Plot? (y/n)', 'n')

    if isinstance(plot_data, str) and plot_data.lower() == 'y':
        #SEC graphs
        plt.figure(figsize=(12, 6))
        #plt.subplot(row, column, order)
        plt.subplot(2, 2, 1)
        plt.plot(ve, hv_height_norm)
        plt.xlabel('Ve / mL'); plt.ylabel('h(V)')
        plt.subplot(2, 2, 2)
        plt.plot(m, xm)
        plt.xscale('log')
        plt.xlabel('M / g mol-1'); plt.ylabel('x(M)')
        plt.subplot(2, 2, 3)
        plt.plot(m, wm)
        plt.xscale('log')
        plt.xlabel('M / g mol-1'); plt.ylabel('w(M)')
        plt.subplot(2, 2, 4)
        plt.plot(m, nm)
        plt.xscale('log')
        plt.xlabel('M / g mol-1'); plt.ylabel('n(M)')
        plt.suptitle('SEC-GPC graphs of  ' + file_name)

    save_plot = row.get('Save Plot? (y/n)', 'n')

    if isinstance(save_plot, str) and save_plot.lower() == 'y':
        plt.savefig((file_name_wout_extension + '.jpg').format(os.path.basename(folder)),dpi=600)
plt.show()


results_df = pd.DataFrame(results)
report_filename = 'report.csv'
results_df.to_csv(report_filename, index=False)
print(f'Report saved as {report_filename}')