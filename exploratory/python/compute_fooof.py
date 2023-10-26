import os
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt

# Import the FOOOF object
from fooof import FOOOF

# Data path
current_file_path = os.path.abspath(os.path.dirname(__file__))
data_path = os.path.abspath(os.path.join(current_file_path,'../../../results/features/power/PFC'))
# Output paths
out_path = os.path.abspath(os.path.join(current_file_path,'../../../results/fooof/'))
figures_path = os.path.abspath(os.path.join(current_file_path,'../../../results/fooof/figures'))

# Frequency range to fit the fooof model
freq_range = [2, 40]

for file in os.listdir(data_path):
     filename = os.fsdecode(file)
     bidsID = os.path.splitext(filename)[0]

     # Import power files from matlab
     mat = sio.loadmat(os.path.join(data_path,filename))
     # Extract power and frequency from matlab structure
     pow = np.squeeze(mat['pow'])
     freq = np.squeeze(mat['freq'])
     # Average across all PFC ROIS
     avgpow = np.average(pow,0)

     # Initialize a FOOOF object
     fm = FOOOF()
     # Fit the model
     fm.fit(freq, avgpow, freq_range)
     # Save report
     fm.save_report(bidsID, file_path=out_path)
     # Extract FOOOF results from object
     fooof_results = fm.get_results()
     # Convert FOOOF results to a dictionary
     fooof_results_dict = fooof_results._asdict()
     # Save FOOOF results out to a mat file
     fnameout = os.path.join(out_path, bidsID + '_fooof.mat')
     sio.savemat(fnameout, fooof_results_dict)

     # Plotting
     fig_name = bidsID + '_foof.svg'
     fig, ax = plt.subplots(figsize=[7, 5])
     # plt.Figure(figsize = [7,5])
     # ax = plt.axis()
     plt.title(bidsID);
     fm.plot(save_fig=True,file_name=fig_name,file_path=figures_path, ax=ax)