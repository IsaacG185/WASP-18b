import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

# Define directory and list of FITS files
dir_path = "/Users/isaacgutierrez/Desktop/Ampersand/Class Project/WASP-18b"
files = [f for f in os.listdir(dir_path) if f.endswith('.fits')]

flux, time = np.array([]), np.array([])

for lc in files:
    with fits.open(os.path.join(dir_path, lc)) as LC:
        current_flux = LC[1].data['PDCSAP_FLUX']
        current_time = LC[1].data['TIME']
        current_quality = LC[1].data['QUALITY']
        
        # Masking bad data
        valid_indices = (~np.isnan(current_flux)) & (np.bitwise_and(current_quality, 0b0101001010111111) == 0)
        masked_flux = current_flux[valid_indices]
        masked_time = current_time[valid_indices]
        
        if masked_flux.size > 0:  # Check for valid flux data
            masked_flux /= np.median(masked_flux)
            
            # Concatenate data
            flux = np.concatenate((flux, masked_flux))
            time = np.concatenate((time, masked_time))

# Plotting
plt.scatter(time, flux, marker='.', alpha=0.9)
plt.xlabel('Time (s)')
plt.ylabel('Flux (W/mˆ2)')
plt.title('Light Curve of WASP-18b (Concatenated)')
plt.show()