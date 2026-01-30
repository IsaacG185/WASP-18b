import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.timeseries as at

# Step 1: Reload the light curve
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

# Step 2: Plot the light curve
plt.figure(figsize=(10, 6))
plt.scatter(time, flux, marker='.', alpha=0.7)
plt.xlabel('Time (days)')
plt.ylabel('Normalized Flux')
plt.title('Light Curve of WASP-18b (Concatenated)')
plt.show()

# Step 3: Use BLS to find periodic transits
# Define the time and flux for BLS
time -= np.min(time)  # Normalize time
flux_err = np.std(flux)  # Estimate flux error as standard deviation

# Initialize the BLS model
bls = at.BoxLeastSquares(time, flux)

# Define the period range for search
periods = np.linspace(0.5, 10, 10000)
bls_result = bls.power(periods, duration=0.1)

# Step 4: Plot the periodogram
plt.figure(figsize=(10, 6))
plt.plot(bls_result.period, bls_result.power, color='b')
plt.xlabel('Period (days)')
plt.ylabel('Power (W)')
plt.title('BLS Periodogram of WASP-18b Light Curve (Concatenated)')
plt.show()

# Step 5: Identify the strongest period
best_period = bls_result.period[np.argmax(bls_result.power)]
print(f'(Best-Fit) Period: {best_period} days')

# Step 6: Fold the light curve at the best period
folded_time = (time % best_period) / best_period  # Phase
sorted_indices = np.argsort(folded_time)
folded_time = folded_time[sorted_indices]
folded_flux = flux[sorted_indices]

# Step 7: Calculate the depth of the transit
# Bin the data to reduce noise and isolate the transit
bin_width = 0.01  # Phase width for binning
bins = np.arange(0, 1 + bin_width, bin_width)
binned_flux, bin_edges = np.histogram(folded_time, bins=bins, weights=folded_flux)
binned_counts, _ = np.histogram(folded_time, bins=bins)
binned_flux /= binned_counts  # Average flux per bin

# Find the minimum flux in the binned light curve
transit_depth = 1 - np.min(binned_flux)
print(f"Transit Depth: {transit_depth:.6f}")

# Step 8: Calculate the planet radius
Rs = 1.319  # Stellar radius in solar radii
Rp = Rs * (transit_depth ** 0.5)
print(f"Radius of Planet: {Rp:.2f} Solar Radii")

a = 0.02024 * 215.032
i = 1.45735
b = np.cos(i) * (a / Rs)
print(f'Impact Parameter: {b}')
#Explain math + physics in Latex doc.
    #Desired results:
    #Period: 0.94 days
    #Radius: 0.11 Solar radii
    #Impact Parameter: ~0.37