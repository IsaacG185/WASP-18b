import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.timeseries as at
import lightkurve as lk

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
plt.scatter(time, flux, marker='.', alpha=0.5)
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
periods = np.linspace(0.8, 1.2, 10000)

# Adjust the duration for each period to avoid the error
# Transit duration is 5% of period, ensure it's always smaller than the period
durations = np.minimum(periods * 0.05, periods * 0.9)

bls_result = bls.power(periods, duration=durations)  # Use dynamic durations

# Step 4: Plot the periodogram
plt.figure(figsize=(10, 6))
plt.plot(bls_result.period, bls_result.power, color='b')
plt.xlabel('Period (days)')
plt.ylabel('Power')
plt.title('BLS Periodogram')
plt.show()

# Step 5: Identify the strongest period
model = at.BoxLeastSquares(time, flux)
results = model.power(periods, durations)
index = np.argmax(results.power)

# Get the best period, transit time, and duration
best_period = bls_result.period[np.argmax(bls_result.power)]
print(f"Best-fit Period: {best_period:.4f} days")

t0 = results.transit_time[index]
duration = results.duration[index]

lc = lk.LightCurve(time=masked_time, flux=masked_flux)
folded_lc = lc.fold(period=best_period, epoch_time=t0)

# Shift the phase for three separate transits
phase = folded_lc.time.value
flux = folded_lc.flux.value

# Step 1: Identify the regions of the two primary transits and mask them
# The first primary transit is around phase 0
# The second primary transit is around phase best_period, so around phase ~ 0.94
mask_primary_transits = ((phase > -0.05) & (phase < 0.05)) | ((phase > best_period - 0.05) & (phase < best_period + 0.05))

# Step 2: Exclude the flux values corresponding to the primary transits
phase_no_primary = phase[~mask_primary_transits]
flux_no_primary = flux[~mask_primary_transits]

# Step 3: Plot only the secondary transits (after excluding the two primary transits)
plt.scatter(phase_no_primary, flux_no_primary, s=5, color='orange', alpha=0.7)
plt.xlabel('Phase')
plt.ylabel('Normalized Flux')
plt.title('Secondary Transits of WASP-18b (without primary transits)')
plt.xlim(-1.5 * best_period, 3 * best_period)  # Adjust x-axis to fit the remaining transits
plt.show()

# Step 4: Optionally, adjust the scaling for the secondary transits
scaled_phase_no_primary = phase_no_primary * 4  # Adjust the scale factor as needed

# Plot the horizontally stretched secondary transits
plt.figure(figsize=(10, 6))
plt.scatter(scaled_phase_no_primary, flux_no_primary, s=10, alpha=0.5, label="Raw Data (Secondary Transits)")

# Optional: Apply binning for clarity
bins = np.linspace(scaled_phase_no_primary.min(), scaled_phase_no_primary.max(), 500)  # Define bins
bin_indices = np.digitize(scaled_phase_no_primary, bins)  # Assign phases to bins
binned_flux = [np.mean(flux_no_primary[bin_indices == i]) for i in range(len(bins))]
plt.plot(bins, binned_flux, color='r', linewidth=1.5, label="Binned Data")

# Label and format the plot
plt.xlabel('Scaled Phase (Radians)')
plt.ylabel('Normalized Flux (W/mˆ2)')
plt.title('Horizontally Stretched Secondary Transits of WASP-18b')
plt.xlim(-1.5 * best_period * 1, 2)  # Adjust x-axis range based on scaling
plt.legend()
plt.show()

# Step 1: Identify the flux values during the secondary transit
# We've already excluded the primary transits and are now working with secondary transits
# We assume the secondary transit flux corresponds to flux values around phase ~0.5, phase ~1.5, etc.
# Let's define the phase windows for the secondary transits. We expect the secondary transit around phase ~0.5.
# Adjust based on your results, if needed.
# Define phase window for secondary transit (e.g., phase ~0.5 ± 0.05)
mask_secondary_transits = ((phase_no_primary > 0.45) & (phase_no_primary < 0.55)) | \
                         ((phase_no_primary > 1.45) & (phase_no_primary < 1.55))

# Extract flux values during the secondary transit
secondary_transit_flux = flux_no_primary[mask_secondary_transits]

# Step 2: Calculate the out-of-transit baseline flux
# Baseline flux is calculated from the data that doesn't fall within the secondary transit phase windows
mask_out_of_transit = ~mask_secondary_transits
out_of_transit_flux = flux_no_primary[mask_out_of_transit]

# Step 3: Calculate the baseline flux (average of the out-of-transit flux)
baseline_flux = np.mean(out_of_transit_flux)

# Step 4: Calculate the mean flux during the secondary transit
secondary_transit_mean_flux = np.mean(secondary_transit_flux)

# Step 5: Calculate the secondary transit depth
secondary_depth = 1 - (secondary_transit_mean_flux / baseline_flux)

# Output the secondary transit depth
print(f"Secondary Transit Depth: {secondary_depth:.4f}")

# Given values
R_sun = 6.96e8  # Solar radius in meters
R_jupiter = 7.1492e7  # Jupiter radius in meters
R_star = 1.378 * R_sun  # Star radius in meters (WASP-18)
R_planet = 1.2 * R_jupiter  # Planet radius in meters (WASP-18b)
T_star = 6400  # Star temperature in K (WASP-18)
delta_F = 0.0011  # Secondary eclipse depth

# Calculate the dayside temperature of the planet
T_planet = T_star * (R_star / R_planet) * (delta_F)**(6/13)
print(f'Temperature of WASP-18b (Day-Side) (K): {T_planet}')