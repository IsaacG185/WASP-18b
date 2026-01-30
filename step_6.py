# Just edit this to show primary and secondary transits. Explain phase modulations (whatever that means)
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
plt.ylabel('Power (W)')
plt.title('BLS Periodogram')
plt.show()

# Step 5: Identify the strongest period
model = at.BoxLeastSquares(time, flux)
results = model.power(periods, durations)
index = np.argmax(results.power)

# Desired results:
# Period: 0.94 days
# Radius: 1.24 Jupiter radii
# Impact Parameter: ~0.32
best_period = bls_result.period[np.argmax(bls_result.power)]
print(f"Best-fit Period: {best_period:.4f} days")
t0 = results.transit_time[index]
duration = results.duration[index]

lc = lk.LightCurve(time=masked_time, flux=masked_flux)
folded_lc = lc.fold(period=best_period, epoch_time=t0)

# Shift the phase for three separate transits
phase = folded_lc.time.value
flux = folded_lc.flux.value

# Shift phases for three distinct transits
phase_shifted_pos1 = phase + best_period  # First shift
phase_shifted_pos2 = phase + 2 * best_period  # Second shift

# Concatenate the phases for three different transits
all_phases = np.concatenate((phase, phase_shifted_pos1, phase_shifted_pos2))
all_fluxes = np.concatenate((flux, flux, flux))

# Plot the three distinct transits
plt.scatter(all_phases, all_fluxes, s=5)
plt.xlabel('Phase (Radians)')
plt.ylabel('Normalized Flux (W/mˆ2')
plt.title('Three Separate Transits of WASP-18b')
plt.xlim(-1.5 * best_period, 3 * best_period)  # Adjust x-axis to fit three transits
plt.show()

# Scale the phase to stretch transits horizontally
scaled_phase = all_phases * 2.5  # Adjust the scale factor as needed

# Plot the horizontally stretched transits
plt.figure(figsize=(10, 6))
plt.scatter(scaled_phase, all_fluxes, s=10, alpha=0.5, label="Raw Data")

# Optional: Apply binning for clarity
bins = np.linspace(scaled_phase.min(), scaled_phase.max(), 500)  # Define bins
bin_indices = np.digitize(scaled_phase, bins)  # Assign phases to bins
binned_flux = [np.mean(all_fluxes[bin_indices == i]) for i in range(len(bins))]
plt.plot(bins, binned_flux, color='r', linewidth=1.5, label="Binned Data")

# Label and format the plot
plt.xlabel('Scaled Phase (Radians)')
plt.ylabel('Normalized Flux (W/mˆ2')
plt.title('(Horizontally Stretched) Transits of WASP-18b')
plt.xlim(-1, 5.8)  # Adjust x-axis range based on scaling
plt.legend()
plt.show()