# WASP-18 b Transit Analysis Pipeline

## Overview
This project implements a simplified exoplanet analysis pipeline using TESS photometric data for the hot Jupiter WASP-18 b. The pipeline follows the steps of a typical exoplanet discovery and characterization workflow, from downloading light curve data to measuring primary and secondary transits and constructing the full phase curve. All code is written in Python, and figures and results are summarized in a LaTeX paper.

## Code Structure

The code is organized according to the main steps of the analysis:

1. **Data Acquisition (Step 1)**
   - Downloads all 2-minute cadence TESS data for WASP-18 b from MAST.
   - Saves each sector's `.fits` file for subsequent processing.

2. **Data Preprocessing (Step 2)**
   - Reads the `TIME`, `PDCSAP_FLUX`, and `QUALITY` fields using `astropy.io.fits`.
   - Masks data flagged as low-quality.
   - Normalizes flux by its median and concatenates all sectors into a single light curve.
   - Plots and saves the concatenated light curve (`Figure 1.png`).

3. **Rediscovering WASP-18 b (Step 3)**
   - Runs the Box-Least Squares (BLS) periodogram on the processed light curve.
   - Identifies transit periods and depths.
   - Calculates the planetary radius and transit impact parameter using host star properties.
   - Saves the periodogram plot (`Figure 2.png`).

4. **Primary Transits (Step 4)**
   - Folds the light curve on the best-fit orbital period.
   - Plots the phase curve across three transit durations, highlighting the primary transit.
   - Saves the phase curve plot (`Figure 3.png`).

5. **Secondary Eclipse (Step 5)**
   - Masks primary transits and reruns BLS to search for secondary eclipses.
   - Measures the depth and phase of secondary eclipses.
   - Estimates the dayside temperature of WASP-18 b.
   - Saves the secondary eclipse plot (`Figure 4.png`).

6. **Full Phase Curve (Step 6)**
   - Generates the complete phase curve from BLS and phase folding.
   - Highlights primary and secondary transits.
   - Saves the full phase curve plot (`Figure 5.png`).

## Dependencies
- Python 3.x
- `numpy`
- `matplotlib`
- `astropy`

## Usage
1. Ensure TESS `.fits` files for WASP-18 b are downloaded (Step 1 handles this automatically).
2. Run the Python scripts in order (Steps 2–6) to generate all plots and perform calculations.
3. Compile the LaTeX paper (`WASP18b_Project.tex`) to produce the final report with figures and results.

## Figures
- **Figure 1:** Concatenated light curve of WASP-18 b.
- **Figure 2:** BLS periodogram showing the transit period.
- **Figure 3:** Phase curve highlighting primary transits.
- **Figure 4:** Phase curve highlighting secondary eclipses.
- **Figure 5:** Full phase curve with primary and secondary transits labeled.

## Notes
- All uncertainties in measurements are propagated from flux error estimates.
- Axis labels, units, and figure captions are included in the LaTeX paper.
- The code is modular, with each step corresponding to a major analysis stage.
