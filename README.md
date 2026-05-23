# In vitro droplet / condensate analysis

This repository contains a Python script for quantifying condensates from microscopy `.nd2` files.

The script detects condensates, saves quality-control overlays, groups replicate files by condition, and plots condensate fraction and partition ratio.

---

# What the script does

For every `.nd2` file, the script:

1. Opens the image
2. Selects the chosen channel
3. Optionally performs max projection if the image is a Z-stack
4. Subtracts smooth background
5. Enhances small bright condensate-like objects
6. Detects condensates using local thresholding
7. Filters detected objects by:
   - size
   - circularity
   - solidity
   - intensity
8. Quantifies condensate metrics
9. Saves CSV tables and summary plots

The script automatically creates an `analysis` folder inside the selected input folder.

---

# File naming and replicate grouping

Replicates are grouped automatically by removing the final underscore number from the filename.

Example:

```
20260522-3.75uMBcat+159+0uM-1433_001.nd2
20260522-3.75uMBcat+159+0uM-1433_005.nd2
```

Both files are grouped as:

```
20260522-3.75uMBcat+159+0uM-1433
```

This allows the script to calculate condition averages and error bars automatically.

---

# Main outputs

Inside the `analysis` folder, the script saves:

```
per_file_summary.csv
per_condensate_measurements.csv
condition_summary.csv
condensate_fraction_by_condition.png
partition_ratio_by_condition.png
overlays/
```

---

# Quality control overlays

The `overlays` folder contains microscopy-style images showing:

- green raw fluorescence signal
- red outlines around detected condensates

Always inspect these overlays before interpreting the quantification results.

The overlays are the main quality-control step of the pipeline.

---

# Metrics calculated

## Condensate fraction

Condensate fraction measures how much of the total image signal is localized inside detected condensates.

Formula:

```
condensate fraction =total intensity inside condensates / total image intensity
```

This metric combines both:
- condensate size
- condensate brightness

A lower condensate fraction means less total signal is concentrated inside condensates.

This metric is relatively robust to differences in total image brightness or exposure.

---

## Partition ratio

Partition ratio measures how enriched the signal is inside condensates compared with the surrounding dilute phase.

Formula:

```
partition ratio = mean intensity inside condensates / mean intensity outside condensates
```

A high partition ratio means the signal is strongly enriched inside condensates.
A low partition ratio means the signal is more evenly distributed between condensates and the dilute phase.

---

# Difference between condensate fraction and partition ratio

Condensate fraction asks: ``` How much of the total signal is inside condensates?```

Partition ratio asks: ``` How concentrated is the signal inside condensates compared with outside? ```

Both metrics are useful and complementary.

Condensate fraction is usually more robust to imaging variability.

Partition ratio is often more sensitive to changes in condensate enrichment.

---

# Installation

Install required Python packages:

```bash
pip install nd2 numpy pandas scipy scikit-image matplotlib
```

---

# Running the analysis

Run the script:

```bash
python condensate_analysis.py
```

A folder-selection window will open.

Select the folder containing your `.nd2` files.

The results will automatically be saved in:

``` selected_folder/analysis/ ```

---

# Recommended workflow

1. Run the script
2. Open the `analysis/overlays` folder
3. Verify that red outlines match real condensates
4. Only then interpret the plots and CSV outputs
5. If detection is too strict or too permissive, tune the settings at the top of the script

---

# Important settings

The most important settings to tune are:

```python
threshold_method = "local"
local_block_size = 51
local_offset = -0.01

min_condensate_area_px = 3
max_condensate_area_px = 500

min_mean_intensity_factor = 1.20
```

If too few condensates are detected, lower:

```python
min_mean_intensity_factor
```

Example:

```python
min_mean_intensity_factor = 1.05
```

If too many background spots are detected, increase:

```python
min_mean_intensity_factor
```

Example:

```python
min_mean_intensity_factor = 1.50
```

---

# Notes

This script was designed for in vitro droplet / biomolecular condensate assays with bright puncta or droplets on a relatively dim background.

Different microscopes, fluorophores, acquisition settings, or sample types may require parameter tuning.

Always validate segmentation quality using the overlay outputs.
