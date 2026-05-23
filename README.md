# In vitro droplet / condensate analysis

This repository contains a Python script for quantifying condensates from microscopy `.nd2` files.

The script detects condensates, saves quality-control overlays, groups replicate files by condition, and plots condensate fraction and partition ratio.

## What the script does

For every `.nd2` file, the script:

1. Opens the image.
2. Selects the chosen channel.
3. Optionally performs max projection if the image is a Z-stack.
4. Subtracts smooth background.
5. Enhances small bright condensate-like objects.
6. Detects condensates using local thresholding.
7. Filters detected objects by size, circularity, solidity, and intensity.
8. Quantifies condensate metrics.
9. Saves CSV tables and summary plots.

The script automatically creates an `analysis` folder inside the selected input folder.

## File naming and replicate grouping

Replicates are grouped by removing the final underscore number from the filename.

Example:

```text
20260522-3.75uMBcat+159+0uM-1433_001.nd2
20260522-3.75uMBcat+159+0uM-1433_005.nd2
