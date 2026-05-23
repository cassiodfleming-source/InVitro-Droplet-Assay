# ============================================================
# Generic in vitro droplet / condensate analysis from ND2 files
# ============================================================

############################
# USER SETTINGS
############################

# ND2 reading
channel_to_analyze = 0
use_max_projection = True
frame_to_analyze = 0

# Preprocessing
subtract_background = True
background_sigma = 50
smooth_sigma = 1.0

# Spot enhancement
use_spot_enhancement = True
spot_sigma_small = 1
spot_sigma_large = 4

# Segmentation
threshold_method = "local"      # "local", "otsu", or "yen"
threshold_multiplier = 0.75
local_block_size = 51
local_offset = -0.01

# Condensate filters
min_condensate_area_px = 3
max_condensate_area_px = 500
min_circularity = 0.05
min_solidity = 0.30
min_mean_intensity_factor = 1.20

# Analysis
exclude_image_border_px = 2
outside_erode_px = 1

# Optional outputs
save_overlay_images = True
save_individual_masks = False
save_plots_png = True
save_plots_pdf = True
run_statistics = True

# Plot settings
error_bar = "sem"
dpi = 300


############################
# IMPORTS
############################

from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog

from scipy import ndimage as ndi
from skimage import filters, measure, morphology, segmentation, exposure
from skimage.io import imsave

try:
    import nd2
except ImportError:
    raise ImportError("Install nd2 first with: pip install nd2")


############################
# FOLDER SELECTION
############################

def select_input_folder():
    root = tk.Tk()
    root.withdraw()

    folder = filedialog.askdirectory(
        title="Select folder containing ND2 files"
    )

    root.destroy()

    if not folder:
        raise ValueError("No folder selected.")

    return Path(folder)


############################
# HELPER FUNCTIONS
############################

def get_condition_name(file_path):
    """
    Removes final replicate number.
    Example:
    sample_condition_001.nd2 -> sample_condition
    sample_condition_005.nd2 -> sample_condition
    """

    return re.sub(r"_\d+$", "", file_path.stem)


def read_nd2_image(file_path):

    data = nd2.imread(str(file_path))
    arr = np.squeeze(np.asarray(data))

    if arr.ndim == 2:
        return arr.astype(float)

    if arr.ndim == 3:

        if arr.shape[0] <= 4:
            img = arr[channel_to_analyze]

        else:
            img = (
                np.max(arr, axis=0)
                if use_max_projection
                else arr[0]
            )

        return img.astype(float)

    while arr.ndim > 2:

        if arr.shape[0] <= 4:
            arr = arr[channel_to_analyze]

        else:
            arr = arr[
                frame_to_analyze
                if frame_to_analyze < arr.shape[0]
                else 0
            ]

    return arr.astype(float)


def preprocess_image(img):

    img = img.astype(float)

    if subtract_background:

        background = ndi.gaussian_filter(
            img,
            sigma=background_sigma
        )

        img = img - background
        img[img < 0] = 0

    if smooth_sigma > 0:

        img = ndi.gaussian_filter(
            img,
            sigma=smooth_sigma
        )

    return img


def enhance_spots(img):

    small = ndi.gaussian_filter(
        img,
        sigma=spot_sigma_small
    )

    large = ndi.gaussian_filter(
        img,
        sigma=spot_sigma_large
    )

    enhanced = small - large
    enhanced[enhanced < 0] = 0

    return enhanced


def segment_condensates(processed_img, raw_img):

    valid_mask = np.ones(
        processed_img.shape,
        dtype=bool
    )

    if exclude_image_border_px > 0:

        valid_mask[:exclude_image_border_px, :] = False
        valid_mask[-exclude_image_border_px:, :] = False
        valid_mask[:, :exclude_image_border_px] = False
        valid_mask[:, -exclude_image_border_px:] = False

    work_img = processed_img.copy()

    if use_spot_enhancement:
        work_img = enhance_spots(work_img)

    if threshold_method.lower() == "local":

        threshold = filters.threshold_local(
            work_img,
            block_size=local_block_size,
            offset=local_offset
        )

        mask = work_img > threshold

    elif threshold_method.lower() == "otsu":

        pixels = work_img[valid_mask]

        threshold = (
            filters.threshold_otsu(pixels)
            * threshold_multiplier
        )

        mask = work_img > threshold

    elif threshold_method.lower() == "yen":

        pixels = work_img[valid_mask]

        threshold = (
            filters.threshold_yen(pixels)
            * threshold_multiplier
        )

        mask = work_img > threshold

    else:
        raise ValueError(
            "threshold_method must be "
            "'local', 'otsu', or 'yen'."
        )

    mask &= valid_mask

    mask = morphology.remove_small_objects(
        mask,
        min_size=min_condensate_area_px
    )

    mask = ndi.binary_fill_holes(mask)

    labeled = measure.label(mask)

    clean_mask = np.zeros_like(
        mask,
        dtype=bool
    )

    background_intensity = np.median(
        raw_img[valid_mask]
    )

    for region in measure.regionprops(
        labeled,
        intensity_image=raw_img
    ):

        area = region.area
        perimeter = region.perimeter

        circularity = 0

        if perimeter > 0:
            circularity = (
                4 * np.pi * area
                / (perimeter ** 2)
            )

        keep = (
            area >= min_condensate_area_px and
            area <= max_condensate_area_px and
            circularity >= min_circularity and
            region.solidity >= min_solidity and
            region.mean_intensity >= (
                background_intensity
                * min_mean_intensity_factor
            )
        )

        if keep:
            clean_mask[
                labeled == region.label
            ] = True

    return clean_mask


def quantify_image(raw_img, condensate_mask):

    total_intensity = np.sum(raw_img)

    inside_pixels = raw_img[
        condensate_mask
    ]

    outside_mask = ~condensate_mask

    if outside_erode_px > 0:

        dilated = morphology.binary_dilation(
            condensate_mask,
            morphology.disk(outside_erode_px)
        )

        outside_mask = ~dilated

    outside_pixels = raw_img[
        outside_mask
    ]

    condensate_intensity = np.sum(
        inside_pixels
    )

    condensate_fraction = np.nan

    if total_intensity > 0:

        condensate_fraction = (
            condensate_intensity
            / total_intensity
        )

    mean_inside = (
        np.mean(inside_pixels)
        if inside_pixels.size > 0
        else np.nan
    )

    mean_outside = (
        np.mean(outside_pixels)
        if outside_pixels.size > 0
        else np.nan
    )

    partition_ratio = np.nan

    if mean_outside > 0:

        partition_ratio = (
            mean_inside
            / mean_outside
        )

    labeled = measure.label(
        condensate_mask
    )

    summary = {
        "n_condensates": int(np.max(labeled)),
        "total_image_intensity": total_intensity,
        "total_condensate_intensity": condensate_intensity,
        "condensate_fraction": condensate_fraction,
        "mean_inside_condensates": mean_inside,
        "mean_outside_condensates": mean_outside,
        "partition_ratio": partition_ratio,
        "condensate_area_px": int(
            np.sum(condensate_mask)
        ),
        "image_area_px": raw_img.size,
        "condensate_area_fraction": (
            np.sum(condensate_mask)
            / raw_img.size
        )
    }

    props = measure.regionprops_table(
        labeled,
        intensity_image=raw_img,
        properties=[
            "label",
            "area",
            "mean_intensity",
            "max_intensity",
            "centroid",
            "perimeter",
            "solidity"
        ]
    )

    condensate_table = pd.DataFrame(props)

    if len(condensate_table) > 0:

        condensate_table[
            "integrated_intensity"
        ] = (
            condensate_table["area"]
            * condensate_table["mean_intensity"]
        )

        condensate_table[
            "circularity"
        ] = (
            4 * np.pi
            * condensate_table["area"]
            / (condensate_table["perimeter"] ** 2)
        )

    return summary, condensate_table


def save_overlay(raw_img, mask, output_path):

    # --------------------------------------------------------
    # FIXED MICROSCOPY-LIKE CONTRAST
    # --------------------------------------------------------

    p_low = np.percentile(raw_img, 1)
    p_high = np.percentile(raw_img, 99.8)

    img_norm = exposure.rescale_intensity(
        raw_img,
        in_range=(p_low, p_high),
        out_range=(0, 1)
    )

    img_norm = np.clip(img_norm, 0, 1)

    # --------------------------------------------------------
    # CREATE GREEN MICROSCOPY LOOK
    # --------------------------------------------------------

    rgb = np.zeros(
        (
            img_norm.shape[0],
            img_norm.shape[1],
            3
        )
    )

    # Green LUT
    rgb[..., 1] = img_norm

    # Small white contribution
    rgb[..., 0] = img_norm * 0.15
    rgb[..., 2] = img_norm * 0.15

    # --------------------------------------------------------
    # RED CONDENSATE OUTLINES
    # --------------------------------------------------------

    boundaries = segmentation.find_boundaries(
        mask,
        mode="outer"
    )

    rgb[boundaries, 0] = 1
    rgb[boundaries, 1] = 0
    rgb[boundaries, 2] = 0

    # --------------------------------------------------------
    # SAVE
    # --------------------------------------------------------

    imsave(
        output_path,
        (rgb * 255).astype(np.uint8)
    )


def summarize_conditions(df):

    if error_bar == "sem":

        error_function = (
            lambda x:
            x.std(ddof=1)
            / np.sqrt(x.count())
        )

    elif error_bar == "sd":

        error_function = (
            lambda x:
            x.std(ddof=1)
        )

    else:
        raise ValueError(
            "error_bar must be "
            "'sem' or 'sd'."
        )

    return (
        df.groupby(
            "condition",
            as_index=False
        )
        .agg(
            n_files=("file", "nunique"),

            mean_condensate_fraction=(
                "condensate_fraction",
                "mean"
            ),

            error_condensate_fraction=(
                "condensate_fraction",
                error_function
            ),

            mean_partition_ratio=(
                "partition_ratio",
                "mean"
            ),

            error_partition_ratio=(
                "partition_ratio",
                error_function
            ),

            mean_n_condensates=(
                "n_condensates",
                "mean"
            )
        )
    )


def plot_bar(
    summary_df,
    value_col,
    error_col,
    ylabel,
    output_name,
    output_folder
):

    fig, ax = plt.subplots(
        figsize=(
            max(6, len(summary_df) * 0.8),
            5
        )
    )

    x = np.arange(len(summary_df))

    ax.bar(
        x,
        summary_df[value_col],
        yerr=summary_df[error_col],
        capsize=4
    )

    ax.set_xticks(x)

    ax.set_xticklabels(
        summary_df["condition"],
        rotation=45,
        ha="right"
    )

    ax.set_ylabel(ylabel)
    ax.set_xlabel("")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    if save_plots_png:

        fig.savefig(
            output_folder / f"{output_name}.png",
            dpi=dpi
        )

    if save_plots_pdf:

        fig.savefig(
            output_folder / f"{output_name}.pdf"
        )

    plt.close(fig)


############################
# MAIN
############################

def main():

    input_folder = select_input_folder()

    output_folder = (
        input_folder / "analysis"
    )

    output_folder.mkdir(
        exist_ok=True
    )

    overlay_folder = (
        output_folder / "overlays"
    )

    mask_folder = (
        output_folder / "masks"
    )

    if save_overlay_images:
        overlay_folder.mkdir(exist_ok=True)

    if save_individual_masks:
        mask_folder.mkdir(exist_ok=True)

    nd2_files = sorted(
        input_folder.glob("*.nd2")
    )

    if len(nd2_files) == 0:

        raise FileNotFoundError(
            f"No ND2 files found in: "
            f"{input_folder}"
        )

    all_summaries = []
    all_condensates = []

    for file_path in nd2_files:

        print(
            f"Analyzing: "
            f"{file_path.name}"
        )

        condition = get_condition_name(
            file_path
        )

        raw_img = read_nd2_image(
            file_path
        )

        processed_img = preprocess_image(
            raw_img
        )

        condensate_mask = segment_condensates(
            processed_img,
            raw_img
        )

        summary, condensate_table = quantify_image(
            raw_img,
            condensate_mask
        )

        summary["file"] = file_path.name
        summary["condition"] = condition

        all_summaries.append(summary)

        if len(condensate_table) > 0:

            condensate_table["file"] = (
                file_path.name
            )

            condensate_table["condition"] = (
                condition
            )

            all_condensates.append(
                condensate_table
            )

        if save_overlay_images:

            save_overlay(
                raw_img,
                condensate_mask,
                overlay_folder /
                f"{file_path.stem}_overlay.png"
            )

        if save_individual_masks:

            imsave(
                mask_folder /
                f"{file_path.stem}_mask.png",

                (
                    condensate_mask * 255
                ).astype(np.uint8)
            )

    per_file_summary = pd.DataFrame(
        all_summaries
    )

    per_file_summary = (
        per_file_summary
        .sort_values(
            ["condition", "file"]
        )
    )

    per_file_summary.to_csv(
        output_folder /
        "per_file_summary.csv",
        index=False
    )

    if len(all_condensates) > 0:

        per_condensate_measurements = (
            pd.concat(
                all_condensates,
                ignore_index=True
            )
        )

        per_condensate_measurements.to_csv(
            output_folder /
            "per_condensate_measurements.csv",
            index=False
        )

    condition_summary = summarize_conditions(
        per_file_summary
    )

    condition_summary.to_csv(
        output_folder /
        "condition_summary.csv",
        index=False
    )

    plot_bar(
        condition_summary,
        value_col="mean_condensate_fraction",
        error_col="error_condensate_fraction",
        ylabel="Condensate fraction",
        output_name="condensate_fraction_by_condition",
        output_folder=output_folder
    )

    plot_bar(
        condition_summary,
        value_col="mean_partition_ratio",
        error_col="error_partition_ratio",
        ylabel="Partition ratio",
        output_name="partition_ratio_by_condition",
        output_folder=output_folder
    )

    print("\nDone.")
    print(f"Input folder: {input_folder}")
    print(f"Results saved in: {output_folder}")


if __name__ == "__main__":
    main()