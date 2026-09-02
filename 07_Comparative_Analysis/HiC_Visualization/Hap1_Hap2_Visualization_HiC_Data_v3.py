import os
import cooler
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize, LinearSegmentedColormap
import matplotlib.ticker as ticker

# EXACT JUICEBOX COLOR PALETTE (White -> Light Red -> Deep Cherry Red)
JUICEBOX_CMAP = LinearSegmentedColormap.from_list(
    "juicebox", ["#ffffff", "#ff3333", "#a60000"]
)

def load_matrix_and_bounds(mcool_path, resolution=500000, balance=True):
    """Loads cooler matrix and calculates cumulative chromosome boundaries in Megabases (Mb)."""
    cooler_uri = f"{mcool_path}::/resolutions/{resolution}"
    clr = cooler.Cooler(cooler_uri)
    
    # Check for balancing weights (ICE normalization)
    is_balanced = 'weight' in clr.bins()
    if balance and is_balanced:
        print(f"Loading balanced matrix (ICE) for {os.path.basename(mcool_path)}...")
        matrix = clr.matrix(balance=True)[:]
        norm_type = "Normalized (ICE)"
    else:
        print(f"Loading raw matrix for {os.path.basename(mcool_path)}...")
        matrix = clr.matrix(balance=False)[:]
        norm_type = "Raw Counts"

    chromnames = clr.chromnames
    chromsizes = clr.chromsizes.values
    
    # Calculate chromosome boundary positions in Megabases (Mb)
    chrom_bounds_bp = np.concatenate([[0], np.cumsum(chromsizes)])
    chrom_bounds_mb = chrom_bounds_bp / 1e6
    
    return clr, matrix, chromnames, chrom_bounds_mb, norm_type

def plot_publication_hic_linear(clr, matrix, chromnames, chrom_bounds_mb, norm_type, output_dir, file_prefix, figure_title, vmax=None, mb_step=50):
    print(f"--- Generating Standardized Linear Plot for: {figure_title} ---")
    
    total_genome_mb = chrom_bounds_mb[-1]
    
    # Setup high-DPI figure canvas
    fig, ax = plt.subplots(figsize=(10, 10), dpi=300)
    
    # Replace NaN values with 0.0 for clean white background
    matrix_clean = np.nan_to_num(matrix, copy=True, nan=0.0)
    
    # Set Linear Intensity Cap
    if vmax is None:
        vmax = np.nanpercentile(matrix_clean[matrix_clean > 0], 99.5)
        
    norm = Normalize(vmin=0, vmax=vmax)

    # Render matrix with physical extent mapped directly to Megabases (Mb)
    img = ax.imshow(
        matrix_clean,
        cmap=JUICEBOX_CMAP,
        norm=norm,
        extent=[0, total_genome_mb, total_genome_mb, 0],
        interpolation="nearest"
    )
    
    # Draw chromosome boundary grid lines
    for bound in chrom_bounds_mb[1:-1]:
        ax.axhline(bound, color="black", linestyle="-", linewidth=0.5, alpha=0.8)
        ax.axvline(bound, color="black", linestyle="-", linewidth=0.5, alpha=0.8)

    # 1. Primary Axes (Bottom & Left): Centered Chromosome Labels
    chrom_centers = [(chrom_bounds_mb[i] + chrom_bounds_mb[i+1]) / 2 for i in range(len(chromnames))]
    ax.set_xticks(chrom_centers)
    ax.set_xticklabels(chromnames, rotation=45, ha='right', fontsize=12, fontweight='bold')
    ax.set_yticks(chrom_centers)
    ax.set_yticklabels(chromnames, fontsize=12, fontweight='bold')

    # 2. Secondary Axis (Top): Fixed 50 Mb Tick Steps across both plots
    secax_x = ax.secondary_xaxis('top')
    secax_x.set_xlabel('Genomic Distance (Mb)', fontsize=12, fontweight='bold', labelpad=10)
    secax_x.xaxis.set_major_locator(ticker.MultipleLocator(mb_step))
    secax_x.tick_params(labelsize=10)

    # Title & Colorbar
    ax.set_title(figure_title, fontsize=16, fontweight='bold', pad=30)
    
    cbar = fig.colorbar(img, ax=ax, shrink=0.8, pad=0.04)
    cbar.set_label(f"Contact Frequency ({norm_type}, Linear Scale)", fontsize=10)
    
    # Standardize colorbar ticks to neat steps
    cbar.locator = ticker.MaxNLocator(nbins=6)
    cbar.update_ticks()

    plt.tight_layout()
    
    # Save outputs in both SVG and PNG formats
    os.makedirs(output_dir, exist_ok=True)
    svg_out = os.path.join(output_dir, f"{file_prefix}.svg")
    png_out = os.path.join(output_dir, f"{file_prefix}.png")
    
    plt.savefig(svg_out, format="svg", bbox_inches="tight")
    print(f"Saved SVG to: {svg_out}")
    
    plt.savefig(png_out, format="png", bbox_inches="tight", dpi=300)
    print(f"Saved PNG to: {png_out}\n")
    
    plt.close()

def main():
    hap1_path = '/path/to/input/directory Test/Hap1_verification.mcool'
    hap2_path = '/path/to/input/directory Test/Hap2_verification.mcool'
    output_dir = '/path/to/output/directory Test/HiC_publication_figures'
    resolution = 500000  # 500 kb resolution
    
    # Load matrices
    clr1, mat1, chroms1, bounds1, norm1 = load_matrix_and_bounds(hap1_path, resolution, balance=True)
    clr2, mat2, chroms2, bounds2, norm2 = load_matrix_and_bounds(hap2_path, resolution, balance=True)

    # Calculate unified vmax across BOTH haplotypes, rounded up to nearest 100 for clean ticks
    vmax_hap1 = np.nanpercentile(mat1[mat1 > 0], 99.5) if np.any(mat1 > 0) else 1.0
    vmax_hap2 = np.nanpercentile(mat2[mat2 > 0], 99.5) if np.any(mat2 > 0) else 1.0
    raw_vmax = max(vmax_hap1, vmax_hap2)
    shared_vmax = float(np.ceil(raw_vmax / 100.0) * 100.0)

    # Plot Hap1 with fixed 50 Mb tick step
    plot_publication_hic_linear(
        clr=clr1, matrix=mat1, chromnames=chroms1, chrom_bounds_mb=bounds1, 
        norm_type=norm1, output_dir=output_dir, file_prefix="Hap1_publication_linear", 
        figure_title="HyPR01 Haplotype 1", vmax=shared_vmax, mb_step=50
    )

    # Plot Hap2 with fixed 50 Mb tick step
    plot_publication_hic_linear(
        clr=clr2, matrix=mat2, chromnames=chroms2, chrom_bounds_mb=bounds2, 
        norm_type=norm2, output_dir=output_dir, file_prefix="Hap2_publication_linear", 
        figure_title="HyPR01 Haplotype 2", vmax=shared_vmax, mb_step=50
    )

if __name__ == "__main__":
    main()