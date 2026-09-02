import os
import sys
import glob
import re
import matplotlib.pyplot as plt
import numpy as np

def generate_db_plots(base_dir):
    db_keys = ["eukaryota", "embryophyta", "eudicots", "malpighiales"]
    run_types = ["Subgenome_A", "Subgenome_B", "Subgenome_C", "Subgenome_D"]
    
    # Publication-quality color scheme
    c_single = '#34b3e4'  # Light Blue
    c_dup = '#047495'     # Deep Slate Blue
    c_frag = '#f7d838'    # Warm Yellow
    c_miss = '#e93c26'    # Crimson Red

    for db in db_keys:
        print(f"Aggregating data for database: {db}")
        db_dir = os.path.join(base_dir, db)
        if not os.path.isdir(db_dir):
            print(f"Skipping database {db}: Folder not found.")
            continue
            
        data = {}
        for r_type in run_types:
            run_folder = f"{r_type}_{db}"
            search_path = os.path.join(db_dir, run_folder, "short_summary.specific.*.txt")
            summary_files = glob.glob(search_path)
            
            if not summary_files:
                print(f"  Warning: No summary files detected in {run_folder}")
                continue
                
            filepath = summary_files[0]
            
            with open(filepath, 'r') as f:
                content = f.read()
                
            try:
                s = int(re.search(r'\s+(\d+)\s+Complete and single-copy', content).group(1))
                d = int(re.search(r'\s+(\d+)\s+Complete and duplicated', content).group(1))
                frag = int(re.search(r'\s+(\d+)\s+Fragmented', content).group(1))
                miss = int(re.search(r'\s+(\d+)\s+Missing', content).group(1))
                total = int(re.search(r'\s+(\d+)\s+Total BUSCO groups searched', content).group(1))
                
                # Format label (e.g., Subgenome_A -> Subgenome A)
                pretty_label = r_type.replace("_", " ")
                
                data[pretty_label] = {
                    'Single': (s / total) * 100,
                    'Duplicated': (d / total) * 100,
                    'Fragmented': (frag / total) * 100,
                    'Missing': (miss / total) * 100
                }
            except AttributeError as e:
                print(f"  Error: Failed to regex-parse BUSCO values from {filepath}: {e}")
                continue
        
        if not data:
            print(f"  Skipping plot generation for {db} - empty data.")
            continue
            
        # Parse datasets
        labels = list(data.keys())
        single = [data[k]['Single'] for k in labels]
        dup = [data[k]['Duplicated'] for k in labels]
        frag = [data[k]['Fragmented'] for k in labels]
        miss = [data[k]['Missing'] for k in labels]
        
        fig, ax = plt.subplots(figsize=(11, 4.5))
        y_pos = np.arange(len(labels))
        
        # Plot bars
        ax.barh(y_pos, single, color=c_single, edgecolor='none', height=0.55, label='Complete (single-copy)')
        ax.barh(y_pos, dup, left=single, color=c_dup, edgecolor='none', height=0.55, label='Complete (duplicated)')
        ax.barh(y_pos, frag, left=np.array(single)+np.array(dup), color=c_frag, edgecolor='none', height=0.55, label='Fragmented')
        ax.barh(y_pos, miss, left=np.array(single)+np.array(dup)+np.array(frag), color=c_miss, edgecolor='none', height=0.55, label='Missing')
        
        # Plot Styling
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=11, fontweight='bold')
        ax.set_xlabel('% BUSCOs', fontsize=12, fontweight='semibold')
        ax.set_title(f'Tetraploid Subgenome BUSCO Comparison - {db.upper()} Lineage', fontsize=13, fontweight='bold', pad=15)
        ax.set_xlim(0, 100)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#cccccc')
        ax.spines['bottom'].set_color('#cccccc')
        ax.xaxis.grid(True, linestyle='--', alpha=0.5, color='#cccccc')
        ax.set_axisbelow(True)
        
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.18), ncol=4, frameon=False, fontsize=9)
        
        plt.tight_layout()
        out_plot_path = os.path.join(base_dir, f"{db}_subgenome_busco_comparison.png")
        plt.savefig(out_plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Successfully saved comparison plot to: {out_plot_path}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python generate_final_plots.py <base_dir>")
        sys.exit(1)
    generate_db_plots(sys.argv[1])
