#!/usr/bin/env python3
import sys
import os
import datetime
from collections import defaultdict

def parse_attributes(attr_string):
    attr = {}
    for item in attr_string.strip().split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            attr[k.strip()] = v.strip()
    return attr

class Logger:
    def __init__(self, log_path):
        self.terminal = sys.stdout
        self.log = open(log_path, "w")
        
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()
        
    def flush(self):
        self.terminal.flush()
        self.log.flush()

def process_single_haplotype(ref_gff, target_gff, output_gff):
    print(f"\n========================================================================")
    print(f" PROCESSING: {os.path.basename(target_gff)}")
    print(f"========================================================================")
    
    print("[-] Step 1: Parsing reference GFF3 to index original confidence tiers...")
    ref_id_to_tier = {}
    ref_child_to_parent = {}
    ref_mrna_to_cds = defaultdict(list)
    ref_mrna_info = {}

    with open(ref_gff, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            chrom, _, feature_type, start, end, _, strand, _, attributes = parts
            attrs = parse_attributes(attributes)
            fid = attrs.get('ID')
            parent = attrs.get('Parent')
            tier = attrs.get('confidence_tier')
            
            if fid and tier:
                ref_id_to_tier[fid] = tier
            if fid and parent:
                ref_child_to_parent[fid] = parent
            if feature_type == 'CDS' and parent:
                ref_mrna_to_cds[parent].append((int(start), int(end)))
                ref_mrna_info[parent] = (chrom, strand)

    exact_footprints = {}
    overlap_list = defaultdict(list)

    for mrna_id, cds_list in ref_mrna_to_cds.items():
        if not cds_list:
            continue
        chrom, strand = ref_mrna_info[mrna_id]
        sorted_cds = tuple(sorted(cds_list))
        
        tier = ref_id_to_tier.get(mrna_id)
        if not tier:
            parent_gene = ref_child_to_parent.get(mrna_id)
            if parent_gene:
                tier = ref_id_to_tier.get(parent_gene)
        
        if tier:
            exact_footprints[(chrom, strand, sorted_cds)] = tier
            overlap_list[chrom].append((sorted_cds[0][0], sorted_cds[-1][1], strand, tier, mrna_id))

    print("[-] Step 2: Extracting structural elements from UTR-updated target GFF3...")
    target_mrna_to_cds = defaultdict(list)
    target_mrna_info = {}
    target_child_to_parent = {}
    target_mrna_to_gene = {}

    with open(target_gff, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            chrom, _, feature_type, start, end, _, strand, _, attributes = parts
            attrs = parse_attributes(attributes)
            fid = attrs.get('ID')
            parent = attrs.get('Parent')
            
            if feature_type == 'mRNA' and fid and parent:
                target_mrna_to_gene[fid] = parent
            if fid and parent:
                target_child_to_parent[fid] = parent
            if feature_type == 'CDS' and parent:
                target_mrna_to_cds[parent].append((int(start), int(end)))
                target_mrna_info[parent] = (chrom, strand)

    print("[-] Step 3: Mapping identifiers using invariant CDS footprints...")
    target_id_to_tier = {}
    exact_matches = 0
    overlap_matches = 0
    failed_matches = 0

    for mrna_id, cds_list in target_mrna_to_cds.items():
        if not cds_list:
            continue
        chrom, strand = target_mrna_info[mrna_id]
        sorted_cds = tuple(sorted(cds_list))
        
        tier = exact_footprints.get((chrom, strand, sorted_cds))
        if tier:
            exact_matches += 1
        else:
            min_c = sorted_cds[0][0]
            max_c = sorted_cds[-1][1]
            best_overlap = 0
            best_tier = None
            matched_ref_id = "None"
            
            for o_start, o_end, o_strand, o_tier, o_id in overlap_list.get(chrom, []):
                if strand == o_strand:
                    overlap_len = max(0, min(max_c, o_end) - max(min_c, o_start))
                    if overlap_len > best_overlap:
                        best_overlap = overlap_len
                        best_tier = o_tier
                        matched_ref_id = o_id
            if best_overlap > 0:
                tier = best_tier
                overlap_matches += 1
                print(f"    [FALLBACK] Model {mrna_id} matched via overlap to Ref {matched_ref_id} ({best_overlap} bp shared)")
            else:
                failed_matches += 1
                print(f"    [WARNING] Model {mrna_id} could not be resolved to any reference tier.")
        
        if tier:
            target_id_to_tier[mrna_id] = tier
            gene_id = target_mrna_to_gene.get(mrna_id)
            if gene_id:
                target_id_to_tier[gene_id] = tier

    print(f"    -> Summary: {exact_matches} exact matches | {overlap_matches} fallback overlaps | {failed_matches} unresolved.")

    print("[-] Step 4: Writing tier-restored structural GFF3 file...")
    with open(target_gff, 'r') as infile, open(output_gff, 'w') as outfile:
        for line in infile:
            if line.startswith('#') or not line.strip():
                outfile.write(line)
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                outfile.write(line)
                continue
            
            attributes = parts[8]
            attrs = parse_attributes(attributes)
            
            assigned_tier = None
            fid = attrs.get('ID')
            if fid and fid in target_id_to_tier:
                assigned_tier = target_id_to_tier[fid]
            else:
                curr_parent = attrs.get('Parent')
                while curr_parent:
                    if curr_parent in target_id_to_tier:
                        assigned_tier = target_id_to_tier[curr_parent]
                        break
                    curr_parent = target_child_to_parent.get(curr_parent)
            
            if assigned_tier and 'confidence_tier' not in attrs:
                new_line = line.strip() + f";confidence_tier={assigned_tier}\n"
                outfile.write(new_line)
            else:
                outfile.write(line)
                
    print(f"[+] Output written cleanly to: {output_gff}")

if __name__ == "__main__":
    work_dir = "/path/to/your/directory"
    log_file = os.path.join(work_dir, "pipeline_restoration.log")
    
    sys.stdout = Logger(log_file)
    
    print(f"[*] Pipeline run initialized on {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"[*] Writing detailed executions to log file: {log_file}")
    
    haplotype_runs = [
        {
            "ref": os.path.join(work_dir, "harmonized_consensus_Hap1.gff3"),
            "target": os.path.join(work_dir, "HyPR01_Hap1_cleaned.gff3"),
            "output": os.path.join(work_dir, "HyPR01_Hap1_ordered.gff3") # Adjusted
        },
        {
            "ref": os.path.join(work_dir, "harmonized_consensus_Hap2.gff3"),
            "target": os.path.join(work_dir, "HyPR01_Hap2_cleaned.gff3"),
            "output": os.path.join(work_dir, "HyPR01_Hap2_ordered.gff3") # Adjusted
        }
    ]

    for run in haplotype_runs:
        if os.path.exists(run["ref"]) and os.path.exists(run["target"]):
            process_single_haplotype(run["ref"], run["target"], run["output"])
        else:
            print(f"[ERROR] Missing files for run setup. Check paths:\nRef: {run['ref']}\nTarget: {run['target']}")
            
    print(f"\n[*] Pipeline finished successfully on {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.")