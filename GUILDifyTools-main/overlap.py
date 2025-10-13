#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
Finds genes in the top 5% across all phenotypes **INCLUDING Astrocytes Underexpressed**.
Labels genes as "seed" if their **Gene Symbol** appears in `seeds.txt` of **ANY disease**, including Astrocytes.

Output:
  - `top5_overlap_results.txt` with **correct Seed Status & Source Disease**.
"""

import os
import sys
import csv

# **Include Astrocytes in Disease Directories for Seed Checking**
DISEASE_DIRS = {
    "Alzheimer": "astrocytes/alzheimer/",
    "Cardiovascular": "astrocytes/cardiovascular/",
    "Diabetes": "astrocytes/diabetes/",
    "Inflammation": "astrocytes/inflammation/",
    "Oxidative_Stress": "astrocytes/oxidative_stress/",
    "Cognitive_Impairment": "astrocytes/cognitive_impairment/",
    "Astrocytes": "astrocytes/astrocytes_overexpressed"  # ✅ Astrocytes now included
}

### **1️⃣ Load Disease Seeds (Including Astrocytes)**
def load_disease_seeds(disease_dirs):
    """
    Reads `seeds.txt` from all disease directories **including Astrocytes**.
    Stores **Gene Symbol** for each disease.
    Returns a dictionary mapping each seed gene symbol to the disease(s) where it appears.
    """
    seed_genes = {}

    for disease_name, disease_path in disease_dirs.items():
        seed_file = os.path.join(disease_path, "seeds.txt")

        # **Skip if no `seeds.txt` found**
        if not os.path.exists(seed_file):
            print("⚠️ WARNING: `seeds.txt` not found in {}".format(disease_path))
            continue

        with open(seed_file, 'r') as f:
            for line in f:
                cols = line.strip().split('\t')

                # **Ensure valid `seeds.txt` format**
                if len(cols) < 3:
                    continue  # Skip bad lines
                
                gene_symbol = cols[2].strip()  # **Column 3: Gene Symbol**

                # **Store gene symbol as a seed**
                if gene_symbol:
                    if gene_symbol not in seed_genes:
                        seed_genes[gene_symbol] = set()
                    seed_genes[gene_symbol].add(disease_name)  # Store disease name

    return seed_genes

### **2️⃣ Load GUILD Scores**
def load_scores(file_path):
    """
    Reads `guild_scores.txt` and returns a list of gene dictionaries.
    """
    scores = []
    with open(file_path, 'r') as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            row = {k.strip(): v.strip() for k, v in row.iteritems()}  # Preserve original case
            scores.append(row)
    return scores

### **3️⃣ Get Top N Genes by Score**
def get_top_proteins(scores, p_top):
    """
    Extracts the top p_top% genes based on 'GUILD Score'.
    """
    scores.sort(key=lambda x: float(x["GUILD Score"]), reverse=True)  # Sort descending by GUILD Score
    n = int(round((p_top / 100.0) * len(scores)))
    return scores[:n]

### **4️⃣ Main Execution**
def main():
    if len(sys.argv) < 2:
        print("Usage: python overlap.py <top_percentage>")
        sys.exit(1)

    top_percentage = float(sys.argv[1])

    # **Automatically locate protein files from DISEASE_DIRS**
    protein_files = []
    for disease_name, disease_path in DISEASE_DIRS.items():
        protein_file = os.path.join(disease_path, "guild_scores.txt")
        if not os.path.exists(protein_file):
            print("⚠️ WARNING: 'guild_scores.txt' not found in {} for {}".format(disease_path, disease_name))
            continue
        protein_files.append(protein_file)

    if len(protein_files) < 2:
        print("Error: At least two protein files are required, but only found {}.".format(len(protein_files)))
        sys.exit(1)

    # **Load all disease seed genes (INCLUDING ASTROCYTES)**
    seed_genes = load_disease_seeds(DISEASE_DIRS)

    # **Load all protein scores**
    all_scores = [load_scores(f) for f in protein_files]

    # **Get top proteins for each phenotype**
    top_proteins_list = [get_top_proteins(scores, top_percentage) for scores in all_scores]

    # **Find intersection of top genes across ALL phenotypes (INCLUDING ASTROCYTES)**
    top_sets = [set([gene["Gene Symbol"] for gene in df]) for df in top_proteins_list]
    common_proteins = set.intersection(*top_sets)

    # **Map Gene Symbol to UniProt ID**
    gene_info = {}
    for df in top_proteins_list:
        for gene in df:
            gene_info[gene["Gene Symbol"]] = {
                "Gene ID": gene["Gene ID"],
                "UniProt ID": gene["UniProt ID"]
            }

    # **Save corrected overlap results**
    output_file = "top5_overlap_results.txt"
    with open(output_file, "w") as f:
        f.write("Gene ID\tUniProt ID\tGene Symbol\tSeed Status\tSeed Source\n")
        for gene_symbol in common_proteins:
            gene_symbol_original = gene_symbol.strip()  # Keep original case
            gene_id = gene_info[gene_symbol_original]["Gene ID"] if gene_symbol_original in gene_info else "NA"
            uniprot_id = gene_info[gene_symbol_original]["UniProt ID"] if gene_symbol_original in gene_info else "NA"

            # **Check if Gene Symbol is in Seeds (INCLUDING ASTROCYTES)**
            matching_diseases = seed_genes.get(gene_symbol_original, set())

            if matching_diseases:
                seed_status = "seed"
                seed_source = ";".join(sorted(matching_diseases))  # List of diseases where it's a seed
            else:
                seed_status = "non-seed"
                seed_source = "None"

            f.write("{}\t{}\t{}\t{}\t{}\n".format(gene_id, uniprot_id, gene_symbol_original, seed_status, seed_source))

    print("✅ Top {}% overlap results saved to `{}`".format(top_percentage, output_file))

if __name__ == "__main__":
    main()
