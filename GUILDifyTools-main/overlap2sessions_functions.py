#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

import pandas as pd
from scipy.stats import fisher_exact, combine_pvalues
import sys
import os

##############################
# Define disease directories #
##############################
DISEASE_DIRS = {
    "Alzheimer": "astrocytes/alzheimer/",
    "Cardiovascular": "astrocytes/cardiovascular/",
    "Diabetes": "astrocytes/diabetes/",
    "Inflammation": "astrocytes/inflammation/",
    "Oxidative_Stress": "astrocytes/oxidative_stress/",
    "Cognitive_Impairment": "astrocytes/cognitive_impairment/",
    "Astrocytes": "astrocytes/astrocytes_overexpressed"  # Astrocytes now included
}

##############################
# File names expected in each directory
##############################
PROTEIN_FILENAME = "guild_scores.txt"
ENRICHMENT_FILENAME = "enrichment.GOmf.fdr_bh.5.txt"
SEEDS_FILENAME = "seeds.txt"  # used for seed annotation

#########################################
# Load disease seeds from seeds.txt files
#########################################
def load_disease_seeds(disease_dirs):
    """
    Reads `seeds.txt` from all disease directories (including Astrocytes).
    Returns a dictionary mapping each seed gene symbol to a set of disease names.
    """
    seed_genes = {}
    for disease_name, disease_path in disease_dirs.items():
        seed_file = os.path.join(disease_path, SEEDS_FILENAME)
        if not os.path.exists(seed_file):
            print("⚠️ WARNING: `seeds.txt` not found in {}".format(disease_path))
            continue
        with open(seed_file, 'r') as f:
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) < 3:
                    continue  # Skip improperly formatted lines
                gene_symbol = cols[2].strip()  # Column 3: Gene Symbol
                if gene_symbol:
                    if gene_symbol not in seed_genes:
                        seed_genes[gene_symbol] = set()
                    seed_genes[gene_symbol].add(disease_name)
    return seed_genes

#########################################
# Load protein scores
#########################################
def load_scores(file_path):
    scores = pd.read_csv(file_path, sep="\t")
    scores.columns = scores.columns.str.strip()  # Clean column names
    return scores

#########################################
# Load GO enrichment data
#########################################
def load_enrichment(file_path):
    enrichment = pd.read_csv(file_path, sep="\t")
    enrichment.columns = enrichment.columns.str.strip()
    enrichment.columns = [col.lstrip('#').strip() for col in enrichment.columns]  # Clean column names
    return enrichment

#########################################
# Get top N proteins based on GUILD Score
#########################################
def get_top_proteins(scores, p_top):
    scores = scores.sort_values(by="GUILD Score", ascending=False)
    n = int(round((p_top / 100.0) * len(scores)))
    return scores.iloc[:n]

#########################################
# Select significant terms from enrichment data
#########################################
def get_top_enrichment(enrichment, pvalue_cutoff=0.05):
    enrichment = enrichment[enrichment['P-value corrected'] <= pvalue_cutoff]
    return enrichment

#########################################
# Fisher's exact test for overlap significance
#########################################
def calculate_fisher(n_common, n1, n2, n_total):
    contingency = [[n_common, n1 - n_common], [n2 - n_common, n_total - n1 - n2 + n_common]]
    oddsratio, pvalue = fisher_exact(contingency, alternative="greater")
    return oddsratio, pvalue

#########################################
# Annotate common proteins using gene symbols and seed data
#########################################
def annotate_common_proteins(session, common_gene_symbols, seed_genes):
    """
    Iterates over the given session (DataFrame) and, for each protein whose 'Gene Symbol'
    is in the common_gene_symbols set, annotates it with seed status and source.
    Returns a list of rows: [Gene ID, UniProt ID, Gene Symbol, Seed Status, Seed Source].
    """
    # Build a mapping from Gene Symbol to (Gene ID, UniProt ID) from the session
    gene_info = {}
    for _, row in session.iterrows():
        gs = row.get("Gene Symbol", "").strip()
        if gs and gs in common_gene_symbols:
            gene_info[gs] = {
                "Gene ID": row.get("Gene ID", "NA"),
                "UniProt ID": row.get("UniProt ID", "NA")
            }
    annotated = []
    for gs in sorted(common_gene_symbols):
        info = gene_info.get(gs, {"Gene ID": "NA", "UniProt ID": "NA"})
        if gs in seed_genes:
            seed_status = "seed"
            seed_source = ";".join(sorted(seed_genes[gs]))
        else:
            seed_status = "non-seed"
            seed_source = "None"
        annotated.append([info["Gene ID"], info["UniProt ID"], gs, seed_status, seed_source])
    return annotated

#########################################
# Main execution function
#########################################
def main():
    if len(sys.argv) < 2:
        print("Usage: python {} <top_percentage>".format(sys.argv[0]))
        sys.exit(1)

    top_percentage = float(sys.argv[1])
    
    # Automatically build file lists from DISEASE_DIRS
    protein_files = []
    enrichment_files = []
    for disease_name, disease_path in DISEASE_DIRS.items():
        protein_file = os.path.join(disease_path, PROTEIN_FILENAME)
        enrichment_file = os.path.join(disease_path, ENRICHMENT_FILENAME)
        if not os.path.exists(protein_file):
            print("⚠️ WARNING: Protein file not found for {}: {}".format(disease_name, protein_file))
        else:
            protein_files.append(protein_file)
        if not os.path.exists(enrichment_file):
            print("⚠️ WARNING: Enrichment file not found for {}: {}".format(disease_name, enrichment_file))
        else:
            enrichment_files.append(enrichment_file)

    if len(protein_files) < 1:
        print("At least one protein file is required.")
        sys.exit(1)
    if len(enrichment_files) < 1:
        print("At least one enrichment file is required.")
        sys.exit(1)

    # Load all protein scores
    all_scores = [load_scores(f) for f in protein_files]

    # Get top proteins for each dataset
    top_proteins_list = [get_top_proteins(scores, top_percentage) for scores in all_scores]

    # Compute overlap based on UniProt IDs (for statistics)
    top_sets = [set(df["UniProt ID"]) for df in top_proteins_list]
    common_proteins = set.intersection(*top_sets)
    n_common_proteins = len(common_proteins)

    # Compute common gene symbols (used for seed annotation) if available
    gene_sets = [set(df["Gene Symbol"].dropna().str.strip()) for df in top_proteins_list if "Gene Symbol" in df.columns]
    if gene_sets:
        common_gene_symbols = set.intersection(*gene_sets)
    else:
        common_gene_symbols = set()

    # Load seeds from disease directories
    seed_genes = load_disease_seeds(DISEASE_DIRS)

    # Protein overlap analysis
    if len(protein_files) == 2:
        # Fisher test for proteins (only if exactly 2 sessions)
        scores1, scores2 = all_scores
        top_proteins1, top_proteins2 = top_proteins_list
        n1 = len(top_proteins1)
        n2 = len(top_proteins2)
        n_total_proteins = len(set(scores1["UniProt ID"]) | set(scores2["UniProt ID"]))
        oddsratio_proteins, pvalue_proteins = calculate_fisher(n_common_proteins, n1, n2, n_total_proteins)

        overlap_file = "protein_overlap_results.txt"
        stats_file = "protein_overlap_stats.txt"
        with open(overlap_file, "w") as f:
            f.write("Gene ID\tUniProt ID\tGene Symbol\tSeed Status\tSeed Source\n")
            annotated_common = annotate_common_proteins(top_proteins_list[0], common_gene_symbols, seed_genes)
            for row in annotated_common:
                f.write("\t".join(map(str, row)) + "\n")

        with open(stats_file, "w") as f:
            f.write("N common (UniProt IDs)\tN total\tN1\tN2\tOdds Ratio\tP value\n")
            f.write("{}\t{}\t{}\t{}\t{:.3f}\t{:.1E}\n".format(
                n_common_proteins, n_total_proteins, n1, n2, oddsratio_proteins, pvalue_proteins))
        print("Protein overlap results written to", overlap_file)
        print("Protein statistical summary written to", stats_file)
    else:
        # For more than two sessions, output common proteins without Fisher's test
        overlap_file = "protein_overlap_results.txt"
        with open(overlap_file, "w") as f:
            f.write("Gene ID\tUniProt ID\tGene Symbol\tSeed Status\tSeed Source\n")
            annotated_common = annotate_common_proteins(top_proteins_list[0], common_gene_symbols, seed_genes)
            for row in annotated_common:
                f.write("\t".join(map(str, row)) + "\n")
        print("Protein overlap results (no Fisher's test for >2 sessions) written to", overlap_file)

    ##############################
    # Enrichment analysis section
    ##############################
    # Fixed number of GO terms (assuming GObp)
    num_GOs = {
        "GObp": 9390,
        "Reactome": 2150,
        "GOmf": 3550
    }
    n_enrich = len(enrichment_files)
    all_enrichments = [load_enrichment(f) for f in enrichment_files]

    # Optionally filter enrichment files by top proteins if possible
    filter_by_genes = ("Gene ID" in all_enrichments[0].columns) and len(protein_files) == 2
    if filter_by_genes:
        top_genes1 = set(top_proteins_list[0].get("Gene ID", []))
        top_genes2 = set(top_proteins_list[1].get("Gene ID", []))
        if len(protein_files) > 2:
            gene_sets = []
            for tp in top_proteins_list:
                if "Gene ID" in tp.columns:
                    gene_sets.append(set(tp["Gene ID"]))
                else:
                    gene_sets = []
                    break
            if gene_sets:
                top_genes_intersect = set.intersection(*gene_sets)
            else:
                top_genes_intersect = top_genes1 & top_genes2
        else:
            top_genes_intersect = top_genes1 & top_genes2

        filtered_enrichments = []
        for ef in all_enrichments:
            filtered_enrichments.append(ef[ef["Gene ID"].isin(top_genes_intersect)])
        all_enrichments = filtered_enrichments
    else:
        if any("Gene ID" not in ef.columns for ef in all_enrichments):
            print("Warning: 'Gene ID' column not found in some enrichment files. Skipping filtering by top proteins.")

    # Get top significant enrichment terms
    top_enrichments = [get_top_enrichment(ef) for ef in all_enrichments]

    # Intersect terms across all enrichment sets
    term_sets = [set(e["Term ID"]) for e in top_enrichments]
    common_functions = set.intersection(*term_sets)
    n_common_functions = len(common_functions)

    # Count functions in first set for category (assuming "GObp")
    n_total_functions = num_GOs["GObp"]
    if n_enrich == 2 and len(common_functions) > 0:
        e1, e2 = top_enrichments
        n1_functions = len(e1)
        n2_functions = len(e2)
        oddsratio_functions, pvalue_functions = calculate_fisher(n_common_functions, n1_functions, n2_functions, n_total_functions)

        function_overlap_file = "function_overlap_results.txt"
        function_stats_file = "function_overlap_stats.txt"

        with open(function_overlap_file, "w") as f:
            f.write("Term ID\tTerm name\t" + "\t".join("P-value {}".format(i+1) for i in range(n_enrich)) + "\tCombined P-value\n")
            for term_id in common_functions:
                pvals = []
                term_name = None
                for ef in top_enrichments:
                    row = ef[ef["Term ID"] == term_id]
                    pval = row["P-value corrected"].values[0]
                    pvals.append(pval)
                    if term_name is None:
                        term_name = row["Term name"].values[0]
                _, combined_pval = combine_pvalues(pvals, method='fisher')
                f.write("{}\t{}\t{}\t{:.2E}\n".format(
                    term_id,
                    term_name,
                    "\t".join("{:.2E}".format(p) for p in pvals),
                    combined_pval
                ))

        with open(function_stats_file, "w") as f:
            f.write("N common\tN total\tN1\tN2\tOdds Ratio\tP value\n")
            f.write("{}\t{}\t{}\t{}\t{:.3f}\t{:.1E}\n".format(
                n_common_functions, n_total_functions, n1_functions, n2_functions, oddsratio_functions, pvalue_functions))
        print("Function overlap results written to", function_overlap_file)
        print("Function statistical summary written to", function_stats_file)
    elif n_enrich == 2:
        function_overlap_file = "function_overlap_results.txt"
        with open(function_overlap_file, "w") as f:
            f.write("No common functions found.\n")
        print("No common functions, no Fisher test performed.")
    else:
        function_overlap_file = "function_overlap_results.txt"
        with open(function_overlap_file, "w") as f:
            f.write("Term ID\tTerm name\n")
            first_e = top_enrichments[0]
            for term_id in common_functions:
                row = first_e[first_e["Term ID"] == term_id]
                term_name = row["Term name"].values[0]
                f.write("{}\t{}\n".format(term_id, term_name))
        print("Function overlap results (no Fisher's test for >2 enrichment sets) written to", function_overlap_file)

if __name__ == "__main__":
    main()
