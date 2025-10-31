# Astrocytes Paper — Data, Code and Reproducibility

Data and code are publicly available at https://github.com/abotlp/Astrocytes_paper. All analysis outputs are organized under `Astrocytes_project/GUILDifyTools-main/astrocytes` and are built on the BIANA‑derived human interactome in `Astrocytes_project/BIANA_phy`.

## Run GUILDify locally (seeds → ranked network + enrichment)

Run from the repository root with repo‑relative paths. Provide HGNC gene symbols separated by semicolons.

Example (generalized):

```
python GUILDifyTools-main/scripts/run_guildify_local.py \
  -i "SEED1; SEED2; SEED3" \
  -n BIANA_phy \
  -o GUILDifyTools-main/local_guildify/<run_name> \
  -s netcombo -re 3 -it 2 \
  -t 9606 -ns BIANA_phy -ti all \
  -dr BIANA_phy/drug_info.txt
```

Notes
- `-i`: semicolon‑separated HGNC symbols (e.g., for astrocytes overexpressed).
- `-n`: network directory (here, `BIANA_phy`).
- `-o`: output folder created for this run.
- `-s`: scoring method (`netcombo`, `netscore`, `netzcore`, `netshort`, or `diamond`).
- `-re`, `-it`: repetitions/iterations for `netscore`/`netcombo`/`netzcore`.
- `-t`: taxonomy ID (9606 for human).
- `-ns`: network source label (`BIANA_phy`).
- `-ti`: tissue filter (`all`).
- `-dr`: drug–target info file path (repo‑relative).

Dependencies
- Python 3, `numpy`, `networkx`, `scipy`, `statsmodels`.

## Where to find key outputs

Per‑run outputs (phenotypes and astrocyte seeds)
- Folder: `GUILDifyTools-main/astrocytes/<run>` where `<run>` is one of:
  - `alzheimer`, `cardiovascular`, `cognitive_impairment`, `diabetes`, `inflammation`, `oxidative_stress`
  - `astrocytes_overexpressed`, `astrocytes_underexpressed`
- Key files:
  - `seeds.txt` — seed genes actually recognized in the network
  - `guild_scores.txt` — consensus ranking across methods
  - Enrichment: `enrichment.GObp.*`, `enrichment.GOmf.*`, `enrichment.Reactome.*` (Bonferroni and FDR)
  - Top‑k subnetworks: `subnetwork.sif.{1,2,5}` and corresponding `network_plot.png{1,2,5}` (when present)

Cross‑network overlaps (shared cores)
- Folders: `GUILDifyTools-main/astrocytes/network_overlap_over5`, `.../network_overlap_under5`, `.../network_overlap_over2`, `.../network_overlap_under2`
- Kept files (cleaned for publication):
  - `protein_overlap_results.txt` — shared gene list with seed status/source
  - `function_overlap_results.txt` — Reactome pathway overlap for the shared set
  - `subnetwork_genes.sif`, `subnetwork_genes.json` — core overlap network

Mutation analyses (FoldX ΔΔG + ClinVar)
- Folder: `ddg_analysis/`
  - `ddg_docking.csv` — docking‑derived complex ΔΔG summary
  - `ddg_homology.csv` — homology/template‑derived complex ΔΔG summary
  - `clinvar.csv` — ClinVar variant table used for mapping

BIANA network assets
- Folder: `BIANA_phy/` — node/edge tables, annotations, and resources used by GUILDify.
- Some very large files are not tracked in Git; host externally (e.g., Zenodo/OSF) and add the link here when available.

## Citation and contact
Please cite the associated manuscript and this repository if you use the data or code. For questions or issues, open a GitHub issue in this repo.
