# Astrocytes Project Data

This repository contains data and scripts used for the astrocytes-related analyses supporting an upcoming article.

## Contents

- `BIANA_phy/`: Network and annotation data used for analyses.
  - Note: A few extremely large files are intentionally excluded due to GitHub file-size limits:
    - `BIANA_phy/all_shortest_paths.txt` (~6 GB)
    - `BIANA_phy/multi-fields_file.txt` (~154 MB)
    - `BIANA_phy/network.method.txt` (~160 MB)
    - `BIANA_phy/drugbank/` (contains files up to ~1.5 GB)
    - `BIANA_phy/go/` (contains files up to ~8.4 GB)
    - `BIANA_phy/sampled_graphs/` (~1.1 GB total)
  - If you need these files, please see the section below on Large Files.

- `GUILDifyTools-main/`: Selected tools and datasets to reproduce overlap and scoring steps:
  - `astrocytes/` (all contents)
  - `clinvar_top5.txt`
  - `overlap2sessions_functions.py`
  - `README.md`
  - `Datos_Genomica_Astrocitos_Baldo.xlsx`
  - `overlap.py`
  - `run_guildify_local.py`

## Large Files

GitHub limits individual files to 100 MB and strongly discourages multi‑GB assets in normal Git history. The excluded files listed above exceed those limits. They should be hosted externally (e.g., Zenodo/OSF/Figshare) and linked here. Once a DOI/URL is available, add it below:

- BIANA large files download: <ADD-LINK-HERE>

## Getting Started

1. Clone the repository.
2. Review `GUILDifyTools-main/README.md` for tool-specific notes.
3. Example entry points:
   - `GUILDifyTools-main/run_guildify_local.py`
   - `GUILDifyTools-main/overlap.py`

## Notes

- Only the subset of `GUILDifyTools-main` required for reproduction is included. Other folders (e.g., `old_astrocytes_wrong/`, `quim_astro/`, `scripts/`) are intentionally excluded.
- No license has been specified yet. Please add one if redistribution is intended.
