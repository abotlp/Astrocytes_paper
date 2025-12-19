#!/usr/bin/env python3
"""
Compare our astrocyte H2O2 DEGs (from provided Excel tables) vs Limbad et al. 2020
human astrocyte senescence signature (PLOS ONE supplemental XLSX: journal.pone.0227887.s008).

Outputs:
  - Limbad2020_Senescence_Up_SUR.csv
  - Limbad2020_Senescence_Down_SDR.csv
  - overlap_our_DEGs_vs_Limbad2020_senescence.csv
"""

from __future__ import annotations
import argparse
import io
import sys
import textwrap
from pathlib import Path
from typing import Iterable, Tuple, Set, Dict

import pandas as pd
import requests
from scipy.stats import fisher_exact


DEFAULT_LIMBAD_SUPP_DOI = "https://doi.org/10.1371/journal.pone.0227887.s008"


def _norm_gene(x) -> str | None:
    if x is None:
        return None
    s = str(x).strip()
    if not s or s.lower() in {"nan", "none"}:
        return None
    # normalize to HGNC-like uppercase symbols (your files already are)
    return s.upper()


def read_gene_set_from_xlsx(
    xlsx_path: Path,
    sheet: str | int = 0,
    gene_col: str = "Gene Name",
    padj_col: str = "Adj. P.Val",
    padj_threshold: float | None = 0.1,
) -> Set[str]:
    df = pd.read_excel(xlsx_path, sheet_name=sheet)
    if gene_col not in df.columns:
        raise ValueError(f"{xlsx_path}: gene_col='{gene_col}' not found. Columns: {list(df.columns)}")
    genes = df[gene_col].map(_norm_gene).dropna()

    if padj_threshold is not None:
        if padj_col not in df.columns:
            raise ValueError(f"{xlsx_path}: padj_col='{padj_col}' not found. Columns: {list(df.columns)}")
        padj = pd.to_numeric(df[padj_col], errors="coerce")
        keep = padj < padj_threshold
        genes = df.loc[keep, gene_col].map(_norm_gene).dropna()

    return set(genes.tolist())


def download_limbad_supp_xlsx(url_or_doi: str, out_path: Path) -> Path:
    """
    Downloads the PLOS supplemental XLSX. If given a DOI URL, requests will follow redirects.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    resp = requests.get(url_or_doi, allow_redirects=True, timeout=60)
    resp.raise_for_status()
    out_path.write_bytes(resp.content)
    return out_path


def load_limbad_signature(xlsx_path: Path) -> Tuple[Set[str], Set[str], Set[str]]:
    """
    Returns: (SUR_up, SDR_down, universe_genes)
    - SUR and SDR are sheet names in the Limbad supplemental file.
    - 'diff' sheet is used for a background/universe if present.
    """
    xl = pd.ExcelFile(xlsx_path)

    # helper: first column that looks like a gene symbol list
    def read_sheet_genes(sheet_name: str) -> Set[str]:
        df = pd.read_excel(xlsx_path, sheet_name=sheet_name)
        # pick a likely gene column
        candidates = [c for c in df.columns if str(c).lower() in {"gene", "genes", "symbol", "hgnc_symbol", "gene name"}]
        gene_col = candidates[0] if candidates else df.columns[0]
        return set(df[gene_col].map(_norm_gene).dropna().tolist())

    if "SUR" not in xl.sheet_names or "SDR" not in xl.sheet_names:
        raise ValueError(f"Limbad XLSX missing SUR/SDR sheets. Found: {xl.sheet_names}")

    sur = read_sheet_genes("SUR")
    sdr = read_sheet_genes("SDR")

    # background/universe: prefer 'diff' sheet, else fall back to union
    if "diff" in xl.sheet_names:
        diff_df = pd.read_excel(xlsx_path, sheet_name="diff")
        candidates = [c for c in diff_df.columns if str(c).lower() in {"gene", "genes", "symbol", "hgnc_symbol", "gene name"}]
        gene_col = candidates[0] if candidates else diff_df.columns[0]
        universe = set(diff_df[gene_col].map(_norm_gene).dropna().tolist())
    else:
        universe = set(sur) | set(sdr)

    return sur, sdr, universe


def fisher_enrichment(universe: Set[str], our: Set[str], sig: Set[str]) -> Dict[str, float | int | str]:
    """
    One-sided Fisher exact test for enrichment of overlap between 'our' and 'sig'.
    Contingency:
        in_sig   not_in_sig
    in_our     x         K-x
    not_in_our M-x       N-K-M+x
    """
    N = len(universe)
    our_u = our & universe
    sig_u = sig & universe

    K = len(our_u)
    M = len(sig_u)
    x = len(our_u & sig_u)

    # protect against negative cell
    a = x
    b = K - x
    c = M - x
    d = N - (a + b + c)
    if d < 0:
        raise ValueError(f"Invalid contingency table: N={N} K={K} M={M} x={x} gives d={d}")

    odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    return {
        "N_universe": N,
        "K_our": K,
        "M_signature": M,
        "x_overlap": x,
        "odds_ratio": float(odds) if odds == odds else float("nan"),
        "p_value_greater": float(p),
    }


def main():
    ap = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Compare our DEGs (Excel) vs Limbad 2020 astrocyte senescence (SUR/SDR).",
        epilog=textwrap.dedent("""\
        Example:
          python3 compare_DEGs_vs_Limbad2020_read_our_xlsx.py \\
            --our-up Table_1_Overespressed_Genes.xlsx \\
            --our-down Table_2_Repressed_Genes.xlsx
        """),
    )
    ap.add_argument("--our-up", default="Table_1_Overespressed_Genes.xlsx", help="Our UP-regulated Excel file")
    ap.add_argument("--our-down", default="Table_2_Repressed_Genes.xlsx", help="Our DOWN-regulated Excel file")
    ap.add_argument("--our-sheet", default="Sheet1", help="Sheet name (or 0) in our Excel files")
    ap.add_argument("--gene-col", default="Gene Name", help="Column containing gene symbols in our Excel files")
    ap.add_argument("--padj-col", default="Adj. P.Val", help="Column containing adjusted p-values in our Excel files")
    ap.add_argument("--padj", type=float, default=0.1, help="Adjusted p-value cutoff (set to -1 to disable filtering)")
    ap.add_argument("--limbad-doi", default=DEFAULT_LIMBAD_SUPP_DOI, help="Limbad supplemental XLSX DOI/URL")
    ap.add_argument("--limbad-cache", default="Limbad2020_supplement.xlsx", help="Local filename for downloaded Limbad XLSX")
    args = ap.parse_args()

    our_up_path = Path(args.our_up)
    our_down_path = Path(args.our_down)
    our_sheet = args.our_sheet
    padj_thr = None if args.padj < 0 else args.padj

    OUR_UP = read_gene_set_from_xlsx(
        our_up_path, sheet=our_sheet, gene_col=args.gene_col, padj_col=args.padj_col, padj_threshold=padj_thr
    )
    OUR_DOWN = read_gene_set_from_xlsx(
        our_down_path, sheet=our_sheet, gene_col=args.gene_col, padj_col=args.padj_col, padj_threshold=padj_thr
    )
    OUR_ALL = OUR_UP | OUR_DOWN

    # Download Limbad XLSX (cache on disk so reruns are fast/offline)
    limbad_xlsx = Path(args.limbad_cache)
    if not limbad_xlsx.exists():
        download_limbad_supp_xlsx(args.limbad_doi, limbad_xlsx)

    SUR, SDR, UNIV = load_limbad_signature(limbad_xlsx)
    LIMBAD_ALL = (SUR | SDR) & UNIV

    # Save SUR/SDR for transparency
    pd.DataFrame({"gene": sorted(SUR)}).to_csv("Limbad2020_Senescence_Up_SUR.csv", index=False)
    pd.DataFrame({"gene": sorted(SDR)}).to_csv("Limbad2020_Senescence_Down_SDR.csv", index=False)

    # Compute enrichment stats
    any_stats = fisher_enrichment(UNIV, OUR_ALL, LIMBAD_ALL)
    up_stats = fisher_enrichment(UNIV, OUR_UP, SUR)
    down_stats = fisher_enrichment(UNIV, OUR_DOWN, SDR)

    def fmt_stats(title: str, stats: Dict[str, float | int | str], genes: Iterable[str]):
        genes = sorted(set(genes))
        print(f"\n{title}")
        print(f"  Universe N={stats['N_universe']}")
        print(f"  our_set K={stats['K_our']}, signature M={stats['M_signature']}, overlap x={stats['x_overlap']}")
        print(f"  Fisher exact (greater): odds={stats['odds_ratio']:.3g}, p={stats['p_value_greater']:.3g}")
        if genes:
            print("  Overlap genes: " + ", ".join(genes))
        else:
            print("  Overlap genes: (none)")

    # Overlap gene lists
    any_overlap = (OUR_ALL & LIMBAD_ALL & UNIV)
    up_overlap = (OUR_UP & SUR & UNIV)
    down_overlap = (OUR_DOWN & SDR & UNIV)

    fmt_stats("Any-direction overlap (OUR_DEGs ∩ Limbad2020_DEGs)", any_stats, any_overlap)
    fmt_stats("Direction-matched UP (OUR_UP ∩ Limbad2020_UP)", up_stats, up_overlap)
    fmt_stats("Direction-matched DOWN (OUR_DOWN ∩ Limbad2020_DOWN)", down_stats, down_overlap)

    # Write one combined CSV for manuscript/rebuttal use
    out = pd.DataFrame(
        {
            "comparison": ["any_direction", "up_vs_SUR", "down_vs_SDR"],
            "N_universe": [any_stats["N_universe"], up_stats["N_universe"], down_stats["N_universe"]],
            "K_our": [any_stats["K_our"], up_stats["K_our"], down_stats["K_our"]],
            "M_signature": [any_stats["M_signature"], up_stats["M_signature"], down_stats["M_signature"]],
            "x_overlap": [any_stats["x_overlap"], up_stats["x_overlap"], down_stats["x_overlap"]],
            "odds_ratio": [any_stats["odds_ratio"], up_stats["odds_ratio"], down_stats["odds_ratio"]],
            "p_value_greater": [any_stats["p_value_greater"], up_stats["p_value_greater"], down_stats["p_value_greater"]],
            "overlap_genes": [
                ";".join(sorted(any_overlap)),
                ";".join(sorted(up_overlap)),
                ";".join(sorted(down_overlap)),
            ],
        }
    )
    out.to_csv("overlap_our_DEGs_vs_Limbad2020_senescence.csv", index=False)

    print("\nWrote CSV outputs in current directory.")


if __name__ == "__main__":
    main()
