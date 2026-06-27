"""
rbcL Sequence Analysis Pipeline - Advanced
==========================================
Fitur:
  1. Membaca file FASTA multi-sekuens
  2. Menghitung frekuensi nukleotida & GC Content
  3. Mengurutkan sekuens berdasarkan GC Content (top 3)
  4. Pairwise Similarity Matrix (% identity)
  5. Motif / Restriction Site Searching (IUPAC-aware, both strands)
  6. Visualisasi: bar chart GC, heatmap similarity, motif hit chart
  7. Ekspor hasil ke CSV

Kebutuhan:
  pip install biopython pandas matplotlib seaborn

Penggunaan:
  python pipeline_advanced.py
  python pipeline_advanced.py --fasta myfile.fasta --motifs GAATTC AAGCTT
"""

import argparse
import re
import sys
from pathlib import Path
from itertools import combinations

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from Bio import SeqIO

# ─────────────────────────────────────────────
# IUPAC ambiguity code → regex
# ─────────────────────────────────────────────
IUPAC_MAP = {
    'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'U': 'T',
    'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
}

COMPLEMENT = str.maketrans('ATGCRYSWKMBDHVN', 'TACGYRSWMKVHDBN')

PRESET_MOTIFS = {
    'EcoRI'   : 'GAATTC',
    'HindIII' : 'AAGCTT',
    'BamHI'   : 'GGATCC',
    'SmaI'    : 'CCCGGG',
    'XhoI'    : 'CTCGAG',
    'SphI'    : 'GCATGC',
    'TATA_box': 'TATAAA',
}


def iupac_to_regex(motif: str) -> str:
    """Konversi motif IUPAC ke string regex."""
    return ''.join(IUPAC_MAP.get(c, c) for c in motif.upper())


def reverse_complement(seq: str) -> str:
    return seq.upper().translate(COMPLEMENT)[::-1]


# ─────────────────────────────────────────────
# 1. BACA FASTA
# ─────────────────────────────────────────────
def load_fasta(filepath: str) -> list:
    records = list(SeqIO.parse(filepath, 'fasta'))
    if not records:
        sys.exit(f'[ERROR] Tidak ada sekuens dalam file: {filepath}')
    print(f'[INFO] Jumlah sekuens: {len(records)}')
    return records


# ─────────────────────────────────────────────
# 2. ANALISIS NUKLEOTIDA & GC CONTENT
# ─────────────────────────────────────────────
def analyze_gc(records: list) -> pd.DataFrame:
    hasil = []
    for record in records:
        seq = str(record.seq).upper()
        freq = {b: seq.count(b) for b in 'ATGC'}
        length = len(seq)
        gc = (freq['G'] + freq['C']) / length * 100 if length else 0
        hasil.append({
            'ID'          : record.id,
            'Panjang (bp)': length,
            'A'           : freq['A'],
            'T'           : freq['T'],
            'G'           : freq['G'],
            'C'           : freq['C'],
            'GC_Content'  : round(gc, 2),
        })
    df = pd.DataFrame(hasil).sort_values('GC_Content', ascending=False).reset_index(drop=True)
    return df


# ─────────────────────────────────────────────
# 3. PAIRWISE SIMILARITY MATRIX
# ─────────────────────────────────────────────
def pairwise_identity(seqA: str, seqB: str) -> float:
    """% identity = matching positions / min(len) * 100 (simple positional)."""
    length = min(len(seqA), len(seqB))
    if length == 0:
        return 0.0
    matches = sum(a == b for a, b in zip(seqA[:length], seqB[:length]))
    return round(matches / length * 100, 2)


def build_similarity_matrix(records: list) -> pd.DataFrame:
    ids  = [r.id for r in records]
    seqs = [str(r.seq).upper() for r in records]
    n    = len(records)

    # Inisialisasi matrix
    matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        matrix[i][i] = 100.0
    for i, j in combinations(range(n), 2):
        val = pairwise_identity(seqs[i], seqs[j])
        matrix[i][j] = val
        matrix[j][i] = val

    df = pd.DataFrame(matrix, index=ids, columns=ids)
    return df


# ─────────────────────────────────────────────
# 4. MOTIF / RESTRICTION SITE SEARCH
# ─────────────────────────────────────────────
def search_motif(records: list, motif: str, both_strands: bool = True) -> pd.DataFrame:
    motif = motif.upper()
    pattern_fwd = re.compile(iupac_to_regex(motif))
    rc_motif    = reverse_complement(motif)
    pattern_rev = re.compile(iupac_to_regex(rc_motif)) if both_strands else None

    rows = []
    for record in records:
        seq = str(record.seq).upper()

        # Forward strand — non-overlapping search
        fwd_hits = [m.start() + 1 for m in re.finditer(f'(?={iupac_to_regex(motif)})', seq)]

        rev_hits = []
        if both_strands and rc_motif != motif:
            rev_hits = [m.start() + 1 for m in re.finditer(f'(?={iupac_to_regex(rc_motif)})', seq)]

        total = len(fwd_hits) + len(rev_hits)
        rows.append({
            'ID'               : record.id,
            'Motif'            : motif,
            'Total_Hits'       : total,
            'Forward_Hits'     : len(fwd_hits),
            'Reverse_Hits'     : len(rev_hits),
            'Forward_Positions': ', '.join(str(p) for p in fwd_hits[:10]) + ('…' if len(fwd_hits) > 10 else ''),
            'Reverse_Positions': ', '.join(str(p) for p in rev_hits[:10]) + ('…' if len(rev_hits) > 10 else ''),
        })

    df = pd.DataFrame(rows).sort_values('Total_Hits', ascending=False).reset_index(drop=True)
    return df


# ─────────────────────────────────────────────
# 5. VISUALISASI
# ─────────────────────────────────────────────
def plot_gc_bar(df: pd.DataFrame, out: str = 'gc_content_rbcl.png'):
    fig, ax = plt.subplots(figsize=(max(10, len(df) * 0.5), 5))

    colors = ['#39d353' if gc >= 50 else '#58a6ff' if gc >= 40 else '#d29922'
              for gc in df['GC_Content']]
    ax.bar(df['ID'], df['GC_Content'], color=colors, edgecolor='none', width=0.7)

    ax.set_ylabel('GC Content (%)', fontsize=11)
    ax.set_xlabel('Sequence ID', fontsize=11)
    ax.set_title('GC Content per Sekuens (diurutkan descending)', fontsize=13, fontweight='bold')
    ax.set_ylim(0, 100)
    ax.tick_params(axis='x', rotation=90, labelsize=7)
    ax.grid(axis='y', linestyle='--', alpha=0.4)

    legend_patches = [
        mpatches.Patch(color='#39d353', label='GC ≥ 50%'),
        mpatches.Patch(color='#58a6ff', label='40% ≤ GC < 50%'),
        mpatches.Patch(color='#d29922', label='GC < 40%'),
    ]
    ax.legend(handles=legend_patches, loc='upper right', fontsize=9)

    plt.tight_layout()
    plt.savefig(out, dpi=150)
    plt.close()
    print(f'[OK] GC bar chart → {out}')


def plot_similarity_heatmap(sim_df: pd.DataFrame, out: str = 'similarity_heatmap.png', max_n: int = 20):
    # Batasi supaya heatmap tetap terbaca
    if sim_df.shape[0] > max_n:
        print(f'[WARN] Heatmap dibatasi {max_n} sekuens pertama (dari {sim_df.shape[0]})')
        sim_df = sim_df.iloc[:max_n, :max_n]

    n = sim_df.shape[0]
    fig_size = max(8, n * 0.5)
    fig, ax = plt.subplots(figsize=(fig_size, fig_size * 0.85))

    sns.heatmap(
        sim_df.astype(float),
        ax=ax,
        cmap='YlGn',
        vmin=0, vmax=100,
        annot=n <= 15,          # tampilkan angka hanya jika ≤15 sekuens
        fmt='.1f',
        linewidths=0.3,
        linecolor='#1a1a1a',
        cbar_kws={'label': 'Identity (%)'},
        square=True,
    )
    ax.set_title('Pairwise Sequence Identity (%)', fontsize=13, fontweight='bold', pad=12)
    ax.tick_params(axis='x', rotation=90, labelsize=7)
    ax.tick_params(axis='y', rotation=0, labelsize=7)

    plt.tight_layout()
    plt.savefig(out, dpi=150)
    plt.close()
    print(f'[OK] Similarity heatmap → {out}')


def plot_motif_hits(motif_results: dict, out: str = 'motif_hits.png'):
    """Bar chart total hits per motif across all sequences."""
    labels = list(motif_results.keys())
    totals = [df['Total_Hits'].sum() for df in motif_results.values()]

    fig, ax = plt.subplots(figsize=(max(8, len(labels) * 1.2), 5))
    bars = ax.bar(labels, totals, color='#58a6ff', edgecolor='none', width=0.6)

    for bar, val in zip(bars, totals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                str(val), ha='center', va='bottom', fontsize=10, fontweight='bold')

    ax.set_ylabel('Total Hits (semua sekuens)', fontsize=11)
    ax.set_xlabel('Motif / Restriction Site', fontsize=11)
    ax.set_title('Frekuensi Motif di Seluruh Sekuens', fontsize=13, fontweight='bold')
    ax.grid(axis='y', linestyle='--', alpha=0.4)
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    plt.close()
    print(f'[OK] Motif hit chart → {out}')


# ─────────────────────────────────────────────
# 6. EKSPOR CSV
# ─────────────────────────────────────────────
def export_csv(df: pd.DataFrame, path: str):
    df.to_csv(path, index=False)
    print(f'[OK] CSV → {path}')


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description='rbcL Advanced Analysis Pipeline')
    parser.add_argument('--fasta',  default='sequence_rbcl.fasta',
                        help='Path ke file FASTA (default: sequence_rbcl.fasta)')
    parser.add_argument('--motifs', nargs='*', default=None,
                        help='Motif tambahan selain preset (e.g. --motifs GAATTC GGCC)')
    parser.add_argument('--no-heatmap', action='store_true',
                        help='Skip perhitungan similarity matrix (hemat waktu untuk >100 sekuens)')
    args = parser.parse_args()

    fasta_path = args.fasta
    if not Path(fasta_path).exists():
        sys.exit(f'[ERROR] File tidak ditemukan: {fasta_path}')

    print('\n' + '='*52)
    print('  rbcL Sequence Analysis Pipeline — Advanced')
    print('='*52)

    # 1. Load
    records = load_fasta(fasta_path)

    # 2. GC Analysis
    print('\n[1/4] Analisis GC Content...')
    gc_df = analyze_gc(records)
    print('\n--- TOP 3 GC CONTENT ---')
    print(gc_df[['ID', 'GC_Content']].head(3).to_string(index=False))
    export_csv(gc_df, 'hasil_analisis_rbcl.csv')
    plot_gc_bar(gc_df)

    # 3. Similarity Matrix
    if not args.no_heatmap:
        print('\n[2/4] Menghitung Pairwise Similarity Matrix...')
        sim_df = build_similarity_matrix(records)
        export_csv(sim_df.reset_index().rename(columns={'index': 'ID'}), 'similarity_matrix_rbcl.csv')
        plot_similarity_heatmap(sim_df)
    else:
        print('\n[2/4] Similarity matrix di-skip (--no-heatmap)')

    # 4. Motif Search
    print('\n[3/4] Motif & Restriction Site Search...')
    motifs_to_search = dict(PRESET_MOTIFS)
    if args.motifs:
        for m in args.motifs:
            motifs_to_search[m] = m

    motif_results = {}
    all_motif_rows = []
    for name, motif in motifs_to_search.items():
        df_m = search_motif(records, motif, both_strands=True)
        motif_results[name] = df_m
        df_m.insert(0, 'Motif_Name', name)
        all_motif_rows.append(df_m)
        total = df_m['Total_Hits'].sum()
        seqs_found = (df_m['Total_Hits'] > 0).sum()
        print(f'  {name:12s} ({motif}): {total:4d} hits di {seqs_found}/{len(records)} sekuens')

    motif_combined = pd.concat(all_motif_rows, ignore_index=True)
    export_csv(motif_combined, 'motif_results_rbcl.csv')
    plot_motif_hits(motif_results)

    # 5. Summary
    print('\n[4/4] Summary')
    print(f'  Total sekuens  : {len(records)}')
    print(f'  Avg GC Content : {gc_df["GC_Content"].mean():.2f}%')
    print(f'  Min GC         : {gc_df["GC_Content"].min():.2f}% ({gc_df.iloc[-1]["ID"]})')
    print(f'  Max GC         : {gc_df["GC_Content"].max():.2f}% ({gc_df.iloc[0]["ID"]})')

    print('\n' + '='*52)
    print('  Output files:')
    print('    hasil_analisis_rbcl.csv')
    if not args.no_heatmap:
        print('    similarity_matrix_rbcl.csv')
        print('    similarity_heatmap.png')
    print('    motif_results_rbcl.csv')
    print('    gc_content_rbcl.png')
    print('    motif_hits.png')
    print('='*52 + '\n')


if __name__ == '__main__':
    main()
