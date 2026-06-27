# 🧬 rbcL Sequence Analysis Pipeline

Pipeline analisis sekuens DNA berbasis Python untuk gen **rbcL** (*ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit*) — gen kloroplas yang umum digunakan sebagai marka DNA barcode tumbuhan.

---

## 📋 Fitur

| Modul | Deskripsi |
|---|---|
| **GC Content Analysis** | Frekuensi nukleotida (A/T/G/C), GC% per sekuens, top 3 tertinggi |
| **Pairwise Similarity Matrix** | % identity antar semua pasang sekuens, visualisasi heatmap |
| **Motif & Restriction Site Search** | IUPAC-aware regex, 7 preset enzim, both strands, custom motif |
| **Visualisasi** | Bar chart GC, similarity heatmap (seaborn), motif hit chart |
| **CSV Export** | Tiga file CSV: GC analysis, similarity matrix, motif results |
| **Dashboard HTML** | Antarmuka interaktif 3 tab — tidak butuh instalasi tambahan |

---

## 🗂️ Struktur Repository

```
rbcl-pipeline/
├── pipeline_advanced.py     # Pipeline CLI utama
├── dashboard.html           # Dashboard interaktif (buka di browser)
├── sequence_rbcl.fasta      # Contoh file FASTA input (isi sendiri dari NCBI)
├── laporan_mini_project.pdf # Laporan ringkas mini project
├── README.md
└── output/                  # Hasil analisis (dibuat otomatis saat pipeline dijalankan)
    ├── hasil_analisis_rbcl.csv
    ├── similarity_matrix_rbcl.csv
    ├── motif_results_rbcl.csv
    ├── gc_content_rbcl.png
    ├── similarity_heatmap.png
    └── motif_hits.png
```

---

## ⚙️ Instalasi

```bash
# Clone repo
git clone https://github.com/username/rbcl-pipeline.git
cd rbcl-pipeline

# Install dependensi
pip install biopython pandas matplotlib seaborn
```

**Dependensi:**

- Python >= 3.8
- [BioPython](https://biopython.org/) — parsing file FASTA
- [pandas](https://pandas.pydata.org/) — tabulasi dan ekspor data
- [matplotlib](https://matplotlib.org/) — visualisasi grafik
- [seaborn](https://seaborn.pydata.org/) — heatmap similarity matrix

---

## 🚀 Cara Penggunaan

### Pipeline CLI

```bash
# Jalankan dengan file default (sequence_rbcl.fasta)
python pipeline_advanced.py

# Tentukan file FASTA sendiri
python pipeline_advanced.py --fasta data/sekuens_saya.fasta

# Tambah motif custom selain preset
python pipeline_advanced.py --motifs GAATTC GGCC TATAAA

# Skip heatmap (hemat waktu untuk >100 sekuens)
python pipeline_advanced.py --no-heatmap
```

**Argumen CLI:**

| Argumen | Default | Keterangan |
|---|---|---|
| `--fasta` | `sequence_rbcl.fasta` | Path ke file FASTA input |
| `--motifs` | *(opsional)* | Satu atau lebih motif tambahan |
| `--no-heatmap` | *(flag)* | Lewati perhitungan similarity matrix |

### Dashboard HTML

Buka file `dashboard.html` di browser (Chrome/Firefox/Edge), lalu:

1. **Drop atau klik** untuk upload file FASTA
2. Tab **GC Analysis** — lihat statistik dan bar chart
3. Tab **Similarity Heatmap** — hover cell untuk detail % identity, klik export untuk CSV
4. Tab **Motif Finder** — pilih preset restriction site atau ketik motif custom (IUPAC supported), pilih strand, klik Search

---

## 📊 Output

### 1. GC Content Analysis (`hasil_analisis_rbcl.csv`)

```
ID,Panjang (bp),A,T,G,C,GC_Content
NC_001319.1,1428,312,298,421,397,57.42
...
```

### 2. Similarity Matrix (`similarity_matrix_rbcl.csv`)

Matrix n×n pairwise identity (%). Diagonal = 100%.

```
ID,Seq1,Seq2,Seq3,...
Seq1,100.0,94.3,88.1,...
Seq2,94.3,100.0,90.2,...
```

### 3. Motif Results (`motif_results_rbcl.csv`)

```
Motif_Name,ID,Motif,Total_Hits,Forward_Hits,Reverse_Hits,Forward_Positions,Reverse_Positions
EcoRI,NC_001319.1,GAATTC,2,1,1,342,987
...
```

---

## 🔍 Preset Restriction Sites

| Nama | Motif | Keterangan |
|---|---|---|
| EcoRI | `GAATTC` | Enzim restriksi paling umum |
| HindIII | `AAGCTT` | Restriksi situs HindIII |
| BamHI | `GGATCC` | Restriksi situs BamHI |
| SmaI | `CCCGGG` | Restriksi situs SmaI |
| XhoI | `CTCGAG` | Restriksi situs XhoI |
| SphI | `GCATGC` | Restriksi situs SphI |
| TATA_box | `TATAAA` | Elemen promoter TATA box |

> Semua motif mendukung kode ambiguitas IUPAC: `R=[AG]`, `Y=[CT]`, `S=[GC]`, `W=[AT]`, `K=[GT]`, `M=[AC]`, `N=[ACGT]`

---

## 📈 Contoh Output Visualisasi

| GC Content Bar Chart | Similarity Heatmap | Motif Hit Chart |
|---|---|---|
| Batang horizontal, warna hijau (GC≥50%), biru (40–50%), amber (<40%) | Heatmap YlGn, annotasi % identity | Bar chart total hits tiap restriction site |

---

## 🗒️ Catatan Metodologi

**Pairwise Identity** yang digunakan adalah *simple positional identity*:

```
Identity (%) = jumlah posisi identik / min(panjang_A, panjang_B) × 100
```

Metode ini cocok untuk sekuens dengan panjang serupa (seperti rbcL ~1.4 kb). Untuk analisis filogenetik yang lebih akurat, disarankan menggunakan *multiple sequence alignment* (MSA) seperti MUSCLE atau MAFFT terlebih dahulu.

---

## 📁 Mendapatkan Data FASTA dari NCBI

1. Buka [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nucleotide/)
2. Cari `rbcL[Gene Name] AND plants[Organism]`
3. Pilih beberapa entri → **Send to** → **File** → Format: **FASTA**
4. Simpan sebagai `sequence_rbcl.fasta` di folder repo

---

## 📄 Laporan

Laporan ringkas tersedia di [`laporan_mini_project_rbcl.pdf`](laporan_mini_project_rbcl.pdf), mencakup pendahuluan, tujuan, metodologi, dan kesimpulan analisis.

---

## 👤 Author

**Mario** — Mahasiswa Bioinformatika, IPB University  
Mini Project — Eksplorasi & Visualisasi Data (BIF) · 2026
