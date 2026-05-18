import os
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Ensure output directory exists
os.makedirs("results", exist_ok=True)

print("1. Loading and Wrangling TCGA-PRAD Data...")

# Load counts (Genes are rows, Samples are columns)
# Actual file is data/TCGA-PRAD.star_counts.tsv
print("Loading expression count matrix (TCGA-PRAD.star_counts.tsv)...")
counts_df = pd.read_csv("data/TCGA-PRAD.star_counts.tsv", sep="\t", index_col=0)

# REVERSE THE LOG TRANSFORMATION (The Bioinformatic Flex)
# Xena provides log2(count + 1). We must revert to raw integer counts for DESeq2.
print("Reversing log2(count + 1) transformation to raw integer counts...")
counts_raw = (2 ** counts_df) - 1
counts_raw = counts_raw.round().astype(int)

# Transpose so Samples are rows and Genes are columns (PyDESeq2 requirement)
counts_raw = counts_raw.T 

# Load Clinical Phenotype Metadata
# Actual file is data/TCGA-PRAD.clinical.tsv
print("Loading clinical metadata (TCGA-PRAD.clinical.tsv)...")
pheno_df = pd.read_csv("data/TCGA-PRAD.clinical.tsv", sep="\t", index_col=0)

# Filter for only Primary Tumors and Solid Tissue Normal
target_conditions = ["Primary Tumor", "Solid Tissue Normal"]
pheno_filtered = pheno_df[pheno_df['sample_type.samples'].isin(target_conditions)]

# Align counts and clinical indices
print("Aligning counts and clinical indices...")
common_samples = counts_raw.index.intersection(pheno_filtered.index)
counts_raw = counts_raw.loc[common_samples]
pheno_filtered = pheno_filtered.loc[common_samples]

print(f"Aligned dataset contains {len(common_samples)} samples.")

# Filter out low-expression genes to speed up computation and increase power
print("Filtering low-expression genes (total counts across all samples >= 10)...")
genes_to_keep = counts_raw.columns[counts_raw.sum(axis=0) >= 10]
counts_raw = counts_raw[genes_to_keep]
print(f"Kept {len(genes_to_keep)} genes out of {counts_df.shape[0]}.")

print("\n2. Initializing DeseqDataSet...")
# Initialize DeseqDataSet
dds = DeseqDataSet(
    counts=counts_raw,
    metadata=pheno_filtered,
    design_factors="sample_type.samples",
    refit_cooks=True,
    n_cpus=8
)

print("Fitting size factors, dispersions, and log-fold changes...")
dds.deseq2()

print("\n3. Calculating Stats (Primary Tumor vs Solid Tissue Normal)...")
# Determine contrast factor name dynamically based on PyDESeq2 version
# In v0.4.x, factor names containing underscores are normalized to hyphens.
# In v0.5.x+, original names with underscores are preserved.
import pydeseq2
contrast_factor = (
    "sample-type.samples" if pydeseq2.__version__.startswith("0.4") else "sample_type.samples"
)

# Run Wald test statistics
stat_res = DeseqStats(
    dds,
    contrast=[contrast_factor, "Primary Tumor", "Solid Tissue Normal"]
)
stat_res.summary()

# Access results dataframe
res_df = stat_res.results_df.copy()

# Save full stats table to results/
res_df.to_csv("results/tcga_prad_differentials.csv")
print("Results successfully saved to results/tcga_prad_differentials.csv")

print("\n4. Generating Premium Volcano Plot...")

# Handle NaNs and calculate -log10(padj)
res_df = res_df.dropna(subset=['padj', 'log2FoldChange'])
res_df['-log10_padj'] = -np.log10(res_df['padj'])

# Map Ensembl IDs to Gene Symbols via MyGene.info API for professional annotation
def map_ensembl_to_symbols(ensembl_ids):
    clean_ids = [eid.split('.')[0] for eid in ensembl_ids]
    url = "https://mygene.info/v3/query"
    payload = {
        "q": clean_ids,
        "scopes": "ensembl.gene",
        "fields": "symbol",
        "species": "human"
    }
    try:
        response = requests.post(url, json=payload, timeout=10)
        if response.status_code == 200:
            results = response.json()
            mapping = {}
            for item in results:
                query_id = item.get("query")
                symbol = item.get("symbol")
                if query_id and symbol:
                    # Match clean query to original Ensembl ID
                    for original_id in ensembl_ids:
                        if original_id.startswith(query_id):
                            mapping[original_id] = symbol
            return mapping
    except Exception as e:
        print(f"Warning: MyGene.info mapping failed: {e}")
    return {}

# Identify top differentially expressed genes for annotation
sig_genes = res_df[res_df['padj'] < 0.05]
top_up = sig_genes[sig_genes['log2FoldChange'] > 1.5].nsmallest(5, 'padj')
top_down = sig_genes[sig_genes['log2FoldChange'] < -1.5].nsmallest(5, 'padj')
top_genes = pd.concat([top_up, top_down])

print("Querying MyGene.info to map Ensembl IDs to Gene Symbols...")
symbol_mapping = map_ensembl_to_symbols(top_genes.index.tolist())

# Setup plot aesthetics
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],
    'figure.dpi': 300,
    'axes.edgecolor': '#2B2B2B',
    'axes.linewidth': 1.2,
    'xtick.color': '#2B2B2B',
    'ytick.color': '#2B2B2B',
})

fig, ax = plt.subplots(figsize=(10, 8))

# Define color categories
res_df['status'] = 'Non-Significant'
res_df.loc[(res_df['padj'] < 0.05) & (res_df['log2FoldChange'] > 1), 'status'] = 'Up-regulated (FC > 2)'
res_df.loc[(res_df['padj'] < 0.05) & (res_df['log2FoldChange'] < -1), 'status'] = 'Down-regulated (FC < 0.5)'

colors = {
    'Up-regulated (FC > 2)': '#E06666',      # Sleek soft red
    'Down-regulated (FC < 0.5)': '#4C6EF5',  # Premium slate blue
    'Non-Significant': '#ADB5BD'             # Neutral grey
}

# Plot scatter points
sns.scatterplot(
    data=res_df,
    x='log2FoldChange',
    y='-log10_padj',
    hue='status',
    palette=colors,
    alpha=0.7,
    edgecolor=None,
    s=15,
    ax=ax
)

# Reference lines for thresholds
plt.axhline(-np.log10(0.05), color='#495057', linestyle='--', linewidth=1, alpha=0.7)
plt.axvline(1, color='#495057', linestyle='--', linewidth=1, alpha=0.7)
plt.axvline(-1, color='#495057', linestyle='--', linewidth=1, alpha=0.7)

# Annotate top genes
from matplotlib.patheffects import withStroke
for idx, row in top_genes.iterrows():
    symbol = symbol_mapping.get(idx, idx.split('.')[0])  # fallback to stripped Ensembl ID
    x = row['log2FoldChange']
    y = row['-log10_padj']
    
    # Elegant clean text offset
    offset_x = 0.15 if x > 0 else -0.15
    align = 'left' if x > 0 else 'right'
    
    txt = ax.text(
        x + offset_x, y, symbol,
        fontsize=9, weight='bold', color='#1A1A1A',
        horizontalalignment=align, verticalalignment='center'
    )
    txt.set_path_effects([withStroke(linewidth=3, foreground='white')])

# Visual styling details
ax.set_title("TCGA-PRAD Differential Expression: Primary Tumor vs Normal", fontsize=14, weight='bold', pad=15)
ax.set_xlabel("Log2 Fold Change", fontsize=11, weight='semibold', labelpad=8)
ax.set_ylabel("-Log10 Adjusted P-Value", fontsize=11, weight='semibold', labelpad=8)
ax.legend(title="Expression Status", title_fontsize='11', loc='upper right', frameon=True, facecolor='white', edgecolor='#E0E0E0')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()

# Save output plot
plt.savefig("results/tcga_prad_volcano.png", bbox_inches='tight')
print("Premium publication-ready Volcano Plot saved to results/tcga_prad_volcano.png")
