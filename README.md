# AQSE 3-Class Classification Workflow (AQSE 3c)

## Overview

The AQSE 3-Class Classification Workflow trains QSAR models for avoidome proteins using a 3-class activity classification (High/Medium/Low) with similarity-based data expansion. The workflow supports two model architectures: Random Forest and Chemprop.

## Workflow Steps

The workflow consists of 4 sequential steps:

1. **Input Preparation** - Fetches and prepares protein sequences
2. **Similarity Search** - Identifies similar proteins in Papyrus database
3. **ESM-C Descriptor Calculation** - Computes protein embeddings (optional if cached)
4. **Model Training** - Trains QSAR models using scripts/AQSE_3c_clf_th.py

---

## Step 1: Input Preparation

### Purpose
Fetches protein sequences from UniProt and prepares FASTA files for similarity search.

### Script
`scripts/01_input_preparation.py`

### Inputs
- **data/avoidome_prot_list.csv**: CSV file with avoidome proteins containing:
  - `Name_2`: Protein name
  - `UniProt ID`: UniProt identifier
  - `ChEMBL_target`: ChEMBL target identifier

### Outputs
- **fasta_files/**: Directory containing individual FASTA files for each protein
  - Format: `{uniprot_id}_{protein_name}.fasta`
- **avoidome_proteins_combined.fasta**: Combined FASTA file with all sequences
- **avoidome_sequences.csv**: CSV file with sequence information
  - Columns: `uniprot_id`, `protein_name`, `sequence`, `sequence_length`
- **input_preparation_summary.txt**: Summary of preparation process
- **logs/**: Directory with processing logs

### Usage
```bash
python scripts/01_input_preparation.py
```

Or programmatically:
```python
from scripts.input_preparation import AvoidomeInputPreparation

preparer = AvoidomeInputPreparation(
    avoidome_file="/path/to/avoidome_prot_list.csv",
    output_dir="/path/to/output"
)
preparer.prepare_inputs()
```

---

## Step 2: Similarity Search

### Purpose
Performs BLAST-based sequence similarity search against Papyrus database to identify similar proteins at three similarity thresholds (high: 70%, medium: 50%, low: 30%).

### Script
`scripts/02_protein_similarity_search_papyrus_blast.py`

### Inputs
- **fasta_files/**: Directory with FASTA files from Step 1
- **Papyrus database** (version 05.7): Automatically loaded via papyrus-scripts library
- **BLAST**: Requires BLAST+ tools installed

### Outputs
- **02_similarity_search/similarity_search_summary.csv**: Main output file with similarity results
  - Columns: `query_protein`, `threshold`, `similar_proteins`, `n_similar`, `max_identity`, `min_identity`, `avg_identity`
- **02_similarity_search/blast_db/**: BLAST database files
  - `papyrus_proteins.fasta`: FASTA file of Papyrus proteins
  - `papyrus_proteins_id_mapping.txt`: Mapping between BLAST indices and Papyrus IDs
- **02_similarity_search/plots/**: Visualization plots
  - `similar_proteins_count.png`: Bar plot showing number of similar proteins per threshold

### Usage
```python
python scripts/02_protein_similarity_search_papyrus_blast.py
```

Or programmatically:
```python
from scripts.protein_similarity_search_papyrus_blast import ProteinSimilaritySearchPapyrus

searcher = ProteinSimilaritySearchPapyrus(
    input_dir="/path/to/fasta_files",
    output_dir="/path/to/02_similarity_search",
    papyrus_version='05.7'
)
searcher.run_similarity_search()
```

### Requirements
- BLAST+ tools installed and in PATH
- papyrus-scripts library
- Internet connection for initial Papyrus dataset download

---

## Step 3: ESM-C Descriptor Calculation (Optional)

### Purpose
Calculates ESM-C (Evolutionary Scale Modeling - Contact) protein embeddings for use in Model B (PCM models). This step is **optional** if descriptors are already cached.

### Requirements
- **esmc package**: Required for calculating ESM-C descriptors
**Python environment**: Use the `esmc` conda/environment for this step


### Inputs
- Protein sequences from Step 1
- ESM-C model (handled internally)

### Outputs
- **esmc_descriptors_cache/**: Directory containing cached ESM-C descriptors
  - Format: `{uniprot_id}_descriptors.pkl` for each protein
  - Contains 960-dimensional ESM-C embeddings

### Usage
ESM-C descriptors are calculated automatically during Step 4 if not found in cache. The cache directory is specified in `config.yaml`:
```yaml
papyrus_cache_dir: "/path/to/esmc_descriptors_cache"
```

If descriptors are pre-calculated and cached, Step 3 can be skipped entirely.

**Note**: When running Step 3, ensure you activate the `esmc` environment if ESM-C descriptors need to be calculated:
```bash
conda activate esmc  # or your environment manager
python scripts/03_calculate_esmc_embeddings.py
```

---

## Step 4: Model Training

### Purpose
Trains QSAR classification models using the AQSE approach with two model types:
- **Model A**: Simple QSAR for proteins without similar proteins
- **Model B**: PCM (Protein-Compound-Model) for proteins with similar proteins, using data expansion at three thresholds (high/medium/low)

### Script
`scripts/AQSE_3c_clf_th.py`

### Requirements
- **micromamba environment**: The `chemprop_env` micromamba environment is required for Step 4 (both Random Forest and Chemprop models)

### Inputs
- **config.yaml**: Configuration file with all paths and parameters
  - `avoidome_file`: Path to avoidome_prot_list.csv
  - `similarity_file`: Path to similarity_search_summary.csv (from Step 2)
  - `sequence_file`: Path to avoidome_sequences.csv (from Step 1)
  - `activity_thresholds_file`: Path to avoidome_cutoff_table.csv
  - `output_dir`: Output directory for results
  - `papyrus_cache_dir`: Path to ESM-C descriptors cache
  - `model_type`: "random_forest" or "chemprop"
  - Model-specific hyperparameters

- **data/avoidome_cutoff_table.csv**: Activity thresholds for each protein
  - Columns: `name2_entry`, `uniprot_id`, `cutoff_high`, `cutoff_medium`

- **Papyrus database**: Bioactivity data (loaded automatically via papyrus-scripts)

### Outputs

#### Main Results
- **workflow_results.csv**: Summary of all model training attempts
- **results_{model_type}/comprehensive_protein_report.csv**: Comprehensive report with all models and metrics
- **results_{model_type}/model_summary_report.json**: JSON summary with aggregated statistics

#### Per-Model Outputs (in `results_{model_type}/`)
- **models/**: Trained model files
  - Format: `{protein}_{uniprot_id}_{model_type}_{threshold}_model.pkl.gz`
- **metrics/**: Performance metrics (JSON files)
  - Format: `{protein}_{uniprot_id}_{model_type}_{threshold}_{split}_metrics.json`
  - Contains: accuracy, F1-macro, F1-weighted, precision-macro, recall-macro, confusion matrix
- **predictions/**: Test set predictions (CSV files)
  - Format: `{protein}_{uniprot_id}_{model_type}_{threshold}_{split}_predictions.csv`
- **class_distributions/**: Class distribution files (CSV + JSON summary)
  - Format: `{protein}_{uniprot_id}_{model_type}_{threshold}_{split}_distribution.csv`
- **individual_reports/**: Detailed per-protein reports (JSON + CSV)

#### Additional Outputs
- **fingerprints/**: Morgan fingerprints cache
  - `papyrus_morgan_fingerprints.parquet`: Cached fingerprints for all compounds
- **protein_list_all_unique.csv**: List of all unique proteins (avoidome + similar)
- **umap_plots/**: UMAP visualizations (if generated)
- **logs/**: Processing logs

### Usage
```bash
micromamba activate chemprop_env
python scripts/AQSE_3c_clf_th.py
```

The script automatically:
1. Loads configuration from `config.yaml`
2. Loads all input data
3. Processes each protein in the avoidome list
4. Trains Model A (if no similar proteins) or both Model A and Model B (if similar proteins exist)
5. Saves all results and generates reports

### Model Types

#### Random Forest
- Uses precomputed features: Morgan fingerprints + Physicochemical descriptors
- Model B includes ESM-C embeddings
- Configurable parameters: `rf_n_estimators`, `rf_max_depth`, `rf_class_weight`

#### Chemprop
- Extracts features on-the-fly from SMILES
- Model B includes ESM-C embeddings automatically
- Configurable architecture: `chemprop_ffn_num_layers`, `chemprop_hidden_size`, `chemprop_depth`, etc.

### Model Status Codes
- `success`: Model trained successfully
- `insufficient_data`: Less than 30 bioactivity samples
- `insufficient_features`: Insufficient samples after feature extraction
- `no_similar_proteins`: No similar proteins found at threshold
- `no_thresholds`: No activity thresholds defined

---

## Configuration File (config.yaml)

Example configuration:
```yaml
# Input files (relative to config.yaml location)
avoidome_file: "data/avoidome_prot_list.csv"
similarity_file: "02_similarity_search/similarity_search_summary.csv"
sequence_file: "avoidome_sequences.csv"
activity_thresholds_file: "data/avoidome_cutoff_table.csv"

# Output directories
output_dir: "/path/to/output"
papyrus_cache_dir: "/path/to/esmc_descriptors_cache"

# Model configuration
model_type: "random_forest"  # or "chemprop"

# RandomForest parameters
rf_n_estimators: 500
rf_max_depth: null
rf_class_weight: "balanced_subsample"
random_state: 42

# Chemprop parameters
chemprop_max_epochs: 100
n_classes: 3
chemprop_batch_size: 200
# ... additional parameters
```

---

## Complete Workflow Execution

### Sequential Execution
```bash
# Step 1: Input preparation (base environment)
conda activate aqse_base
python scripts/01_input_preparation.py

# Step 2: Similarity search (base environment)
conda activate aqse_base
python scripts/02_protein_similarity_search_papyrus_blast.py

# Step 3: ESM-C descriptors (optional, auto-calculated if missing)
# Skip if descriptors are pre-cached
# If calculating descriptors, activate esmc environment:
conda activate esmc
python scripts/03_calculate_esmc_embeddings.py

# Step 4: Model training (chemprop_env required for both Random Forest and Chemprop)
micromamba activate chemprop_env
python scripts/AQSE_3c_clf_th.py


```

### Prerequisites

#### General Requirements
- Python 3.9-3.12
- Required packages (see `requirements.txt` or `environment.yml`):
  - pandas, numpy
  - Bio (Biopython)
  - papyrus-scripts
  - scikit-learn (for Random Forest)
  - chemprop (for Chemprop models, optional)
  - matplotlib, seaborn (for visualizations)
  - pyyaml
  - rdkit (for molecular descriptors - install via conda-forge)
- BLAST+ tools (for Step 2)
- Internet connection (for UniProt API and Papyrus download)
