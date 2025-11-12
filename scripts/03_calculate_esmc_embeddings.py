#!/usr/bin/env python3
"""
AQSE Workflow - Step 3: Calculate ESM-C Embeddings

This script calculates ESM-C embeddings for:
1. All avoidome proteins
2. All similar proteins found in the similarity search

Cache Location:
    Embeddings are saved as pickle files in the cache directory specified in config.yaml:
    Format: {uniprot_id}_descriptors.pkl (e.g., P05177_descriptors.pkl)
    Content: NumPy array with shape (960,) - sequence-level ESM-C embedding

Needs esmc environment activated:
    conda activate esmc

Usage:
    Normal mode (skips existing cache):
        python 03_calculate_esmc_embeddings.py
    
    Force recalculate all (overwrites existing cache):
        python 03_calculate_esmc_embeddings.py --force-recalculate-all

Note: Use --force-recalculate-all to fix embeddings with incorrect dimensions
      (e.g., per-residue embeddings that were incorrectly saved).


"""

import pandas as pd
import numpy as np
import pickle
import logging
import re
import requests
import argparse
import os
import sys
from pathlib import Path
from typing import Set, Dict, Optional
import yaml
import time

# Set up logging first
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Import ESM-C embedding function
# Try to import from common locations or use environment variable

# Try to find single_esmc_embeddings module
esm_module_path = None

# Check environment variable first
if os.getenv('ESM_MODULE_PATH'):
    esm_module_path = Path(os.getenv('ESM_MODULE_PATH')).expanduser().resolve()
elif os.getenv('ESM_EMBEDDINGS_DIR'):
    esm_module_path = Path(os.getenv('ESM_EMBEDDINGS_DIR')).expanduser().resolve()
else:
    # Try relative path: look for analyses/run_esm_embeddings relative to script
    script_dir = Path(__file__).parent.absolute()
    # Try going up to parent directory and looking for analyses/run_esm_embeddings
    possible_paths = [
        script_dir.parent.parent / "analyses" / "run_esm_embeddings",
        script_dir.parent / "analyses" / "run_esm_embeddings",
        script_dir / "analyses" / "run_esm_embeddings",
    ]
    for path in possible_paths:
        if path.exists() and (path / "single_esmc_embeddings.py").exists():
            esm_module_path = path
            break

if esm_module_path and esm_module_path.exists():
    sys.path.insert(0, str(esm_module_path))
    logger.info(f"Using ESM module path: {esm_module_path}")
else:
    logger.warning("ESM module path not found. Trying direct import...")
    logger.info("If import fails, set ESM_MODULE_PATH or ESM_EMBEDDINGS_DIR environment variable")

try:
    from single_esmc_embeddings import get_single_esmc_embedding
except ImportError:
    logger.error("Could not import get_single_esmc_embedding from single_esmc_embeddings")
    logger.error("Please ensure the esmc package is installed and single_esmc_embeddings module is accessible")
    logger.error("You can set ESM_MODULE_PATH environment variable to point to the directory containing single_esmc_embeddings.py")
    raise


def load_config(config_path: str) -> dict:
    """Load configuration from YAML file"""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def get_avoidome_proteins(avoidome_file: str) -> Set[str]:
    """
    Load avoidome protein UniProt IDs from the avoidome file
    
    Args:
        avoidome_file: Path to avoidome protein list CSV
        
    Returns:
        Set of UniProt IDs
    """
    logger.info(f"Loading avoidome proteins from {avoidome_file}")
    df = pd.read_csv(avoidome_file)
    
    # Extract UniProt IDs from the 'UniProt ID' column
    uniprot_ids = set(df['UniProt ID'].dropna().unique())
    logger.info(f"Found {len(uniprot_ids)} unique avoidome proteins")
    
    return uniprot_ids


def get_similar_proteins(similarity_file: str) -> Set[str]:
    """
    Extract similar protein IDs from similarity search summary
    
    Args:
        similarity_file: Path to similarity search summary CSV
        
    Returns:
        Set of UniProt IDs (without _WT suffix)
    """
    logger.info(f"Loading similar proteins from {similarity_file}")
    df = pd.read_csv(similarity_file)
    
    similar_proteins = set()
    
    # Parse the 'similar_proteins' column which contains entries like:
    # "P05177 (100.0%), P04799 (75.1%), ..." (UniProt IDs without _WT suffix)
    # OR "P05177_WT (100.0%), P04799_WT (75.1%), ..." (with _WT suffix - legacy format)
    for _, row in df.iterrows():
        similar_proteins_str = row.get('similar_proteins', '')
        if pd.notna(similar_proteins_str) and similar_proteins_str.strip():
            # Try to match UniProt ID format: P/Q/O followed by 5 alphanumeric characters
            # Pattern matches: "P05177 (100.0%)" or "P05177_WT (100.0%)"
            # Extract the UniProt ID part (before _WT if present, or before space)
            matches = re.findall(r'([OPQ][0-9A-Z]{5})(?:_WT)?\s*\(', similar_proteins_str)
            for protein_id in matches:
                if protein_id:
                    similar_proteins.add(protein_id)
    
    logger.info(f"Found {len(similar_proteins)} unique similar proteins")
    return similar_proteins


def get_sequence_from_uniprot(uniprot_id: str, max_retries: int = 3) -> Optional[str]:
    """
    Fetch protein sequence from UniProt API
    
    Args:
        uniprot_id: UniProt ID
        max_retries: Maximum number of retry attempts
        
    Returns:
        Protein sequence string or None if not found
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    
    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=30)
            
            if response.status_code == 200:
                # Parse FASTA format
                lines = response.text.strip().split('\n')
                if len(lines) > 1:
                    # Skip header line (first line) and join sequence lines
                    sequence = ''.join(lines[1:])
                    logger.debug(f"Retrieved sequence for {uniprot_id} from UniProt (length: {len(sequence)})")
                    return sequence
                else:
                    logger.warning(f"Empty FASTA response from UniProt for {uniprot_id}")
                    return None
            elif response.status_code == 404:
                logger.warning(f"Protein {uniprot_id} not found in UniProt")
                return None
            else:
                logger.warning(f"UniProt API returned status {response.status_code} for {uniprot_id}")
                if attempt < max_retries - 1:
                    time.sleep(2 ** attempt)  # Exponential backoff
                    continue
                return None
                
        except requests.exceptions.RequestException as e:
            logger.warning(f"Error fetching sequence from UniProt for {uniprot_id} (attempt {attempt + 1}/{max_retries}): {e}")
            if attempt < max_retries - 1:
                time.sleep(2 ** attempt)  # Exponential backoff
                continue
            return None
    
    return None


def get_protein_sequence(uniprot_id: str, sequences_df: pd.DataFrame) -> Optional[str]:
    """
    Get protein sequence from avoidome sequences or UniProt API
    
    Args:
        uniprot_id: UniProt ID
        sequences_df: DataFrame with avoidome sequences
        
    Returns:
        Protein sequence string or None if not found
    """
    # First try avoidome sequences
    avoidome_seq = sequences_df[sequences_df['uniprot_id'] == uniprot_id]
    if not avoidome_seq.empty:
        sequence = avoidome_seq.iloc[0]['sequence']
        logger.debug(f"Found sequence for {uniprot_id} in avoidome sequences")
        return sequence
    
    # If not found, try UniProt API
    logger.info(f"Fetching sequence for {uniprot_id} from UniProt API...")
    sequence = get_sequence_from_uniprot(uniprot_id)
    
    if sequence is None:
        logger.warning(f"Sequence not found for {uniprot_id}")
    
    return sequence


def calculate_and_save_embedding(uniprot_id: str, sequence: str, 
                                 cache_dir: Path, force_recalculate: bool = False) -> bool:
    """
    Calculate ESM-C embedding for a protein and save to cache
    
    Args:
        uniprot_id: UniProt ID
        sequence: Protein sequence
        cache_dir: Cache directory path
        force_recalculate: If True, recalculate even if cache exists
        
    Returns:
        True if successful, False otherwise
    """
    cache_file = cache_dir / f"{uniprot_id}_descriptors.pkl"
    
    # Skip if already cached and not forcing recalculation
    #if cache_file.exists() and not force_recalculate:
    #    logger.info(f"ESM-C descriptors already cached for {uniprot_id}, skipping")
    #    return True
    
    try:
        logger.info(f"Calculating ESM-C embedding for {uniprot_id} (sequence length: {len(sequence)})")
        
        # Calculate embedding
        embedding = get_single_esmc_embedding(
            protein_sequence=sequence,
            model_name="esmc_300m",
            device=None,  # Auto-detect
            return_per_residue=False,
            verbose=False
        )
        
        # Save to cache
        with open(cache_file, 'wb') as f:
            pickle.dump(embedding, f)
        
        logger.info(f"Saved ESM-C descriptors for {uniprot_id} to {cache_file}")
        logger.info(f"  Embedding shape: {embedding.shape}")
        
        return True
        
    except Exception as e:
        logger.error(f"Error calculating ESM-C embedding for {uniprot_id}: {e}")
        return False


def main():
    """Main function to calculate ESM-C embeddings"""
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Calculate ESM-C embeddings for avoidome and similar proteins')
    parser.add_argument('--force-recalculate-all', action='store_true',
                       help='Force recalculation of all embeddings, even if they already exist in cache')
    parser.add_argument('--config', type=str,
                       default=None,
                       help='Path to config.yaml file. If not provided, looks for config.yaml in script directory.')
    args = parser.parse_args()
    
    # Get project root directory (parent of scripts directory)
    script_dir = Path(__file__).parent.absolute()
    project_root = script_dir.parent.absolute()
    
    # Determine config file path
    if args.config:
        config_path = Path(args.config).expanduser().resolve()
    elif os.getenv('CONFIG_FILE'):
        config_path = Path(os.getenv('CONFIG_FILE')).expanduser().resolve()
    else:
        # Default: look for config.yaml in project root
        config_path = project_root / "config.yaml"
    
    # Validate config file exists
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    logger.info(f"Using config file: {config_path}")
    
    # Load configuration
    config = load_config(str(config_path))
    
    # Get paths from config
    avoidome_file = config.get('avoidome_file')
    similarity_file = config.get('similarity_file')
    sequence_file = config.get('sequence_file')
    cache_dir = Path(config.get('papyrus_cache_dir'))
    
    # Create cache directory if it doesn't exist
    cache_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Cache directory: {cache_dir}")
    
    # Load avoidome proteins
    avoidome_proteins = get_avoidome_proteins(avoidome_file)
    
    # Load similar proteins
    similar_proteins = get_similar_proteins(similarity_file)
    
    # Combine all proteins (avoidome + similar)
    all_proteins = avoidome_proteins.union(similar_proteins)
    logger.info(f"Total unique proteins to process: {len(all_proteins)}")
    
    # Load avoidome sequences
    logger.info(f"Loading sequences from {sequence_file}")
    sequences_df = pd.read_csv(sequence_file)
    
    # Process each protein
    successful = 0
    failed = 0
    skipped = 0
    already_cached = 0
    
    for i, uniprot_id in enumerate(sorted(all_proteins), 1):
        logger.info(f"\n[{i}/{len(all_proteins)}] Processing {uniprot_id}")
        
        # Check if already cached (unless forcing recalculation)
        cache_file = cache_dir / f"{uniprot_id}_descriptors.pkl"
        if cache_file.exists() and not args.force_recalculate_all:
            logger.info(f"ESM-C descriptors already cached for {uniprot_id}, skipping calculation")
            already_cached += 1
            continue
        
        # If forcing recalculation and file exists, log it
        if cache_file.exists() and args.force_recalculate_all:
            logger.info(f"Force recalculating ESM-C descriptors for {uniprot_id} (existing cache will be overwritten)")
        
        # Get sequence
        sequence = get_protein_sequence(uniprot_id, sequences_df)
        
        if sequence is None:
            logger.warning(f"Skipping {uniprot_id} - sequence not found")
            skipped += 1
            continue
        
        # Calculate and save embedding
        if calculate_and_save_embedding(uniprot_id, sequence, cache_dir, force_recalculate=args.force_recalculate_all):
            successful += 1
        else:
            failed += 1
    
    # Summary
    logger.info("\n" + "="*60)
    logger.info("ESM-C EMBEDDING CALCULATION SUMMARY")
    logger.info("="*60)
    if args.force_recalculate_all:
        logger.info("Mode: FORCE RECALCULATE ALL (existing cache files were overwritten)")
    logger.info(f"Total proteins: {len(all_proteins)}")
    logger.info(f"  - Avoidome proteins: {len(avoidome_proteins)}")
    logger.info(f"  - Similar proteins: {len(similar_proteins)}")
    logger.info(f"Already cached (skipped): {already_cached}")
    logger.info(f"Successfully calculated: {successful}")
    logger.info(f"Failed: {failed}")
    logger.info(f"Skipped (no sequence): {skipped}")
    logger.info(f"Cache directory: {cache_dir}")
    logger.info("="*60)


if __name__ == "__main__":
    main()

