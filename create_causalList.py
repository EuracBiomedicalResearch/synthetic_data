#!/usr/bin/env python3

"""
GWAS Catalog to Causal Variant List Parser
==========================================

This script parses GWAS Catalog data and converts it to a format suitable for 
genetic simulation tools that require user-specified variant effects.

Input: GWAS Catalog TSV file with standard columns
Output: Causal variant list with SNP IDs and effect sizes for homozygous/heterozygous genotypes

Usage: python create_causalList.py --input gwas_catalog.txt --output causal_variants.txt
"""

import pandas as pd
import numpy as np
import argparse
import re
import logging
from pathlib import Path
import sys
import requests

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def parse_effect_size(or_beta_str, ci_str=None):
    """
    Parse effect size from OR or BETA column and confidence interval.
    
    Args:
        or_beta_str (str): The OR or BETA value from GWAS catalog
        ci_str (str): The 95% CI text (optional)
    
    Returns:
        float: Parsed effect size (converted to beta scale)
    """
    if pd.isna(or_beta_str) or or_beta_str == '' or or_beta_str == 'NR':
        return None
    
    try:
        # Clean the string
        effect_str = str(or_beta_str).strip()
        
        # Extract numeric value using regex
        # Handle formats like: "1.23", "1.23 (unit)", "0.2 unit increase", etc.
        numeric_match = re.search(r'([0-9]+\.?[0-9]*)', effect_str)
        
        if numeric_match:
            return float(numeric_match.group(1))            
        else:
            logger.warning(f"Could not parse effect size: {effect_str}")
            return None
            
    except (ValueError, TypeError) as e:
        logger.warning(f"Error parsing effect size '{or_beta_str}': {e}")
        return None

def calculate_genotype_effects(beta, model='additive', het_ratio=0.5):
    """
    Calculate genotype-specific effects from a single beta coefficient.
    
    Args:
        beta (float): Main effect size (typically for heterozygous effect or additive model)
        model (str): Genetic model ('additive', 'dominant', 'recessive')
        het_ratio (float): Ratio for heterozygous effect (0.5 = additive, 1.0 = dominant)
    
    Returns:
        tuple: (homozygous_effect, heterozygous_effect)
    """
    if beta is None or beta == 0:
        return 0, 0
    
    if model == 'additive':
        # Additive model: hetero = beta/2, homo = beta
        het_effect = beta * het_ratio
        homo_effect = beta
    elif model == 'dominant':
        # Dominant model: hetero = homo = beta
        het_effect = beta
        homo_effect = beta
    elif model == 'recessive':
        # Recessive model: hetero = 0, homo = beta
        het_effect = 0
        homo_effect = beta
    else:
        # Default to additive
        het_effect = beta * het_ratio
        homo_effect = beta
    
    return homo_effect, het_effect

def extract_snp_id(snp_str, snp_current_str=None):
    """
    Extract the primary SNP ID from various SNP string formats.
    
    Args:
        snp_str (str): SNP string from STRONGEST SNP-RISK ALLELE column
        snp_current_str (str): Current SNP ID if available
    
    Returns:
        str: Clean SNP ID (rs number)
    """
    if pd.isna(snp_str) and pd.isna(snp_current_str):
        return None
    
    # Prefer current SNP ID if available
    if not pd.isna(snp_current_str) and snp_current_str != '':
        snp_id = str(snp_current_str).strip()
        if snp_id.startswith('rs'):
            return snp_id
    
    # Parse from strongest SNP-risk allele column
    if not pd.isna(snp_str):
        snp_str = str(snp_str).strip()
        # Extract rs number (handles formats like "rs123456-A", "rs123456", etc.)
        rs_match = re.search(r'(rs\d+)', snp_str)
        if rs_match:
            return rs_match.group(1)
    
    return None


def rsids_to_variantids(rsids, batch_size=200):
    """
    Convert rsIDs to variant IDs using Ensembl REST API batch queries.
    
    Args:
        rsids (list): List of rsIDs (e.g. ["rs7412", "rs429358"])
        batch_size (int): Max number of rsIDs per API call (Ensembl allows up to 1000)
    
    Returns:
        dict: Mapping {rsid: variant_id or None}
    """
    url = "https://rest.ensembl.org/variation/human"
    headers = {"Content-Type": "application/json"}
    results = {}

    for i in range(0, len(rsids), batch_size):
        batch = rsids[i:i + batch_size]
        response = requests.post(url, headers=headers, json={"ids": batch})
        
        if not response.ok:
            raise Exception(f"Failed to fetch data for batch {i//batch_size}: {response.text}")
        
        data = response.json()
        for rsid, info in data.items():
            if "mappings" not in info or not info["mappings"]:
                results[rsid] = None
                continue
            mapping = info["mappings"][0]
            chrom = mapping["seq_region_name"]
            pos = mapping["start"]
            ref = mapping["allele_string"].split("/")[0]
            alt = mapping["allele_string"].split("/")[-1]
            results[rsid] = f"chr{chrom}:{pos}:{ref}:{alt}"
    
    return results


def process_gwas_catalog(input_file, output_file, 
                        genetic_model='additive', het_ratio=0.5):
    """
    Main processing function to convert GWAS Catalog to causal variant list.
    """
    logger.info(f"Reading GWAS Catalog file: {input_file}")
    
    try:
        df = pd.read_csv(input_file, sep='\t', low_memory=False)
        logger.info(f"Loaded {len(df)} entries from GWAS catalog")
        
        df = df[df['P-VALUE'] <= 1E-10]
        causal_variants = []
        processed_snps = set()
        
        for idx, row in df.iterrows():
            snp_id = extract_snp_id(row.get('STRONGEST SNP-RISK ALLELE'), 
                                    row.get('SNP_ID_CURRENT'))
            if not snp_id or snp_id in processed_snps:
                continue

            effect_size = parse_effect_size(row.get('OR or BETA'), 
                                            row.get('95% CI (TEXT)'))
            if effect_size is None:
                effect_size = 0
            
            homo_effect, het_effect = calculate_genotype_effects(
                effect_size, genetic_model, het_ratio
            )
            
            causal_variants.append({
                'SNP_ID': snp_id,   # temporarily keep rsID
                'HOMO_EFFECT': homo_effect,
                'HET_EFFECT': het_effect
            })
            processed_snps.add(snp_id)
        
        if not causal_variants:
            logger.error("No valid causal variants found")
            return
        
        rsids = [v['SNP_ID'] for v in causal_variants]
        logger.info(f"Converting {len(rsids)} rsIDs to variant IDs via Ensembl API...")
        rsid_to_variant = rsids_to_variantids(rsids)
        
        # Replace rsIDs with variantIDs
        for v in causal_variants:
            v['VARIANT_ID'] = rsid_to_variant.get(v['SNP_ID'], v['SNP_ID'])
        
        
        logger.info(f"Writing causal variant list to: {output_file}")
        with open(output_file, 'w') as f:
            for variant in causal_variants:
                if variant['VARIANT_ID'] is None:
                    continue  # skip unresolved ones
                f.write(f"{variant['VARIANT_ID']}\t{variant['HOMO_EFFECT']:.2f},{variant['HET_EFFECT']:.2f}\n")
        
        logger.info(f"Successfully processed {len(causal_variants)} causal variants")
    
    except Exception as e:
        logger.error(f"Error processing GWAS catalog: {e}")
        raise
    
def main():
    parser = argparse.ArgumentParser(
        description='Convert GWAS Catalog data to causal variant list for genetic simulation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python create_causalList.py --input gwas_catalog.txt --output causal_variants.txt
        """
    )
    
    parser.add_argument('--input', '-i', required=True,
                       help='Input GWAS Catalog TSV file')
    parser.add_argument('--output', '-o', required=True,
                       help='Output causal variant list file')
    parser.add_argument('--model', '-m', choices=['additive', 'dominant', 'recessive'],
                       default='additive', help='Genetic model (default: additive)')
    parser.add_argument('--het-ratio', '-r', type=float, default=0.5,
                       help='Heterozygous effect ratio for additive model (default: 0.5)')
    
    args = parser.parse_args()
    
   
    if not Path(args.input).exists():
        logger.error(f"Input file does not exist: {args.input}")
        sys.exit(1)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        process_gwas_catalog(
            input_file=args.input,
            output_file=args.output,
            genetic_model=args.model,
            het_ratio=args.het_ratio,
        )
        
        logger.info("Processing completed successfully!")
        
    except Exception as e:
        logger.error(f"Processing failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()