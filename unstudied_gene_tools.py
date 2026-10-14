"""
Utility functions for analyzing unstudied genes protein structures.

This module contains reusable functions for:
- Downloading AlphaFold structure predictions and metadata
- Parsing CIF (Crystallographic Information File) structures
- Computing pLDDT (predicted Local Distance Difference Test) statistics
"""

import os
import sys
import re
import csv
import json
from pathlib import Path
import pandas as pd
import requests
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Polypeptide import is_aa


# ============================================================================
# Gene ID <-> Gene Name Mapping (loaded once at import time)
# ============================================================================

GENE_IDS_NAMES_PRODUCTS_PATH = Path(r"G:\My Drive\Crick_LMS\projects\a_unstudied_genes\gene_IDs_names_products.tsv")


def _load_gene_id_to_name_map(path):
    """Load a gene_systematic_id -> gene_name mapping from the PomBase gene IDs/names file."""
    try:
        gene_names_df = pd.read_csv(path, sep="\t")
        return dict(zip(gene_names_df["gene_systematic_id"], gene_names_df["gene_name"]))
    except Exception as e:
        print(f"Warning: could not load gene ID/name map from {path}: {e}", file=sys.stderr)
        return {}


# Module-level map, available as soon as unstudied_gene_tools is imported
GENE_ID_TO_NAME = _load_gene_id_to_name_map(GENE_IDS_NAMES_PRODUCTS_PATH)


def _load_gene_id_to_uniprot_map(path):
    """Load a gene_systematic_id -> uniprot_id mapping from the PomBase gene IDs/names file.

    Genes without a UniProt ID are omitted.
    """
    try:
        gene_df = pd.read_csv(path, sep="\t").dropna(subset=["uniprot_id"])
        return dict(zip(gene_df["gene_systematic_id"], gene_df["uniprot_id"]))
    except Exception as e:
        print(f"Warning: could not load gene ID/UniProt map from {path}: {e}", file=sys.stderr)
        return {}


# Module-level map of systematic ID -> UniProt ID
GENE_ID_TO_UNIPROT_ID = _load_gene_id_to_uniprot_map(GENE_IDS_NAMES_PRODUCTS_PATH)


# ============================================================================
# Ortholog Parsing Functions
# ============================================================================

def parse_orthologs(ortholog_str):
    """
    Parse a comma-separated ortholog field (e.g. from PomBase 'Human orthologs')
    into a list of gene symbols.  Returns an empty list if the value is missing
    or blank.

    Parameters
    ----------
    ortholog_str : str or float
        Raw field value, e.g. "APEH, ARSB" or NaN.

    Returns
    -------
    list of str
    """
    if pd.isna(ortholog_str) or str(ortholog_str).strip() == "":
        return []
    return [o.strip() for o in str(ortholog_str).split(",") if o.strip()]


def parse_ortholog_gene_field(gene_field):
    """
    Parse a gene field from a tab-delimited ortholog file and return a list of
    clean gene names.  Handles the notation used by PomBase ortholog exports for
    both human and S. cerevisiae orthologs.

    Rules applied
    -------------
    - "NONE"        → []
    - "A|B|C"       → ['A', 'B', 'C']     (alternative orthologs, pick each)
    - "A(C)+B(N)"   → ['A', 'B']          (gene fusion: split on +, strip (C)/(N))
    - "A(C)"        → ['A']               (single gene with terminal notation)

    Parameters
    ----------
    gene_field : str
        Raw value from the second column of the ortholog file.

    Returns
    -------
    list of str
    """
    if gene_field == "NONE":
        return []
    names = []
    for alternative in gene_field.split("|"):
        alternative = alternative.strip()
        if "+" in alternative:
            # Gene fusion: split partners and strip (C)/(N) terminal notation
            for partner in alternative.split("+"):
                clean = re.sub(r'\([CN]\)$', '', partner.strip())
                if clean:
                    names.append(clean)
        else:
            # Single name — strip (C)/(N) if present
            clean = re.sub(r'\([CN]\)$', '', alternative)
            if clean:
                names.append(clean)
    return names


# ============================================================================
# AlphaFold Download Functions
# ============================================================================

def download_alphafold_metadata(uniprot_id, systematic_id, gene_name, model_metadata_dir):
    """
    Download metadata for a protein structure from the AlphaFold DB API.
    
    Parameters
    ----------
    uniprot_id : str
        The UniProt ID of the protein
    systematic_id : str
        The systematic ID (e.g., SPAC ID for S. pombe)
    gene_name : str
        The gene name
    model_metadata_dir : str
        Directory to save the metadata JSON file
        
    Returns
    -------
    None
        Saves metadata JSON file to disk
    """
    # AlphaFold DB API URL for metadata
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
        
    # Output file path
    output_file = os.path.join(model_metadata_dir, f"metadata_{systematic_id}_{gene_name}_{uniprot_id}.json")
    
    try:
        print(f"Downloading metadata for {uniprot_id} ({systematic_id} - {gene_name})...")
        response = requests.get(url)
        response.raise_for_status()  # Raise an error for bad status codes
        
        # Convert response to JSON format
        metadata = response.json()
        
        # Save the JSON file
        with open(output_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        print(f"  ✓ Saved metadata to {output_file}")
        
    except requests.exceptions.RequestException as e:
        print(f"  ✗ Failed to download metadata for {uniprot_id}: {e}")
    except json.JSONDecodeError as e:
        print(f"  ✗ Failed to decode JSON for {uniprot_id}: {e}")

    return


def download_pdb(uniprot_id, systematic_id, gene_name, model_metadata_dir, pdb_dir):
    """
    Download PDB structure file from AlphaFold DB using metadata.
    
    Parameters
    ----------
    uniprot_id : str
        The UniProt ID of the protein
    systematic_id : str
        The systematic ID (e.g., SPAC ID for S. pombe)
    gene_name : str
        The gene name
    model_metadata_dir : str
        Directory containing the metadata JSON file
    pdb_dir : str
        Directory to save the PDB structure file
        
    Returns
    -------
    None
        Saves PDB structure file to disk
    """
    output_file = os.path.join(pdb_dir, f"{systematic_id}_{gene_name}_{uniprot_id}.pdb")

    metadata = json.load(open(os.path.join(model_metadata_dir, f"metadata_{systematic_id}_{gene_name}_{uniprot_id}.json")))
    pdb_url = metadata[0]['pdbUrl']

    try:
        print(f"Downloading pdb for {uniprot_id} ({systematic_id} - {gene_name})...")
        response = requests.get(pdb_url)
        response.raise_for_status()  # Raise an error for bad status codes
        
        pdb_data = response.text
        
        # Save the PDB file
        with open(output_file, 'w') as f:
            f.write(pdb_data)
        print(f"  ✓ Saved PDB to {output_file}")
        
    except requests.exceptions.RequestException as e:
        print(f"  ✗ Failed to download PDB for {uniprot_id}: {e}")


# ============================================================================
# CIF File Parsing Functions
# ============================================================================

def parse_cif_atom_site(cif_file_path):
    """
    Parse atom site data from a CIF file and return as a pandas DataFrame.
    
    Parameters
    ----------
    cif_file_path : str
        Path to the CIF file
        
    Returns
    -------
    pd.DataFrame
        DataFrame with atom site information, columns named after the field definitions
    """
    
    with open(cif_file_path, 'r') as f:
        lines = f.readlines()
    
    # Find the atom_site loop section
    atom_site_start = None
    for i, line in enumerate(lines):
        if line.strip().startswith('loop_') and i + 1 < len(lines):
            # Check if the next line starts with _atom_site
            if lines[i + 1].strip().startswith('_atom_site.'):
                atom_site_start = i
                break
    
    if atom_site_start is None:
        raise ValueError("Could not find _atom_site loop section in CIF file")
    
    # Extract column headers
    headers = []
    data_start = None
    
    for i in range(atom_site_start + 1, len(lines)):
        line = lines[i].strip()
        
        # Column headers start with _atom_site.
        if line.startswith('_atom_site.'):
            # Extract the part after the dot
            header = line.split('_atom_site.')[1]
            headers.append(header)
        
        # Data lines start with ATOM or HETATM (not starting with _)
        elif line.startswith('ATOM') or line.startswith('HETATM'):
            data_start = i
            break
        
        # Stop if we hit another section
        elif line.startswith('#') and i > atom_site_start + 2:
            break
    
    if data_start is None:
        raise ValueError("Could not find data lines in _atom_site section")
    
    # Parse data rows
    data_rows = []
    for i in range(data_start, len(lines)):
        line = lines[i].strip()
        
        # Stop at the next section (empty line or # comment)
        if not line or line.startswith('#'):
            break
        
        # Split by whitespace, but be careful with quoted fields
        # CIF format uses spaces as delimiters for unquoted fields
        values = []
        current_field = ""
        in_quotes = False
        
        for char in line:
            if char == '"':
                in_quotes = not in_quotes
            elif char == ' ' and not in_quotes:
                if current_field:
                    values.append(current_field)
                    current_field = ""
            else:
                current_field += char
        
        if current_field:
            values.append(current_field)
        
        if len(values) == len(headers):
            data_rows.append(values)
    
    # Create DataFrame
    df = pd.DataFrame(data_rows, columns=headers)
    
    return df


# ============================================================================
# pLDDT Statistics Functions
# ============================================================================

def calculate_plddt_stats(df, complex_name):
    """
    Calculate average pLDDT (predicted Local Distance Difference Test) scores 
    for each protein and the whole complex.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing atom site data with 'auth_asym_id' and 'B_iso_or_equiv' columns
    complex_name : str
        Name of the complex for display in output
        
    Returns
    -------
    tuple
        (stats_df, overall_avg_plddt)
        - stats_df: DataFrame with per-protein statistics
        - overall_avg_plddt: Overall average pLDDT for the complex
    """
    
    # Convert B_iso_or_equiv to numeric (in case there are any non-numeric values)
    df['B_iso_or_equiv_numeric'] = pd.to_numeric(df['B_iso_or_equiv'], errors='coerce')
    
    # Get unique proteins
    proteins = df['auth_asym_id'].unique()
    
    print(f"\n{'='*60}")
    print(f"pLDDT Statistics for {complex_name}")
    print(f"{'='*60}")
    
    # Calculate per-protein average
    per_protein_stats = []
    for protein in sorted(proteins):
        protein_df = df[df['auth_asym_id'] == protein]
        avg_plddt = protein_df['B_iso_or_equiv_numeric'].mean()
        n_atoms = len(protein_df)
        per_protein_stats.append({
            'Protein': protein,
            'Atoms': n_atoms,
            'Avg_pLDDT': avg_plddt
        })
        print(f"\nProtein {protein}:")
        print(f"  Number of atoms: {n_atoms}")
        print(f"  Average pLDDT:   {avg_plddt:.2f}")
    
    # Calculate overall average
    overall_avg_plddt = df['B_iso_or_equiv_numeric'].mean()
    print(f"\nOverall Complex:")
    print(f"  Total atoms: {len(df)}")
    print(f"  Average pLDDT: {overall_avg_plddt:.2f}")
    
    stats_df = pd.DataFrame(per_protein_stats)
    stats_df['Complex'] = complex_name
    
    return stats_df, overall_avg_plddt


# ============================================================================
# AlphaFold Download Orchestration
# ============================================================================

def alphafold_download(input_file, output_dir='.'):
    """
    Download metadata and PDB structures from AlphaFold DB for proteins in a TSV file.
    
    Reads a TSV file with protein information, queries the AlphaFold DB API for each UniProt ID,
    and downloads both metadata (JSON) and structure files (PDB).
    
    Parameters
    ----------
    input_file : str
        Path to the tab-delimited input file with protein information.
        Required columns: "Systematic ID", "Gene name", "UniProt ID"
    output_dir : str, optional
        Directory where output files will be saved. Default is current directory.
        Creates subdirectories: "model_metadata/" for JSON files and "pdb/" for PDB files.
        
    Returns
    -------
    None
        Saves metadata and PDB files to disk in the output directory.
        
    Raises
    ------
    ValueError
        If input file does not contain required columns
    FileNotFoundError
        If input file does not exist
    """
    
    # Create output directory if it doesn't exist
    if output_dir != ".":
        try:
            os.makedirs(output_dir, exist_ok=True)
        except Exception as e:
            print(f"Error creating output directory '{output_dir}': {str(e)}", file=sys.stderr)
            return
    
    print("=" * 60)
    print("AlphaFold PDB Downloader")
    print("=" * 60)
    print(f"Reading sequences from: {input_file}")
    print(f"Output directory:       {output_dir}")
    print()
    
    # Create model_metadata directory if it doesn't exist
    model_metadata_dir = os.path.join(output_dir, "model_metadata")
    try:
        os.makedirs(model_metadata_dir, exist_ok=True)
    except Exception as e:
        print(f"Error creating model_metadata directory '{model_metadata_dir}': {str(e)}", file=sys.stderr)
        return

    # Create pdb directory
    pdb_dir = os.path.join(output_dir, "pdb")
    try:
        os.makedirs(pdb_dir, exist_ok=True)
    except Exception as e:
        print(f"Error creating pdb directory '{pdb_dir}': {str(e)}", file=sys.stderr)
        return

   
    protein_metadata = []
    with open(input_file, 'r') as f_in: 
        # Read header line
        header = f_in.readline().strip().split('\t')
        
        # Find column indices
        try:
            systematic_id_idx = header.index('Systematic ID')
            gene_name_idx = header.index('Gene name')
            uniprot_id_idx = header.index('UniProt ID')
        except ValueError as e:
            print(f"Error: Could not find required column in header: {e}", file=sys.stderr)
            return
        
        # Read data rows
        for line in f_in: 
            if line.strip():  # Skip empty lines
                linesp = line.strip().split('\t')
                systematic_id = linesp[systematic_id_idx]
                gene_name = linesp[gene_name_idx]
                uniprot_id = linesp[uniprot_id_idx]
                
                protein_metadata.append({
                    'systematic_id': systematic_id,
                    'gene_name': gene_name,
                    'uniprot_id': uniprot_id
                })

    for protein in protein_metadata:       
        download_alphafold_metadata(protein['uniprot_id'], protein['systematic_id'], protein['gene_name'], model_metadata_dir)
        download_pdb(protein['uniprot_id'], protein['systematic_id'], protein['gene_name'], model_metadata_dir, pdb_dir)
    
    print(f"\nMetadata and PDB download complete!")
    print(f"\nMetadata Files saved in: {model_metadata_dir}")
    print(f"\nPDB Files saved in: {pdb_dir}")


# ============================================================================
# PDB File Manipulation Functions (using BioPython)
# ============================================================================

def extract_sequence_from_pdb(pdb_content):
    """
    Extract amino acid sequence from PDB file content.
    
    Reads a PDB file content and extracts the protein sequence from chain A
    by converting 3-letter residue codes to 1-letter amino acid codes.
    
    Parameters
    ----------
    pdb_content : str
        PDB file content as a string
        
    Returns
    -------
    str
        Protein sequence in 1-letter amino acid code.
        Unknown residues are represented as 'X'.
        
    Examples
    --------
    >>> with open('protein.pdb', 'r') as f:
    ...     sequence = extract_sequence_from_pdb(f.read())
    >>> print(sequence)
    'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV...'
    """
    from io import StringIO
    
    # Three-letter to one-letter amino acid code mapping
    three_to_one_code = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', StringIO(pdb_content))
    
    # Get all residues from chain A (primary chain)
    residues = []
    for model in structure:
        for chain in model:
            # Only process chain A to get single sequence
            if chain.get_id() == 'A':
                for residue in chain:
                    # Get the 3-letter residue name
                    res_name = residue.get_resname().strip()
                    # Convert 3-letter to 1-letter code
                    one_letter = three_to_one_code.get(res_name, 'X')
                    residues.append(one_letter)
    
    return ''.join(residues)


def extract_pdb_subset(pdb_file, chain_label, residue_range, output_file):
    """
    Extract a subset of a PDB file by chain and residue range using BioPython.
    
    Creates a new PDB file containing only atoms from a specified chain and 
    residue number range. Uses BioPython's PDBParser for robust PDB parsing
    and PDBIO for output.
    
    Parameters
    ----------
    pdb_file : str
        Path to the input PDB file
    chain_label : str
        Chain identifier (e.g., 'A', 'B', or '*' for all chains)
    residue_range : tuple, list, or str
        Residue range to extract. Can be specified as:
        - Tuple or list: (start_residue, end_residue), e.g., (10, 50)
        - String: "start-end", e.g., "10-50"
        - None: Include all residues in the selected chain
    output_file : str
        Path to write the output PDB file
        
    Returns
    -------
    int
        Number of atoms written to the output file
        
    Raises
    ------
    FileNotFoundError
        If the input PDB file does not exist
    ValueError
        If residue_range format is invalid
        
    Examples
    --------
    >>> # Extract chain A, residues 1-100
    >>> extract_pdb_subset('input.pdb', 'A', (1, 100), 'output.pdb')
    
    >>> # Extract chain B, residues 50-150 using string format
    >>> extract_pdb_subset('input.pdb', 'B', '50-150', 'output.pdb')
    
    >>> # Extract all residues from chain A
    >>> extract_pdb_subset('input.pdb', 'A', None, 'output.pdb')
    """
    
    if not os.path.exists(pdb_file):
        raise FileNotFoundError(f"Input PDB file not found: {pdb_file}")
    
    # Parse residue range
    if residue_range is not None:
        if isinstance(residue_range, str):
            try:
                parts = residue_range.split('-')
                start_res = int(parts[0].strip())
                end_res = int(parts[1].strip())
            except (ValueError, IndexError):
                raise ValueError(f"Invalid residue range format: {residue_range}. Use 'start-end' format.")
        elif isinstance(residue_range, (tuple, list)):
            start_res, end_res = residue_range
        else:
            raise ValueError("residue_range must be a tuple (start, end), list, or string 'start-end'")
    else:
        start_res = None
        end_res = None
    
    try:
        # Parse the PDB file
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file)
        
        # Create a new structure for output
        io = PDBIO()
        io.set_structure(structure)
        
        # Define selection function for residues to include
        class ChainRangeSelect:
            def __init__(self, chain_id, start, end):
                self.chain_id = chain_id
                self.start = start
                self.end = end
            
            def accept_model(self, model):
                return True
            
            def accept_chain(self, chain):
                if self.chain_id == '*':
                    return True
                return chain.id == self.chain_id
            
            def accept_residue(self, residue):
                # Get residue number
                res_num = residue.id[1]
                
                # Check residue range
                if self.start is not None and self.end is not None:
                    return self.start <= res_num <= self.end
                else:
                    return True
            
            def accept_atom(self, atom):
                return True
        
        # Create selector
        selector = ChainRangeSelect(chain_label, start_res, end_res)
        
        # Save the subset
        io.save(output_file, selector)
        
        # Count atoms written
        atoms_written = 0
        parser_count = PDBParser(QUIET=True)
        structure_out = parser_count.get_structure('output', output_file)
        for model in structure_out:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        atoms_written += 1
        
        print(f"✓ Extracted {atoms_written} atoms from chain {chain_label}", end="")
        if residue_range is not None:
            print(f", residues {start_res}-{end_res}")
        else:
            print()
        print(f"✓ Output saved to: {output_file}")
        
        return atoms_written
    
    except Exception as e:
        print(f"✗ Error processing PDB file: {e}")
        raise


# ============================================================================
# AlphaFold3 PPI Analysis Functions
# ============================================================================

def extract_af3_ppi_metrics(filename):
    """
    Extract useful metrics from an AlphaFold3 protein-protein interaction prediction.
    
    Parses a summary_confidences JSON file from AlphaFold3 and extracts key metrics
    including protein identifiers, confidence scores, and inter-chain PAE (Predicted
    Aligned Error) values.
    
    Parameters
    ----------
    filename : str
        Path to the summary_confidences JSON file from AlphaFold3.
        Expected filename format: fold_<Protein_A>_<Protein_B>_<timestamp>_summary_confidences_<number>.json
        
    Returns
    -------
    dict
        Dictionary containing extracted metrics:
        - Protein_A, Protein_B: Protein identifiers extracted from filename
        - fraction_disordered, has_clash, iptm, num_recycles, ptm, ranking_score: 
          Confidence metrics from the JSON file
        - chain_pair_pae_min_A_A, chain_pair_pae_min_A_B, chain_pair_pae_min_B_A, 
          chain_pair_pae_min_B_B: Inter-chain PAE values from chain_pair_pae_min field
        - chain_ptm_A, chain_ptm_B: Per-chain PTM (predicted TM-score) values from chain_ptm field
        
    Raises
    ------
    FileNotFoundError
        If the JSON file does not exist
    ValueError
        If the filename format is invalid or required fields are missing from JSON
    json.JSONDecodeError
        If the JSON file is malformed
        
    Examples
    --------
    >>> metrics = extract_af3_ppi_metrics('fold_rea1_yta12_2026_04_24_11_03_summary_confidences_0.json')
    >>> print(metrics['Protein_A'])
    'rea1'
    """
    
    # Check if file exists
    if not os.path.exists(filename):
        raise FileNotFoundError(f"JSON file not found: {filename}")
    
    # Extract protein identifiers from filename
    # Expected format: fold_<Protein_A>_<Protein_B>_<timestamp>_summary_confidences_<number>.json
    basename = os.path.basename(filename)
    parts = basename.replace('.json', '').split('_')
    
    # Should start with 'fold', then have Protein_A and Protein_B
    if len(parts) < 3 or parts[0] != 'fold':
        raise ValueError(f"Invalid filename format: {basename}")
    
    # Extract protein identifiers and convert to uppercase
    protein_a = parts[1].upper()
    protein_b = parts[2].upper()
    
    # Read JSON file
    with open(filename, 'r') as f:
        json_data = json.load(f)
    
    # Extract metrics from JSON
    metrics = {
        'Protein_A': protein_a,
        'Protein_B': protein_b,
        'fraction_disordered': json_data.get('fraction_disordered'),
        'has_clash': json_data.get('has_clash'),
        'iptm': json_data.get('iptm'),
        'num_recycles': json_data.get('num_recycles'),
        'ptm': json_data.get('ptm'),
        'ranking_score': json_data.get('ranking_score')
    }
    
    # Extract chain-specific metrics if available
    if 'chain_pair_pae_min' in json_data:
        chain_pae = json_data['chain_pair_pae_min']
        # Assuming chain_pae is a 2x2 matrix or dict with appropriate keys
        metrics['chain_pair_pae_min_A_A'] = chain_pae[0][0] if isinstance(chain_pae, list) else chain_pae.get('A_A')
        metrics['chain_pair_pae_min_A_B'] = chain_pae[0][1] if isinstance(chain_pae, list) else chain_pae.get('A_B')
        metrics['chain_pair_pae_min_B_A'] = chain_pae[1][0] if isinstance(chain_pae, list) else chain_pae.get('B_A')
        metrics['chain_pair_pae_min_B_B'] = chain_pae[1][1] if isinstance(chain_pae, list) else chain_pae.get('B_B')
    
    if 'chain_ptm' in json_data:
        chain_ptm = json_data['chain_ptm']
        metrics['chain_ptm_A'] = chain_ptm[0] if isinstance(chain_ptm, list) else chain_ptm.get('A')
        metrics['chain_ptm_B'] = chain_ptm[1] if isinstance(chain_ptm, list) else chain_ptm.get('B')
    
    return metrics


def extract_msas_from_af3_data(data_directory):
    """
    Extract Multiple Sequence Alignments (MSAs) from AlphaFold3 data.json files.
    
    Reads the data.json file in the specified directory, extracts paired and unpaired
    MSAs for both proteins, and saves them as .a3m files in a new 'msas' subdirectory.
    
    Parameters
    ----------
    data_directory : str
        Path to directory containing the *_data.json file.
        Expected format: <path>/pair_<index>_<uniprot_A>_<uniprot_B>/
        
    Returns
    -------
    None
        Creates 'msas' subdirectory and saves 4 .a3m files:
        - pair_<index>_<uniprot_A>_<uniprot_B>_unpaired_msa_chains_a.a3m
        - pair_<index>_<uniprot_A>_<uniprot_B>_unpaired_msa_chains_b.a3m
        - pair_<index>_<uniprot_A>_<uniprot_B>_paired_msa_chains_a.a3m
        - pair_<index>_<uniprot_A>_<uniprot_B>_paired_msa_chains_b.a3m
        
    Raises
    ------
    FileNotFoundError
        If no *_data.json file is found in the directory
    KeyError
        If expected keys are missing from the JSON structure
    json.JSONDecodeError
        If the JSON file is malformed
        
    Examples
    --------
    >>> extract_msas_from_af3_data('/path/to/pair_0_O60082_Q6LA55')
    ✓ Extracted MSAs for pair_0_O60082_Q6LA55
    ✓ Saved: pair_0_O60082_Q6LA55_unpaired_msa_chains_a.a3m
    ✓ Saved: pair_0_O60082_Q6LA55_unpaired_msa_chains_b.a3m
    ✓ Saved: pair_0_O60082_Q6LA55_paired_msa_chains_a.a3m
    ✓ Saved: pair_0_O60082_Q6LA55_paired_msa_chains_b.a3m
    """
    
    # Find the _data.json file in the directory
    data_file = None
    for file in os.listdir(data_directory):
        if file.endswith('_data.json'):
            data_file = os.path.join(data_directory, file)
            break
    
    if data_file is None:
        raise FileNotFoundError(f"No *_data.json file found in {data_directory}")
    
    # Extract base name from directory or JSON file
    dir_basename = os.path.basename(data_directory.rstrip(os.sep))
    file_basename = os.path.basename(data_file).replace('_data.json', '')
    base_name = dir_basename if dir_basename.startswith('pair_') else file_basename
    
    # Create msas subdirectory
    msas_dir = os.path.join(data_directory, 'msas')
    os.makedirs(msas_dir, exist_ok=True)
    
    # Read JSON data
    try:
        with open(data_file, 'r') as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        raise json.JSONDecodeError(f"Failed to parse {data_file}: {e.msg}", e.doc, e.pos)
    
    # Extract and save MSAs
    protein_labels = ['a', 'b']
    msa_types = [('unpairedMsa', 'unpaired'), ('pairedMsa', 'paired')]
    
    try:
        sequences = data['sequences']
        num_proteins = len(sequences)
        
        # Check if this is a homodimer (only one unique sequence provided)
        if num_proteins == 1:
            print(f"⚠ Note: Homodimer detected - only one sequence provided")
        
        for protein_idx, protein_label in enumerate(protein_labels):
            # For homodimers, use the single sequence for both proteins
            if protein_idx >= num_proteins:
                if num_proteins == 1:
                    print(f"⚠ Warning: Using protein A MSA for protein B (homodimer)")
                    protein_data = sequences[0]['protein']
                else:
                    print(f"⚠ Warning: Protein {protein_label} not available in sequence data")
                    continue
            else:
                protein_data = sequences[protein_idx]['protein']
            
            for json_msa_type, file_msa_type in msa_types:
                # Get MSA content
                msa_content = protein_data.get(json_msa_type)
                
                if msa_content is None:
                    print(f"⚠ Warning: {json_msa_type} not found for protein {protein_label} in {base_name}")
                    continue
                
                # Create output filename
                output_filename = f"{base_name}_{file_msa_type}_msa_chains_{protein_label}.a3m"
                output_path = os.path.join(msas_dir, output_filename)
                
                # Save MSA to file
                with open(output_path, 'w') as f:
                    f.write(msa_content)
                
                print(f"✓ Saved: {output_filename}")
    
    except KeyError as e:
        raise KeyError(f"Expected key not found in JSON structure: {e}")
    
    print(f"✓ Extracted MSAs for {base_name}")


def process_alphastream_output(input_dir, output_dir):
    """
    Process Alphastream cofolding output folders into a single dataframe.

    The function scans pair directories in the input directory, extracts the
    summary confidence metrics from each seed/sample JSON file, assembles the
    rows into a dataframe, and saves the combined output as
    ``alphastream_cofolding_data.csv`` in the output directory.

    Parameters
    ----------
    input_dir : str
        Directory containing ``pair_<index>_<uniprot_A>_<uniprot_B>`` folders.
    output_dir : str
        Directory where ``alphastream_cofolding_data.csv`` will be written.

    Returns
    -------
    pd.DataFrame
        Combined Alphastream cofolding dataframe.
    """

    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    if not input_path.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    folder_pattern = re.compile(r'^pair_(\d+)_([A-Za-z0-9]+)_([A-Za-z0-9]+)$')
    subfolder_pattern = re.compile(r'^seed-1_sample-(\d+)$')

    cofolding_data = []
    processed_pairs = []

    print(f"\nScanning {input_path}...")

    for folder in sorted(input_path.iterdir()):
        if not folder.is_dir():
            continue

        folder_name = folder.name
        match = folder_pattern.match(folder_name)

        if not match:
            print(f"⚠ Skipping folder (incorrect format): {folder_name}")
            continue

        pair_index, uniprot_a, uniprot_b = match.groups()

        # Convert to uppercase for consistency across runs and comparisons.
        uniprot_a = uniprot_a.upper()
        uniprot_b = uniprot_b.upper()
        processed_pairs.append((uniprot_a, uniprot_b))

        for subfolder in sorted(folder.iterdir()):
            if not subfolder.is_dir():
                continue

            subfolder_name = subfolder.name
            subfolder_match = subfolder_pattern.match(subfolder_name)

            if not subfolder_match:
                print(f"⚠ Skipping subfolder (incorrect format): {subfolder_name}")
                continue

            run_number = int(subfolder_match.group(1))

            json_filename = (
                f"pair_{pair_index}_{uniprot_a}_{uniprot_b}_"
                f"seed-1_sample-{run_number}_summary_confidences.json"
            )
            json_path = subfolder / json_filename
            if not json_path.exists():
                print(f"⚠ Missing JSON file: {json_filename}")
                continue

            try:
                with open(json_path, 'r') as f:
                    json_data = json.load(f)

                row_data = {
                    'uniprot_ID_A': uniprot_a,
                    'uniprot_ID_B': uniprot_b,
                    'run_number': run_number,
                    'fraction_disordered': json_data.get('fraction_disordered', None),
                    'has_clash': json_data.get('has_clash', None),
                    'iptm': json_data.get('iptm', None),
                    'ptm': json_data.get('ptm', None),
                    'ranking_score': json_data.get('ranking_score', None)
                }
                cofolding_data.append(row_data)

            except Exception as e:
                print(f"✗ Error reading {json_filename}: {e}")

    cofolding_df = pd.DataFrame(cofolding_data)

    if not cofolding_df.empty:
        cofolding_df = cofolding_df.sort_values(
            ['uniprot_ID_A', 'uniprot_ID_B', 'run_number']
        ).reset_index(drop=True)

    output_csv = output_path / 'alphastream_cofolding_data.csv'
    cofolding_df.to_csv(output_csv, index=False)

    print(f"\n✓ Processed {len(set(processed_pairs))} unique protein pairs")
    print(f"✓ Extracted data from {len(cofolding_df)} JSON files")
    print(f"✓ Alphastream cofolding data saved to: {output_csv}")
    print(f"\nDataFrame shape: {cofolding_df.shape}")
    print(f"Columns: {cofolding_df.columns.tolist()}")

    return cofolding_df


# ============================================================================
# Protein-Protein Interaction (PPI) Compilation Functions
# ============================================================================

def normalize_pombe_gene_id(gene_id):
    """
    Normalize an S. pombe systematic gene ID to the standard PomBase casing,
    where the trailing strand-designation letter (e.g. 'c' for Crick strand)
    is lowercase and the rest of the identifier is left unchanged.

    Parameters
    ----------
    gene_id : str
        Systematic gene ID, e.g. "SPAC1002.03C" or "SPAC1002.03c".

    Returns
    -------
    str
        Normalized gene ID, e.g. "SPAC1002.03c".
    """
    gene_id = str(gene_id).strip()
    if gene_id and gene_id[-1].isalpha():
        return gene_id[:-1] + gene_id[-1].lower()
    return gene_id


def sort_gene_pair(gene_a, gene_b):
    """
    Sort a pair of gene names alphabetically (case-insensitive) so that the
    first element of the returned tuple always precedes the second.

    This is the project standard for ordering gene pairs in PPI tables so that
    a given interaction always appears in the same column order, regardless of
    which order it was reported in the source data. If exactly one of the two
    genes is missing (NA), the missing value is always placed second, so the
    known gene stays in the first slot.

    Parameters
    ----------
    gene_a : str
        First gene name.
    gene_b : str
        Second gene name.

    Returns
    -------
    tuple of (str, str)
        (first_gene, second_gene) sorted case-insensitively, with a lone NA
        value always placed second.
    """
    a_is_na = pd.isna(gene_a)
    b_is_na = pd.isna(gene_b)

    if a_is_na and not b_is_na:
        return gene_b, gene_a
    if b_is_na or a_is_na:
        return gene_a, gene_b

    if str(gene_a).lower() <= str(gene_b).lower():
        return gene_a, gene_b
    return gene_b, gene_a


def sort_gene_pair_columns(df, col_a, col_b, linked_column_pairs=None):
    """
    Given a dataframe with two gene-name columns, return a copy where the
    values in each row have been reordered so that ``col_a`` alphabetically
    (case-insensitive) precedes ``col_b``.

    Some datasets carry additional columns that are tied to a specific gene
    of the pair (e.g. a per-gene source ID), rather than to the interaction
    as a whole. Use ``linked_column_pairs`` to have those columns swapped in
    tandem with ``col_a``/``col_b`` whenever a row is reordered, so they stay
    correctly associated with the gene now in each slot.

    If exactly one of ``col_a``/``col_b`` is missing (NA) for a row, the
    missing value is always placed in ``col_b``, so a known gene always ends
    up in ``col_a``.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing the two gene columns.
    col_a : str
        Name of the first gene column.
    col_b : str
        Name of the second gene column.
    linked_column_pairs : list of (str, str) tuples, optional
        Additional column pairs tied to col_a/col_b (e.g.
        [("orf_name_1", "orf_name_2"), ("gene_name_1", "gene_name_2")]) that
        should be swapped together whenever col_a/col_b are swapped.

    Returns
    -------
    pd.DataFrame
        Copy of ``df`` with ``col_a``/``col_b`` (and any linked columns)
        values sorted/swapped per row.
    """
    df = df.copy()
    linked_column_pairs = linked_column_pairs or []

    a_na = df[col_a].isna()
    b_na = df[col_b].isna()

    alphabetical_swap = pd.Series(
        [str(a).lower() > str(b).lower() for a, b in zip(df[col_a], df[col_b])],
        index=df.index,
    )

    # A lone NA always belongs in col_b; ties (both/neither NA) fall back to
    # the alphabetical comparison.
    needs_swap = alphabetical_swap.where(~(a_na | b_na), False)
    needs_swap = needs_swap.mask(a_na & ~b_na, True)
    needs_swap = needs_swap.mask(b_na & ~a_na, False)

    for left_col, right_col in [(col_a, col_b)] + list(linked_column_pairs):
        left = df[left_col].copy()
        right = df[right_col].copy()
        df[left_col] = right.where(needs_swap, left)
        df[right_col] = left.where(needs_swap, right)

    return df


def annotate_unknown_and_pug(df, gene_col, unknown_df, gene_id_col='gene_id',
                              class_col='unknown_class', pug_label='Priority unstudied'):
    """
    Add boolean "unknown" and "PUG" (Priority Unstudied Gene) annotation
    columns to a dataframe, based on membership of ``gene_col`` in the
    unstudied genes reference table.

    Adds two columns named ``f"{gene_col}_unknown"`` and ``f"{gene_col}_PUG"``:
    - ``{gene_col}_unknown`` is True if the gene appears in ``unknown_df``.
    - ``{gene_col}_PUG`` is True if the gene appears in ``unknown_df`` with
      ``class_col`` equal to ``pug_label``.
    - Both are NA (rather than False) for rows where ``gene_col`` itself is NA
      (i.e. no corresponding S. pombe ortholog), since "unknown"/"PUG" status
      is undefined without a gene to look up.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing the gene column to annotate.
    gene_col : str
        Name of the column in ``df`` holding gene IDs.
    unknown_df : pd.DataFrame
        Reference dataframe of unstudied genes (e.g. the PUG dataset).
    gene_id_col : str, optional
        Column in ``unknown_df`` holding gene IDs. Default 'gene_id'.
    class_col : str, optional
        Column in ``unknown_df`` holding the unstudied classification.
        Default 'unknown_class'.
    pug_label : str, optional
        Value in ``class_col`` that designates a Priority Unstudied Gene.
        Default 'Priority unstudied'.

    Returns
    -------
    pd.DataFrame
        Copy of ``df`` with the two new nullable-boolean columns appended.
    """
    df = df.copy()
    unknown_genes = set(unknown_df[gene_id_col])
    pug_genes = set(unknown_df.loc[unknown_df[class_col] == pug_label, gene_id_col])

    gene_is_na = df[gene_col].isna()
    df[f"{gene_col}_unknown"] = df[gene_col].isin(unknown_genes).astype('boolean').mask(gene_is_na, pd.NA)
    df[f"{gene_col}_PUG"] = df[gene_col].isin(pug_genes).astype('boolean').mask(gene_is_na, pd.NA)
    return df


def map_interactions_to_ortholog(df, col_a, col_b, ortholog_map):
    """
    Map a two-column interaction dataframe from one species' gene identifiers
    to S. pombe ortholog gene IDs, expanding rows when a gene maps to more
    than one ortholog. All original columns are preserved and carried over to
    the expanded rows, so source gene names/IDs and any other interaction
    data (e.g. screen/publication columns) remain linked to the new rows.

    For each interaction (gene_a, gene_b), every combination of gene_a's
    ortholog(s) and gene_b's ortholog(s) is emitted as a separate row (e.g. a
    2-ortholog gene paired with a 3-ortholog gene yields 6 rows). If a gene
    has no entry in ``ortholog_map`` (or the id itself is missing), it is
    treated as having no ortholog and the corresponding output column is set
    to NA for that row.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with two columns holding gene identifiers for interactor A/B
        in the source species (may include other columns, which are kept).
    col_a : str
        Name of the column with the first interactor's ID.
    col_b : str
        Name of the column with the second interactor's ID.
    ortholog_map : dict or defaultdict(list)
        Mapping from source-species gene ID to a list of one or more S. pombe
        ortholog gene IDs.

    Returns
    -------
    pd.DataFrame
        Copy of ``df`` with 'sp_gene_1' and 'sp_gene_2' columns added as the
        leftmost columns, with one row per (ortholog_a, ortholog_b)
        combination.
    """
    def _orthologs(gene_id):
        if pd.isna(gene_id):
            return [pd.NA]
        orthologs = ortholog_map.get(gene_id)
        return orthologs if orthologs else [pd.NA]

    rows = []
    for record in df.to_dict(orient='records'):
        for ortholog_a in _orthologs(record[col_a]):
            for ortholog_b in _orthologs(record[col_b]):
                rows.append({'sp_gene_1': ortholog_a, 'sp_gene_2': ortholog_b, **record})

    return pd.DataFrame(rows)


def add_gene_name_column(df, gene_col, name_map=None):
    """
    Insert a companion "common name" column immediately to the right of
    ``gene_col``, mapping each S. pombe systematic gene ID to its common gene
    name. Falls back to the systematic ID itself when no name is available.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing the gene ID column to annotate.
    gene_col : str
        Name of the column in ``df`` holding systematic gene IDs.
    name_map : dict, optional
        Mapping from systematic gene ID to gene name. Defaults to the
        module-level ``GENE_ID_TO_NAME`` map (loaded from
        gene_IDs_names_products.tsv).

    Returns
    -------
    pd.DataFrame
        Copy of ``df`` with a new ``f"{gene_col}_name"`` column inserted
        immediately after ``gene_col``.
    """
    name_map = GENE_ID_TO_NAME if name_map is None else name_map
    df = df.copy()

    def _lookup(gene_id):
        if pd.isna(gene_id):
            return gene_id
        name = name_map.get(gene_id)
        return name if name and pd.notna(name) else gene_id

    name_col = f"{gene_col}_name"
    insert_pos = df.columns.get_loc(gene_col) + 1
    df.insert(insert_pos, name_col, df[gene_col].apply(_lookup))
    return df


def get_pombe_gene_interactions(genes, datasets):
    """Collect interactions involving one or more S. pombe systematic IDs.

    ``datasets`` is an iterable of ``(name, dataframe, source_col_1,
    source_col_2)`` tuples. The source columns identify the original
    interactors and are also used to determine self-interactions.
    """
    if isinstance(genes, str):
        genes = [genes]
    requested_genes = {str(gene).casefold() for gene in genes if pd.notna(gene)}

    output_columns = [
        "sp_gene_1", "sp_gene_name_1", "sp_gene_2", "sp_gene_name_2",
        "original_interactor_1",
        "original_interactor_2", "dataset", "self_interaction",
        "sp_gene_1_unknown", "sp_gene_1_PUG", "sp_gene_2_unknown",
        "sp_gene_2_PUG",
    ]
    interaction_tables = []

    def _gene_name(gene_id):
        if pd.isna(gene_id):
            return gene_id
        gene_name = GENE_ID_TO_NAME.get(gene_id)
        return gene_name if gene_name and pd.notna(gene_name) else gene_id

    for dataset_name, dataframe, source_col_1, source_col_2 in datasets:
        required_columns = {
            "sp_gene_1", "sp_gene_2", source_col_1, source_col_2,
            "sp_gene_1_unknown", "sp_gene_1_PUG",
            "sp_gene_2_unknown", "sp_gene_2_PUG",
        }
        missing_columns = required_columns - set(dataframe.columns)
        if missing_columns:
            raise KeyError(
                f"{dataset_name} is missing required columns: {sorted(missing_columns)}"
            )

        has_requested_gene = (
            dataframe["sp_gene_1"].astype("string").str.casefold().isin(requested_genes)
            | dataframe["sp_gene_2"].astype("string").str.casefold().isin(requested_genes)
        )
        subset = dataframe.loc[has_requested_gene]
        if subset.empty:
            continue

        table = pd.DataFrame({
            "sp_gene_1": subset["sp_gene_1"],
            "sp_gene_name_1": subset["sp_gene_1"].map(_gene_name),
            "sp_gene_2": subset["sp_gene_2"],
            "sp_gene_name_2": subset["sp_gene_2"].map(_gene_name),
            "original_interactor_1": subset[source_col_1],
            "original_interactor_2": subset[source_col_2],
            "sp_gene_1_unknown": subset["sp_gene_1_unknown"],
            "sp_gene_1_PUG": subset["sp_gene_1_PUG"],
            "sp_gene_2_unknown": subset["sp_gene_2_unknown"],
            "sp_gene_2_PUG": subset["sp_gene_2_PUG"],
        })
        table["dataset"] = dataset_name
        table["self_interaction"] = (
            table["original_interactor_1"].notna()
            & table["original_interactor_1"].eq(table["original_interactor_2"]).fillna(False)
        ).astype(bool)
        interaction_tables.append(table[output_columns])

    if not interaction_tables:
        return pd.DataFrame(columns=output_columns)
    return pd.concat(interaction_tables, ignore_index=True)[output_columns]


def uniprot_biomart_dataset_tailoring(biomart_path, dataset_path, dataset_name, output_dir,
                                      id_columns=('Protein1', 'Protein2')):
    """
    Build a UniProt ID -> human gene name map from an Ensembl BioMart export,
    restricted to the UniProt IDs present in a dataset of interest.

    Swiss-Prot IDs are used first; IDs still unmapped fall back to the
    TrEMBL column; IDs with no mapping get "NONE" for both fields. If a
    UniProt ID has several gene names in BioMart, the first one is kept.

    Parameters
    ----------
    biomart_path : str or Path
        Tab-separated BioMart export with "Gene name", "UniProtKB/Swiss-Prot ID"
        and "UniProtKB/TrEMBL ID" columns.
    dataset_path : str or Path
        Tab-separated dataset keyed by UniProt IDs. Leading '#' lines are skipped.
    dataset_name : str
        Short name used in the output filename.
    output_dir : str or Path
        Directory for ``uniprot_biomart_mapping_{dataset_name}.csv``.
    id_columns : sequence of str
        Dataset columns holding UniProt IDs.

    Returns
    -------
    pd.DataFrame
        Indexed by "Uniprot ID" with columns "Gene Name" and "Source".
    """
    with open(dataset_path, encoding='utf-8') as fh:
        n_comment_lines = 0
        for line in fh:
            if not line.startswith('#'):
                break
            n_comment_lines += 1

    dataset_df = pd.read_csv(dataset_path, sep='\t', skiprows=n_comment_lines,
                             usecols=list(id_columns), dtype='string')
    dataset_ids = set()
    for column in id_columns:
        dataset_ids |= set(dataset_df[column].dropna().str.strip()) - {''}

    biomart_df = pd.read_csv(
        biomart_path, sep='\t', dtype='string',
        usecols=['Gene name', 'UniProtKB/Swiss-Prot ID', 'UniProtKB/TrEMBL ID'],
    )

    def _id_to_name(id_column):
        pairs = biomart_df[[id_column, 'Gene name']].apply(lambda s: s.str.strip())
        pairs = pairs.replace('', pd.NA).dropna().drop_duplicates(subset=id_column)
        return dict(zip(pairs[id_column], pairs['Gene name']))

    swissprot_map = _id_to_name('UniProtKB/Swiss-Prot ID')
    trembl_map = _id_to_name('UniProtKB/TrEMBL ID')

    # Restrict to IDs present in the dataset.
    swissprot_subset = {i: n for i, n in swissprot_map.items() if i in dataset_ids}
    mapping_df = pd.DataFrame({
        'Gene Name': pd.Series(swissprot_subset),
        'Source': 'UniProtKB/Swiss-Prot ID',
    })

    remaining_ids = dataset_ids - set(swissprot_subset)
    trembl_subset = {i: trembl_map[i] for i in remaining_ids if i in trembl_map}
    trembl_df = pd.DataFrame({
        'Gene Name': pd.Series(trembl_subset, dtype='object'),
        'Source': 'UniProtKB/TrEMBL ID',
    })

    unmapped_ids = sorted(remaining_ids - set(trembl_subset))
    unmapped_df = pd.DataFrame({'Gene Name': 'NONE', 'Source': 'NONE'}, index=unmapped_ids)

    mapping_df = pd.concat([mapping_df, trembl_df, unmapped_df])
    mapping_df.index.name = 'Uniprot ID'

    output_path = Path(output_dir) / f'uniprot_biomart_mapping_{dataset_name}.csv'
    output_path.parent.mkdir(parents=True, exist_ok=True)
    mapping_df.to_csv(output_path)

    print(f"{dataset_name}: {len(dataset_ids)} UniProt IDs | "
          f"Swiss-Prot {len(swissprot_subset)} | TrEMBL {len(trembl_subset)} | "
          f"unmapped {len(unmapped_ids)} -> {output_path}")
    return mapping_df


def load_biogrid_physical_data(organism_name, biogrid_dir, biogrid_version='5.0.261',
                                col_a='sp_gene_1', col_b='sp_gene_2'):
    """
    Load a BIOGRID ORGANISM tab3 file for a given organism and filter it down
    to physical interactions only (``Experimental System Type`` == 'physical').

    The file is read line by line (rather than loaded whole into a dataframe)
    since these files can be large and most rows may be filtered out.

    The 'Systematic Name Interactor A'/'B' columns are renamed to ``col_a``/
    ``col_b`` and moved to the leftmost columns of the table, to match the
    project convention used for other PPI sources.

    Parameters
    ----------
    organism_name : str
        Organism name as it appears in the BIOGRID filename, e.g.
        "Schizosaccharomyces_pombe_972h".
    biogrid_dir : str or Path
        Directory containing the BIOGRID-ORGANISM tab3 file.
    biogrid_version : str, optional
        BIOGRID release version used in the filename. Default '5.0.261'.
    col_a : str, optional
        Name to use for 'Systematic Name Interactor A'. Default 'sp_gene_1'.
        Use e.g. 'sc_gene_1' for a S. cerevisiae dataset, since these
        systematic names are not S. pombe gene IDs.
    col_b : str, optional
        Name to use for 'Systematic Name Interactor B'. Default 'sp_gene_2'.

    Returns
    -------
    pd.DataFrame
        Physical interactions only, with ``col_a``/``col_b`` as the leftmost
        columns.
    """
    file_path = Path(biogrid_dir) / f"BIOGRID-ORGANISM-{organism_name}-{biogrid_version}.tab3.txt"

    with open(file_path, 'r', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        header[0] = header[0].lstrip('#')  # BIOGRID prefixes the first column with '#'

        exp_type_idx = header.index('Experimental System Type')

        physical_rows = [row for row in reader if row[exp_type_idx] == 'physical']

    df = pd.DataFrame(physical_rows, columns=header)
    df = df.rename(columns={
        'Systematic Name Interactor A': col_a,
        'Systematic Name Interactor B': col_b,
    })

    other_cols = [c for c in df.columns if c not in (col_a, col_b)]
    df = df[[col_a, col_b] + other_cols]

    return df


def extract_model_proteins_from_panther_families(BASE_DIR, gene_of_interest, panther_family_list, species_dict):
    """
    Concatenate PantherDB cluster.fasta files for several families of a gene of
    interest, keeping only sequences from the species in species_dict.

    Reads BASE_DIR / gene_dir / family / cluster.fasta, where gene_dir is
    gene_of_interest with "." replaced by "_". The PantherDB species identifier
    is the text between ">" and "|" in each record id (e.g. ">SCHPO|...").
    Writes BASE_DIR / gene_dir / concatenated_families.fasta.

    Returns the output Path and the number of sequences written.
    """
    from Bio import SeqIO

    base_dir = Path(BASE_DIR)
    gene_dir = base_dir / gene_of_interest.replace(".", "_")
    # Each value is (common_name, panther_id, ...); extra fields such as kingdom are ignored
    keep_ids = {values[1] for values in species_dict.values()}

    kept, seen = [], set()
    for family in panther_family_list:
        fasta_path = gene_dir / family / "cluster.fasta"
        if not fasta_path.exists():
            print(f"Warning: {fasta_path} not found, skipping")
            continue
        n_total = n_kept = 0
        for record in SeqIO.parse(str(fasta_path), "fasta"):
            n_total += 1
            if record.id.split("|")[0] in keep_ids and record.id not in seen:
                seen.add(record.id)
                kept.append(record)
                n_kept += 1
        print(f"{gene_of_interest} {family}: kept {n_kept} of {n_total} sequences")

    out_path = gene_dir / "concatenated_families.fasta"
    SeqIO.write(kept, str(out_path), "fasta")
    print(f"Wrote {len(kept)} sequences to {out_path}")
    return out_path, len(kept)
