"""
Sequence manipulation utilities for CRISPR-Cas9 simulation.
"""
from typing import Dict, Optional, List, Tuple
from Bio.Seq import Seq

# Import IUPAC_CODES directly from constants to avoid circular imports
IUPAC_CODES = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'U': 'U',
    'R': 'AG',
    'Y': 'CT',
    'S': 'GC',
    'W': 'AT',
    'K': 'GT',
    'M': 'AC',
    'B': 'CGT',
    'D': 'AGT',
    'H': 'ACT',
    'V': 'ACG',
    'N': 'ACGT',
    '.': '-',
    '-': '-',
    ' ': ' '
}

def is_valid_pam(pam_seq: str, pam_pattern: str) -> bool:
    """
    Check if a PAM sequence matches the given PAM pattern.
    
    Args:
        pam_seq: PAM sequence to check
        pam_pattern: PAM pattern in IUPAC format
        
    Returns:
        True if the sequence matches the pattern, False otherwise
    """
    if len(pam_seq) != len(pam_pattern):
        return False
        
    for base, pat in zip(pam_seq.upper(), pam_pattern.upper()):
        if base not in IUPAC_CODES.get(pat, pat):
            return False
    return True

def reverse_complement(seq: str) -> str:
    """
    Get the reverse complement of a DNA sequence.
    
    Args:
        seq: DNA sequence
        
    Returns:
        Reverse complement of the sequence
    """
    return str(Seq(seq).reverse_complement())

def find_pam_sites(sequence: str, pam_pattern: str) -> list:
    """
    Find all PAM sites in a DNA sequence.
    
    Args:
        sequence: DNA sequence to search in
        pam_pattern: PAM pattern in IUPAC format
        
    Returns:
        List of (position, pam_sequence) tuples
    """
    sites = []
    pam_len = len(pam_pattern)
    sequence = sequence.upper()
    
    for i in range(len(sequence) - pam_len + 1):
        pam_candidate = sequence[i:i+pam_len]
        if is_valid_pam(pam_candidate, pam_pattern):
            sites.append((i, pam_candidate))
    
    return sites
