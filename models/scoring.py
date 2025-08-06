"""
Scoring functions for CRISPR-Cas9 target evaluation.
"""
from typing import Dict, Any
import math

def calculate_on_target_score(guide_seq: str) -> float:
    """
    Calculate an on-target score for the guide RNA.
    
    This is a simplified version inspired by the Doench 2016 algorithm.
    
    Args:
        guide_seq: 20-nucleotide guide RNA sequence (DNA)
        
    Returns:
        Score between 0 (worst) and 1 (best)
    """
    if len(guide_seq) != 20:
        return 0.0
        
    score = 0.0
    
    # Position-dependent nucleotide features (simplified)
    pos_weights = [
        (0, 'G', 0.22), (1, 'A', -0.15), (2, 'C', 0.12), (3, 'T', -0.10),
        (4, 'G', 0.10), (5, 'A', -0.12), (6, 'C', 0.08), (7, 'T', -0.08),
        (8, 'G', 0.09), (9, 'A', -0.09), (10, 'C', 0.07), (11, 'T', -0.07),
        (12, 'G', 0.06), (13, 'A', -0.06), (14, 'C', 0.05), (15, 'T', -0.05),
        (16, 'G', 0.04), (17, 'A', -0.04), (18, 'C', 0.03), (19, 'T', -0.03)
    ]
    
    for pos, nt, wt in pos_weights:
        if guide_seq[pos] == nt:
            score += wt
    
    # GC content penalty (optimal ~40-60%)
    gc = guide_seq.count('G') + guide_seq.count('C')
    gc_content = gc / 20
    if gc_content < 0.4 or gc_content > 0.8:
        score -= 0.2
    
    # T at position 20 (penalty)
    if guide_seq[19] == 'T':
        score -= 0.2
    
    # Clamp to [0,1] range
    return max(0.0, min(1.0, 0.5 + score))

def off_target_score(guide_seq: str, target_seq: str) -> float:
    """
    Calculate an off-target score between a guide and a potential off-target.
    
    Lower scores are better (0 = perfect match).
    
    Args:
        guide_seq: Guide RNA sequence
        target_seq: Potential off-target sequence
        
    Returns:
        Score between 0 (best) and 1 (worst)
    """
    if len(guide_seq) != len(target_seq):
        return 1.0
        
    score = 0.0
    for i, (g, t) in enumerate(zip(guide_seq, target_seq)):
        if g != t:
            # Mismatches near PAM (3') are more penalized
            pos_weight = 1.0 if i >= 15 else 0.5
            score += pos_weight
    
    return min(1.0, score / len(guide_seq))

def pam_score(pam_seq: str) -> float:
    """
    Score a PAM sequence based on known preferences.
    
    Args:
        pam_seq: PAM sequence to score
        
    Returns:
        Score between 0 (worst) and 1 (best)
    """
    # Example: SpCas9 NGG, N=any, G=preferred
    PAM_SCORE_MATRIX = [
        {"A": 1, "C": 1, "G": 1, "T": 1},  # N
        {"A": 0, "C": 0, "G": 1, "T": 0},  # G
        {"A": 0, "C": 0, "G": 1, "T": 0},  # G
    ]
    
    score = 0.0
    for i, base in enumerate(pam_seq):
        if i < len(PAM_SCORE_MATRIX):
            score += PAM_SCORE_MATRIX[i].get(base, 0)
        else:
            score += 1  # For longer PAMs, treat as N
            
    return score / len(pam_seq)

def secondary_structure_penalty(guide_seq: str) -> float:
    """
    Calculate a penalty for potential secondary structure in the guide RNA.
    
    Args:
        guide_seq: Guide RNA sequence
        
    Returns:
        Multiplier between 0 (worst) and 1 (best)
    """
    # Penalize long runs of G/C (proxy for hairpins)
    if "GGGG" in guide_seq or "CCCC" in guide_seq:
        return 0.7
    return 1.0
