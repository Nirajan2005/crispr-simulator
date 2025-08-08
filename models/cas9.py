"""
CRISPR-Cas9 Core Model

Contains the Cas9 class and related functions for target finding and editing simulation.
"""
from typing import List, Dict, Any, Optional, Tuple
from Bio.Seq import Seq
import random

# Import at function level to avoid circular imports
# These will be imported when needed inside the methods
IUPAC_CODES = None
PAM_PATTERNS = None
is_valid_pam = None

class Cas9:
    """
    A class representing the Cas9 protein for CRISPR gene editing simulation.
    """
    
    def __init__(self, pam_pattern: str = "NGG", max_mismatches: int = 2) -> None:
        """
        Initialize the Cas9 model.
        
        Args:
            pam_pattern: The PAM pattern this Cas9 variant recognizes (IUPAC format)
            max_mismatches: Maximum number of mismatches allowed in gRNA binding
        """
        self.pam_pattern: str = pam_pattern.upper()
        self.max_mismatches: int = max_mismatches
        self.pam_sites_log: List[Dict[str, Any]] = []

    def find_pam_sites(self, seq: str) -> List[int]:
        """
        Find all PAM sites in the given sequence.
        
        Args:
            seq: DNA sequence to search for PAM sites
            
        Returns:
            List of start positions of PAM sites
        """
        # Import here to avoid circular imports
        from utils.sequence_utils import is_valid_pam
        
        seq = seq.upper()
        pam_len = len(self.pam_pattern)
        sites = []
        for i in range(len(seq) - pam_len + 1):
            pam_seq = seq[i:i+pam_len]
            if is_valid_pam(pam_seq, self.pam_pattern):
                sites.append(i)
        return sites

    def find_targets(self, seq: str, guide_rna: str) -> List[Dict[str, Any]]:
        """
        Find all potential target sites for the given guide RNA in the sequence.
        
        Args:
            seq: DNA sequence to search in
            guide_rna: Guide RNA sequence (without PAM)
            
        Returns:
            List of target site dictionaries with position, sequence, and strand info
        """
        # Import here to avoid circular imports
        from utils.sequence_utils import is_valid_pam
        
        seq = seq.upper()
        guide_rna = guide_rna.upper()
        pam = self.pam_pattern.upper()
        pam_len = len(pam)
        guide_len = len(guide_rna)
        targets: List[Dict[str, Any]] = []
        pam_sites_log: List[Dict[str, Any]] = []

        # Forward strand: gRNA immediately followed by PAM
        for i in range(len(seq) - guide_len - pam_len + 1):
            candidate = seq[i:i+guide_len]
            pam_candidate = seq[i+guide_len:i+guide_len+pam_len]
            if candidate == guide_rna:
                pam_sites_log.append({"pos": i+guide_len, "pam_seq": pam_candidate, "strand": "+", "valid": is_valid_pam(pam_candidate, pam)})
                if is_valid_pam(pam_candidate, pam):
                    targets.append({
                        "start": i,
                        "end": i+guide_len+pam_len,
                        "pam_pos": i+guide_len,
                        "mismatches": 0,
                        "strand": "+",
                        "target_seq": candidate,
                        "pam_seq": pam_candidate
                    })

        # Reverse strand: gRNA (revcomp) immediately followed by PAM
        rc_seq = str(Seq(seq).reverse_complement())
        rc_guide = str(Seq(guide_rna).reverse_complement())
        for i in range(len(rc_seq) - guide_len - pam_len + 1):
            candidate = rc_seq[i:i+guide_len]
            pam_candidate = rc_seq[i+guide_len:i+guide_len+pam_len]
            if candidate == rc_guide:
                pam_sites_log.append({"pos": i+guide_len, "pam_seq": pam_candidate, "strand": "-", "valid": is_valid_pam(pam_candidate, pam)})
                if is_valid_pam(pam_candidate, pam):
                    # Map back to original coordinates
                    orig_start = len(seq) - (i+guide_len+pam_len)
                    orig_end = len(seq) - i
                    targets.append({
                        "start": orig_start,
                        "end": orig_end,
                        "pam_pos": orig_start+guide_len,
                        "mismatches": 0,
                        "strand": "-",
                        "target_seq": candidate,
                        "pam_seq": pam_candidate
                    })
        self.pam_sites_log = pam_sites_log  # For debugging/educational output
        return sorted(targets, key=lambda x: x['start'])

    def simulate_edit(
        self,
        seq: str,
        target: Dict[str, Any],
        edit_type: str,
        donor_template: Optional[str] = None
    ) -> str:
        """
        Simulate a gene edit at the target site.
        
        Args:
            seq: Original DNA sequence
            target: Target site information
            edit_type: Type of edit ("Knockout (NHEJ)" or "Knock-in/Edit (HDR)")
            donor_template: Donor template for HDR (required for knock-in)
            
        Returns:
            Edited DNA sequence
        """
        seq = str(seq)
        if edit_type == "Knockout (NHEJ)":
            cut_site = target['start'] + len(target['target_seq']) // 2
            del_len = random.randint(1, 3)
            edited = seq[:cut_site] + seq[cut_site+del_len:]
            return edited
        elif edit_type == "Knock-in/Edit (HDR)":
            if donor_template:
                cut_site = target['start'] + len(target['target_seq']) // 2
                edited = seq[:cut_site] + donor_template + seq[cut_site:]
                return edited
            else:
                return seq
        else:
            return seq
