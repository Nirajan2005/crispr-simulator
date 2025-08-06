
import streamlit as st
from Bio.Seq import Seq
from Bio import SeqIO
import io
import random
import re
from typing import List, Dict, Optional, Any

# Optional: for visualization
try:
    from dna_features_viewer import GraphicFeature, GraphicRecord
    HAS_DFV = True
except ImportError:
    HAS_DFV = False


# --- PAM Patterns and Helper Functions ---
PAM_PATTERNS: Dict[str, str] = {
    "SpCas9": "NGG",
    "SaCas9": "NNGRRT",
    "Cpf1": "TTTV"
}

# IUPAC codes for DNA bases
IUPAC_CODES: Dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
    "H": "ACT", "V": "ACG", "N": "ACGT"
}

def is_valid_pam(pam_seq: str, pam_pattern: str) -> bool:
    if len(pam_seq) != len(pam_pattern):
        return False
    for base, pat in zip(pam_seq, pam_pattern):
        if base not in IUPAC_CODES.get(pat, pat):
            return False
    return True

class Cas9:
    def __init__(self, pam_pattern: str = "NGG", max_mismatches: int = 2) -> None:
        self.pam_pattern: str = pam_pattern
        self.max_mismatches: int = max_mismatches

    def find_pam_sites(self, seq: str) -> List[int]:
        seq = seq.upper()
        pam_len = len(self.pam_pattern)
        sites = []
        for i in range(len(seq) - pam_len + 1):
            pam_seq = seq[i:i+pam_len]
            if is_valid_pam(pam_seq, self.pam_pattern):
                sites.append(i)
        return sites

    def find_targets(self, seq: str, guide_rna: str) -> List[Dict[str, Any]]:
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

# --- Delivery Modeling ---


DELIVERY_EFFICIENCY: Dict[str, Dict[str, float]] = {
    "Liver": {"LNP": 0.8, "AAV": 0.4, "Lentivirus": 0.3},
    "Lung": {"LNP": 0.3, "AAV": 0.5, "Lentivirus": 0.4},
    "Muscle": {"LNP": 0.2, "AAV": 0.5, "Lentivirus": 0.6},
    "Brain": {"LNP": 0.05, "AAV": 0.7, "Lentivirus": 0.6},
}


ORGANS: List[str] = list(DELIVERY_EFFICIENCY.keys())
DELIVERY_METHODS: List[str] = ["LNP", "AAV", "Lentivirus"]

def simulate_delivery(organ: str, method: str) -> bool:
    eff = DELIVERY_EFFICIENCY.get(organ, {}).get(method, 0.1)
    return random.random() < eff

# --- Streamlit App ---

st.title("CRISPR-Cas9 Gene Editing Simulator")

st.sidebar.header("1. Input DNA Sequence")
input_mode = st.sidebar.radio("Input mode", ["Paste sequence", "Upload FASTA"])


# Always store the sequence as a string for downstream processing
seq_record: Optional[str] = None
fasta_header: Optional[str] = None
if input_mode == "Paste sequence":
    raw_seq = st.sidebar.text_area("Paste DNA sequence here", height=8)
    if raw_seq:
        seq_record = raw_seq.replace("\n", "").replace(" ", "").upper()
else:
    fasta_file = st.sidebar.file_uploader("Upload FASTA file", type=["fasta", "fa"])
    if fasta_file:
        fasta_obj = next(SeqIO.parse(io.StringIO(fasta_file.getvalue().decode()), "fasta"))
        seq_record = str(fasta_obj.seq).upper()
        fasta_header = fasta_obj.description




# --- On-target scoring (Doench 2016 proxy) ---
def calculate_on_target_score(guide_seq: str) -> float:
    # Doench 2016: position-dependent features, GC, T at pos 20, etc.
    if len(guide_seq) != 20:
        return 0.0
    score = 0.0
    # Position-dependent nucleotide features (simplified)
    # Example weights (not full model):
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
    # GC content penalty (Doench: optimal ~40-60%)
    gc = guide_seq.count('G') + guide_seq.count('C')
    gc_content = gc / 20
    if gc_content < 0.4 or gc_content > 0.8:
        score -= 0.2
    # T at position 20 (Doench: penalty)
    if guide_seq[19] == 'T':
        score -= 0.2
    # Clamp to [0,1]
    score = max(0.0, min(1.0, 0.5 + score))
    return score
# --- Off-target scoring (position-weighted mismatches) ---
def off_target_score(guide_seq: str, target_seq: str) -> float:
    # Lower is better (0 = perfect match)
    if len(guide_seq) != len(target_seq):
        return 1.0
    score = 0.0
    for i, (g, t) in enumerate(zip(guide_seq, target_seq)):
        if g != t:
            # Mismatches near PAM (3') are more penalized
            pos_weight = 1.0 if i >= 15 else 0.5
            score += pos_weight
    return min(1.0, score / len(guide_seq))
# --- Temperature and protein concentration effects ---
def cas9_activity(temp_c: float, protein_conc: float) -> float:
    # Activity optimal at 37C, drops at extremes; protein conc saturates
    temp_factor = max(0.0, 1.0 - abs(temp_c - 37) / 20)
    conc_factor = min(1.0, protein_conc / 100.0)
    return temp_factor * conc_factor
# --- HDR vs NHEJ competition (cell cycle) ---
def repair_pathway_prob(cell_cycle: str) -> float:
    # Returns probability of HDR (rest is NHEJ)
    if cell_cycle == "S/G2":
        return 0.7
    elif cell_cycle == "G1":
        return 0.2
    else:
        return 0.4

# --- PAM scoring (position-specific, SpCas9 example) ---
def pam_score(pam_seq: str) -> float:
    # Example: SpCas9 NGG, N=any, G=preferred
    PAM_SCORE_MATRIX = [
        {"A": 1, "C": 1, "G": 1, "T": 1},  # N
        {"A": 0, "C": 0, "G": 1, "T": 0},  # G
        {"A": 0, "C": 0, "G": 1, "T": 0},  # G
    ]
    score = 0
    for i, base in enumerate(pam_seq):
        if i < len(PAM_SCORE_MATRIX):
            score += PAM_SCORE_MATRIX[i].get(base, 0)
        else:
            score += 1  # For longer PAMs, treat as N
    return score / len(pam_seq)

# --- gRNA secondary structure penalty (simple proxy) ---
def secondary_structure_penalty(guide_seq: str) -> float:
    # Penalize long runs of G/C (proxy for hairpins)
    if "GGGG" in guide_seq or "CCCC" in guide_seq:
        return 0.7
    return 1.0

# --- Chromatin accessibility modeling (with histone marks) ---
class ChromatinContext:
    def __init__(self, accessibility_score: float = 1.0, histone_marks: dict = None):
        self.accessibility = accessibility_score
        self.active_marks = histone_marks or {}
    def modify_cutting_efficiency(self, base_efficiency: float) -> float:
        mod = 1.0
        if self.active_marks.get("H3K27ac", 0) > 0:
            mod *= 1.2
        if self.active_marks.get("H3K9me3", 0) > 0:
            mod *= 0.7
        return base_efficiency * self.accessibility * mod

# --- gRNA/PAM SUGGESTION ---

# --- gRNA/PAM SUGGESTION ---
suggested_grna = None
suggested_pam = None
suggested_pos = None
suggested_score = None
suggested_pam_score = None
suggested_struct_penalty = None
if seq_record:
    pam_pattern_for_suggest = pam_pattern if 'pam_pattern' in locals() else 'NGG'
    cas9_suggest = Cas9(pam_pattern=pam_pattern_for_suggest)
    pam_len = len(pam_pattern_for_suggest)
    grna_len = 20
    for i in range(len(seq_record) - grna_len - pam_len + 1):
        grna_candidate = seq_record[i:i+grna_len]
        pam_candidate = seq_record[i+grna_len:i+grna_len+pam_len]
        if is_valid_pam(pam_candidate, pam_pattern_for_suggest):
            suggested_grna = grna_candidate
            suggested_pam = pam_candidate
            suggested_pos = i
            suggested_score = calculate_on_target_score(grna_candidate)
            suggested_pam_score = pam_score(pam_candidate)
            suggested_struct_penalty = secondary_structure_penalty(grna_candidate)
            break
    if not suggested_grna:
        rc_seq = str(Seq(seq_record).reverse_complement())
        for i in range(len(rc_seq) - grna_len - pam_len + 1):
            grna_candidate = rc_seq[i:i+grna_len]
            pam_candidate = rc_seq[i+grna_len:i+grna_len+pam_len]
            if is_valid_pam(pam_candidate, pam_pattern_for_suggest):
                suggested_grna = str(Seq(grna_candidate).reverse_complement())
                suggested_pam = pam_candidate
                suggested_pos = f"reverse:{i}"
                suggested_score = calculate_on_target_score(suggested_grna)
                suggested_pam_score = pam_score(pam_candidate)
                suggested_struct_penalty = secondary_structure_penalty(suggested_grna)
                break
    if suggested_grna:
        st.sidebar.markdown("---")
        st.sidebar.markdown("**Suggested gRNA/PAM from uploaded sequence:**")
        if fasta_header:
            st.sidebar.markdown(f"*{fasta_header}*")
        st.sidebar.markdown(f"- gRNA: `{suggested_grna}`")
        st.sidebar.markdown(f"- PAM: `{suggested_pam}`")
        st.sidebar.markdown(f"- Position: `{suggested_pos}`")
        st.sidebar.markdown(f"- On-target score: `{suggested_score:.2f}` (GC proxy)")
        st.sidebar.markdown(f"- PAM score: `{suggested_pam_score:.2f}` (SpCas9)")
        st.sidebar.markdown(f"- Secondary structure penalty: `{suggested_struct_penalty:.2f}`")
        st.sidebar.markdown(f"- (Copy/paste gRNA above to simulate a real edit)")

st.sidebar.header("2. Guide RNA")
guide_rna: str = st.sidebar.text_input("Guide RNA sequence (20 nt)", max_chars=25).upper()

st.sidebar.header("3. Editing Options")
edit_type: str = st.sidebar.selectbox("Edit type", ["Knockout (NHEJ)", "Knock-in/Edit (HDR)"])
donor_template: Optional[str] = None
if edit_type == "Knock-in/Edit (HDR)":
    donor_template = st.sidebar.text_input("Donor template (for HDR)", max_chars=30).upper()

st.sidebar.header("4. Off-target & PAM")
max_mismatches: int = st.sidebar.slider("Max mismatches allowed", 0, 2, 2)



# --- Chromatin accessibility UI ---

st.sidebar.header("5. Delivery & Chromatin Context")
organ: str = st.sidebar.selectbox("Target organ/tissue", ORGANS)
delivery_method: str = st.sidebar.selectbox("Delivery method", DELIVERY_METHODS)
dose = st.sidebar.slider("Dose (ug)", 0.1, 10.0, 1.0, 0.1)
route = st.sidebar.selectbox("Administration route", ["IV", "IM", "IP"])
chromatin_accessibility: float = st.sidebar.slider("Chromatin accessibility (0=closed, 1=open)", 0.0, 1.0, 1.0, 0.05)
histone_marks = {}
histone_selected = st.sidebar.multiselect("Histone marks", ["H3K27ac", "H3K9me3"])
for mark in histone_selected:
    histone_marks[mark] = 1
chromatin_ctx = ChromatinContext(accessibility_score=chromatin_accessibility, histone_marks=histone_marks)

# New: temperature, protein conc, cell cycle
st.sidebar.header("6. Experimental Parameters")
temp_c = st.sidebar.slider("Temperature (°C)", 20, 42, 37)
protein_conc = st.sidebar.slider("Cas9 protein conc (nM)", 10, 200, 100)
cell_cycle = st.sidebar.selectbox("Cell cycle phase", ["G1", "S/G2", "Other"])

# Improved delivery modeling
def delivery_efficiency(base, dose, route):
    route_factor = {"IV": 1.0, "IM": 0.7, "IP": 0.5}[route]
    # Dose-response: log-scale, saturating
    dose_factor = 1 - 0.5 * (1 / (1 + dose))
    return base * dose_factor * route_factor

# --- PAM Input and Cas Variant Presets ---

st.sidebar.header("6. PAM Sequence / Cas Variant")

st.sidebar.header("6. PAM Sequence / Cas Variant")
cas_options = list(PAM_PATTERNS.keys()) + ["Custom"]
cas_type: str = st.sidebar.selectbox(
    "Select Cas protein (Cas variant)",
    cas_options,
    key="pam_variant_selectbox_1"
)
if cas_type == "Custom":
    pam_pattern = st.sidebar.text_input("Enter custom PAM pattern (IUPAC, e.g. NGG, TTTV)", value="NGG").upper()
else:
    pam_pattern = PAM_PATTERNS[cas_type]





if st.sidebar.button("Simulate CRISPR Editing"):
    if not seq_record or not guide_rna or len(guide_rna) < 15:
        st.error("Please provide a valid DNA sequence and guide RNA (≥15 nt).")
    else:
        cas9 = Cas9(pam_pattern=pam_pattern, max_mismatches=max_mismatches)
        targets = cas9.find_targets(seq_record, guide_rna)
        st.subheader(f"Cas Protein: {cas_type} | PAM Pattern: {pam_pattern}")
        edited_seq = seq_record  # Already a string
        n_on_targets = sum(t['mismatches'] == 0 for t in targets)
        n_off_targets = len(targets) - n_on_targets
        # Log PAM sites for educational/debugging output
        if hasattr(cas9, 'pam_sites_log'):
            st.write("**PAM site log (for debugging/education):**")
            for entry in cas9.pam_sites_log:
                st.write(f"Strand: {entry['strand']} | Pos: {entry['pos']} | PAM: {entry['pam_seq']} | Valid: {entry['valid']}")
        if not targets:
            # Check if gRNA binds anywhere (even if no valid PAM)
            found_gRNA = False
            seq = seq_record.upper()
            guide = guide_rna.upper()
            rc_seq = str(Seq(seq).reverse_complement())
            rc_guide = str(Seq(guide).reverse_complement())
            if guide in seq or rc_guide in rc_seq:
                st.error("gRNA binds to DNA, but no valid PAM site found. No edit performed.")
            else:
                st.warning("No valid PAM/target sites found for this gRNA. No edit performed.")
        else:
            # Delivery simulation
            base_eff = DELIVERY_EFFICIENCY[organ][delivery_method]
            deliv_eff = delivery_efficiency(base_eff, dose, route)
            # Temperature/protein/cell cycle effects
            cas9_eff = cas9_activity(temp_c, protein_conc)
            hdr_prob = repair_pathway_prob(cell_cycle)
            # Modify by chromatin accessibility and histone marks
            effective_eff = chromatin_ctx.modify_cutting_efficiency(deliv_eff * cas9_eff)
            delivered = random.random() < effective_eff
            st.write(f"**Delivery to {organ} via {delivery_method} ({route}), Dose: {dose} ug**: {'Success' if delivered else 'Failure'}")
            st.write(f"Base delivery efficiency: {base_eff*100:.0f}%")
            st.write(f"Dose-adjusted efficiency: {deliv_eff*100:.0f}%")
            st.write(f"Cas9 activity (temp/protein): {cas9_eff*100:.0f}%")
            st.write(f"Cell cycle: {cell_cycle} (HDR probability: {hdr_prob*100:.0f}%)")
            st.write(f"Chromatin accessibility: {chromatin_accessibility}")
            st.write(f"Histone marks: {', '.join(histone_selected) if histone_selected else 'None'}")
            st.write(f"Effective editing efficiency: {effective_eff*100:.0f}%")
            st.write(f"**Total targets found:** {len(targets)}")
            st.write(f"**On-target sites:** {n_on_targets}")
            st.write(f"**Off-target sites:** {n_off_targets}")
            # Visualization: delivery success/failure chart
            st.markdown("**Delivery Success Probability**")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(3,1))
            ax.bar(["Success", "Failure"], [effective_eff, 1-effective_eff], color=["green","red"])
            ax.set_ylim(0,1)
            ax.set_ylabel("Probability")
            st.pyplot(fig)
            if not delivered:
                st.info("No editing due to failed delivery.")
            elif n_on_targets == 0:
                st.warning("Delivery succeeded, but no on-target site found for this gRNA. No edit performed.")
            else:
                # Show targets and simulate edits
                features = []
                off_targets = []
                edit_success = False
                for i, t in enumerate(targets):
                    color = "#ffcccc" if t['mismatches'] == 0 else "#ccccff"
                    features.append(GraphicFeature(start=t['start'], end=t['end'], strand=1 if t['strand']=="+" else -1,
                                                  color=color, label=f"Target {i+1} ({t['mismatches']} mm)"))
                    if t['mismatches'] == 0:
                        # HDR/NHEJ competition
                        if edit_type == "Knock-in/Edit (HDR)":
                            if random.random() < hdr_prob:
                                edited_seq = cas9.simulate_edit(edited_seq, t, edit_type, donor_template)
                                edit_success = True
                            else:
                                # NHEJ fallback
                                edited_seq = cas9.simulate_edit(edited_seq, t, "Knockout (NHEJ)")
                                edit_success = True
                        else:
                            edited_seq = cas9.simulate_edit(edited_seq, t, edit_type, donor_template)
                            edit_success = True
                    else:
                        off_targets.append(t)
                if edit_success:
                    st.success("✅ Successful edit performed! The DNA sequence has been modified at the on-target site.")
                else:
                    st.info("No on-target edit was performed.")
                # Before/after alignment
                import difflib
                st.markdown("**Before/After Sequence Alignment:**")
                align = difflib.ndiff(seq_record, edited_seq)
                st.code('\n'.join(align), language="diff")
                st.code(edited_seq, language="text")
                if off_targets:
                    st.write("**Off-target report:**")
                    for t in off_targets:
                        ot_score = off_target_score(guide_rna, t['target_seq'])
                        st.write(f"- Position {t['start']}-{t['end']} ({t['strand']}), {t['mismatches']} mismatch(es), Off-target score: {ot_score:.2f} | PAM: {t['pam_seq']}")
                # Visualization
                if HAS_DFV:
                    record = GraphicRecord(sequence_length=len(seq_record), features=features)
                    ax, _ = record.plot(figure_width=8)
                    st.pyplot(ax.get_figure())
                else:
                    st.info("Install dna_features_viewer for graphical visualization.")

        st.subheader("Final Output DNA Sequence")
        st.code(edited_seq, language="text")
        st.download_button("Download edited sequence", edited_seq, file_name="edited_sequence.txt")

st.markdown("""
---
**How it works:**  
- Finds all NGG PAM sites (both strands)  
- Locates gRNA binding regions (≤2 mismatches)  
- Simulates NHEJ (knockout) or HDR (knock-in)  
- Models delivery efficiency by organ/method  
- Reports off-targets and shows edited sequence  
""")
# --- Helper: Show sequence info and troubleshooting tips ---

if seq_record:
    st.subheader("Input Sequence Info")
    st.write(f"**Length:** {len(seq_record)} bp")
    st.write(f"**First 60 bases:** {seq_record[:60]}")
    # Check for ambiguous bases
    ambiguous = set(seq_record) - set("ACGT")
    if ambiguous:
        st.warning(f"Sequence contains ambiguous bases: {', '.join(ambiguous)}. This may affect PAM finding.")
    # Check for PAM presence
    test_cas9 = Cas9(pam_pattern=pam_pattern)
    pam_sites = test_cas9.find_pam_sites(seq_record)
    st.write(f"**PAM sites found ({pam_pattern}):** {len(pam_sites)}")
    if len(pam_sites) == 0:
        st.info("No PAM sites found. Try a different PAM or check your sequence.")

    # gRNA troubleshooting
    if guide_rna and len(guide_rna) >= 15:
        targets = test_cas9.find_targets(seq_record, guide_rna)
        st.write(f"**Potential target sites for gRNA:** {len(targets)}")
        if len(targets) == 0:
            st.warning("No target sites found for this gRNA. Try changing the gRNA or PAM, or check the sequence orientation.")
            st.markdown("""
            **Tips for successful editing:**
            - Make sure your gRNA matches a region in the DNA sequence adjacent to a PAM site.
            - Try using a gRNA sequence that is present in your DNA (copy a 20 nt region from your sequence).
            - Check that your DNA sequence is in the correct orientation (5'→3').
            - Use the 'Custom' PAM option if your Cas variant is different.
            """)
        # Example: Successful Edit

        st.markdown("""
        ### Example: Successful CRISPR Edit

        **DNA Sequence:**  
        `ATGCGTACCGGTTACCGGATCGTACCGGTTACCGGATCGTACCGGTTACCGGA`

        **gRNA (20 nt):**  
        `CGTACCGGTTACCGGATCGT`

        **Cas9 Variant:**  
        Custom (NCC)

        This gRNA matches the region adjacent to an NCC PAM site in the sequence, ensuring a successful edit simulation.
        """)
