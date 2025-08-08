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

# Import from local modules
from constants import IUPAC_CODES
from config import PAM_PATTERNS, DELIVERY_EFFICIENCY, ORGANS, DELIVERY_METHODS
from models.cas9 import Cas9
from models.scoring import calculate_on_target_score, off_target_score, pam_score, secondary_structure_penalty
from models.simulation import ChromatinContext, cas9_activity, repair_pathway_prob, simulate_delivery, delivery_efficiency
from utils.sequence_utils import is_valid_pam

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

# --- gRNA/PAM SUGGESTION ---
suggested_grna = None
suggested_pam = None
suggested_pos = None
suggested_score = None
suggested_pam_score = None
suggested_struct_penalty = None
if seq_record:
    pam_pattern_for_suggest = 'NGG'  # Default for suggestion
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
        st.sidebar.markdown(f"- gRNA: {suggested_grna}")
        st.sidebar.markdown(f"- PAM: {suggested_pam}")
        st.sidebar.markdown(f"- Position: {suggested_pos}")
        st.sidebar.markdown(f"- On-target score: {suggested_score:.2f} (GC proxy)")
        st.sidebar.markdown(f"- PAM score: {suggested_pam_score:.2f} (SpCas9)")
        st.sidebar.markdown(f"- Secondary structure penalty: {suggested_struct_penalty:.2f}")
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

# --- PAM Input and Cas Variant Presets ---
st.sidebar.header("7. PAM Sequence / Cas Variant")
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
                    if HAS_DFV:
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
- Finds all PAM sites (both strands) based on selected Cas variant  
- Locates gRNA binding regions (configurable mismatches)  
- Simulates NHEJ (knockout) or HDR (knock-in) editing  
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
    test_cas9 = Cas9(pam_pattern=pam_pattern if 'pam_pattern' in locals() else 'NGG')
    pam_sites = test_cas9.find_pam_sites(seq_record)
    st.write(f"**PAM sites found ({test_cas9.pam_pattern}):** {len(pam_sites)}")
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
        ATGCGTACCGGTTACCGGATCGTACCGGTTACCGGATCGTACCGGTTACCGGA

        **gRNA (20 nt):**  
        CGTACCGGTTACCGGATCGT

        **Cas9 Variant:**  
        Custom (NCC)

        This gRNA matches the region adjacent to an NCC PAM site in the sequence, ensuring a successful edit simulation.
        """)
