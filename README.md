

# CRISPR-Cas9 Gene Editing Simulator

A comprehensive Streamlit application for simulating CRISPR-Cas9 gene editing experiments with advanced biological modeling.

## Features

### Core CRISPR Functionality
- **PAM Site Detection**: Supports multiple Cas variants (SpCas9, SaCas9, Cpf1) with custom PAM patterns
- **Target Site Finding**: Locates guide RNA binding sites with configurable mismatch tolerance
- **Gene Editing Simulation**: Models both NHEJ (knockout) and HDR (knock-in) editing outcomes

### Advanced Biological Modeling
- **On-target Scoring**: Implements Doench 2016-inspired algorithm for guide RNA efficiency prediction
- **Off-target Analysis**: Position-weighted mismatch scoring for off-target assessment
- **PAM Scoring**: Position-specific PAM sequence optimization
- **Secondary Structure**: Simple hairpin formation penalty assessment

### Experimental Parameters
- **Temperature Effects**: Cas9 activity modeling based on temperature (optimal at 37°C)
- **Protein Concentration**: Saturation kinetics for Cas9 protein levels
- **Cell Cycle**: HDR vs NHEJ pathway competition based on cell cycle phase

### Delivery Modeling
- **Tissue-Specific Delivery**: Organ-specific efficiency for different delivery methods
- **Delivery Methods**: LNP, AAV, and Lentivirus options with route-dependent efficiency
- **Dose-Response**: Saturating dose-response modeling

### Chromatin Context
- **Accessibility Modeling**: Chromatin state effects on editing efficiency
- **Histone Marks**: H3K27ac (activating) and H3K9me3 (repressive) mark modeling

## Project Structure

```
├── app.py                 # Main Streamlit application
├── CRISPR_v3.py          # Original single-file implementation
├── constants.py          # IUPAC codes and constants
├── config.py             # PAM patterns and delivery parameters
├── models/
│   ├── __init__.py
│   ├── cas9.py           # Core Cas9 class and methods
│   ├── scoring.py        # Scoring algorithms (Doench, off-target, etc.)
│   └── simulation.py     # Biological modeling (chromatin, delivery, etc.)
├── utils/
│   ├── __init__.py
│   ├── sequence_utils.py # DNA sequence manipulation utilities
│   └── visualization.py # Plotting and visualization functions
├── requirements.txt      # Python dependencies
└── test_functionality.py # Comprehensive test suite
```

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/Nirajan2005/CRISPR-simulator.git
   cd CRISPR-simulator
   ```

2. Create and activate a virtual environment (recommended):
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

### Running the Application

```bash
streamlit run app.py
```

### Example Simulation

**DNA Sequence:**
```
ATGCGTACCGGTTACCGGATCGTACCGGTTACCGGATCGTACCGGTTACCGGA
```

**Guide RNA (20 nt):**
```
CGTACCGGTTACCGGATCGT
```

**Cas9 Variant:** Custom (NCC)

This example demonstrates a successful CRISPR edit where the gRNA matches a region adjacent to an NCC PAM site.

### Input Options

1. **DNA Sequence Input**
   - Paste sequence directly
   - Upload FASTA file
   - Use provided examples

2. **Guide RNA Design**
   - Manual entry (15-25 nucleotides)
   - Automated suggestions from uploaded sequences

3. **Editing Parameters**
   - Edit type: Knockout (NHEJ) or Knock-in (HDR)
   - Donor template for HDR reactions
   - Mismatch tolerance (0-2 mismatches)

4. **Experimental Conditions**
   - Temperature (20-42°C)
   - Cas9 protein concentration (10-200 nM)
   - Cell cycle phase (G1, S/G2, Other)

5. **Delivery Parameters**
   - Target organ/tissue (Liver, Lung, Muscle, Brain)
   - Delivery method (LNP, AAV, Lentivirus)
   - Dose and administration route (IV, IM, IP)

6. **Chromatin Context**
   - Accessibility score (0-1)
   - Histone marks (H3K27ac, H3K9me3)

## Testing

Run the test suite to verify functionality:

```bash
python test_functionality.py
```

The test suite validates:
- Module imports and dependencies
- Core CRISPR functionality (PAM finding, target identification, editing)
- Delivery and biological modeling
- Compatibility with the original implementation

## Dependencies

- **streamlit**: Web application framework
- **biopython**: Biological sequence analysis
- **dna-features-viewer**: Graphical sequence visualization (optional)
- **matplotlib**: Plotting and visualization
- **numpy**: Numerical operations

## License

MIT License - see [LICENSE](LICENSE) for details.
