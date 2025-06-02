# 🧬 BRCA1 Sequence Analyzer

An interactive bioinformatics tool for analyzing BRCA1 gene sequences, detecting variants, and providing clinical significance annotations through ClinVar integration.

## 📋 Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Algorithms](#algorithms)
- [Visualization](#visualization)
- [API Reference](#api-reference)
- [Contributing](#contributing)
- [License](#license)

## ✨ Features

### Core Functionality
- **Sequence Alignment**: Implements both Smith-Waterman (local) and Needleman-Wunsch (global) algorithms
- **Variant Detection**: Identifies SNVs, insertions, and deletions
- **Clinical Annotation**: Integrates with ClinVar database for variant significance
- **Interactive Visualizations**: Rich plotting with Plotly for alignment and variant analysis
- **Web Interface**: User-friendly Streamlit application

### Analysis Capabilities
- 🔍 **Sequence Quality Assessment**: Comprehensive metrics and quality indicators
- 📊 **Variant Distribution Analysis**: Visual breakdown of variant types and clinical significance
- 🎨 **Interactive Alignment Viewer**: Color-coded nucleotide visualization
- 📈 **Statistical Reporting**: Detailed analysis reports in JSON and CSV formats
- ⚠️ **Clinical Risk Assessment**: Pathogenic variant identification and alerts

## 🚀 Installation

### Prerequisites
- Python 3.8 or higher
- pip package manager

### Quick Start
You can access the web-based tool at 
```bash
link coming soon
```

or for setting up with localhost with these steps

1. **Clone the repository:**
```bash
git clone <repository-url>
cd brca1-analyzer
```

2. **. Setup the environments**
#### Option 1: Using Setup Scripts
these scripts are the compact of the manual one

**For Linux/macOS:**
```bash
chmod +x run_app.sh
./run_app.sh
```

**For Windows:**
```cmd
run_app.bat
```

#### Option 2: Manual Installation

2. **Create virtual environment:**
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. **Install dependencies:**
```bash
pip install --upgrade pip
pip install -r requirements.txt
```

4. **Run the application:**
```bash
streamlit run main.py
```

### Dependencies
all the dependencies are stated in the  `requirements.txt` file.

## 📖 Usage

### Web Interface

1. **Launch the application:**
   ```bash
   streamlit run main.py
   ```

2. **Access the interface:**
   Open your browser to `http://localhost:8501`

3. **Input your sequence:**
   - Upload a FASTA file, or
   - Paste sequence directly into the text area

4. **Configure analysis:**
   - Choose alignment method (Global, Local, or Both)
   - Select visualization options
   - Set analysis parameters

5. **Analyze and explore results:**
   - View alignment visualizations
   - Examine variant distribution charts
   - Review clinical significance annotations
   - Download comprehensive reports

### Command Line Usage

```python
from brca1_analyzer import BRCA1Analyzer

# Initialize analyzer
analyzer = BRCA1Analyzer()

# Analyze sequence
sample_sequence = "ATGGATTTATCTGCTCTTCGCGTT..."
results = analyzer.analyze_sequence(sample_sequence, alignment_type="both")

# Access results
print(f"Variants found: {len(results['variants'])}")
for variant in results['variants']:
    print(f"Position {variant['position']}: {variant['variant_name']}")
```

## 📁 Project Structure

```
brca1_analyzer/
├── brca1_analyzer/
│   ├── __init__.py           # Package initialization
│   ├── aligner.py            # Sequence alignment algorithms
│   ├── analyzer.py           # Main analyzer class
│   ├── clinvar.py           # ClinVar integration
│   └── visuals.py           # Visualization functions
├── sample/
│   └── patient_sample.fasta  # Sample input file
├── main.py                   # Streamlit application
├── requirements.txt          # Python dependencies
├── run_app.sh               # Linux/macOS launcher
├── run_app.bat              # Windows launcher
└── README.md                # This file
```

## 🧮 Algorithms

### Smith-Waterman Algorithm (Local Alignment)
- **Purpose**: Finds optimal local alignments between sequences
- **Use Case**: Identifying regions of high similarity
- **Scoring**: Configurable match/mismatch/gap penalties
- **Output**: Best local alignment with position information

### Needleman-Wunsch Algorithm (Global Alignment)
- **Purpose**: Optimal global alignment of entire sequences
- **Use Case**: Complete sequence comparison
- **Scoring**: Dynamic programming approach
- **Output**: End-to-end alignment with comprehensive statistics

### Variant Detection Pipeline
1. **Sequence Alignment**: Align sample against BRCA1 reference
2. **Difference Identification**: Scan aligned sequences for mismatches
3. **Variant Classification**: Categorize as SNV, insertion, or deletion
4. **Clinical Annotation**: Query ClinVar for significance data

## 📊 Visualization

### Alignment Visualization
- **Color-coded nucleotides**: A (red), T (teal), C (blue), G (green)
- **Mismatch highlighting**: Red X markers for variant positions
- **Interactive tooltips**: Position and nucleotide information
- **Scalable display**: Adjustable length for large sequences

### Variant Distribution Charts
- **Pie charts**: Variant type and clinical significance breakdown
- **Scatter plots**: Variant positions and significance scores
- **Multi-panel dashboard**: Comprehensive variant overview

### Quality Metrics
- **Alignment scores**: Comparative algorithm performance
- **Identity percentages**: Sequence similarity measures
- **Gap distribution**: Alignment quality indicators
- **Quality gauge**: Overall sequence assessment score

## 🔧 API Reference

### Core Classes

#### `BRCA1Analyzer`
Main analysis orchestrator.

```python
analyzer = BRCA1Analyzer()
results = analyzer.analyze_sequence(sequence, alignment_type="both")
```

#### `SequenceAligner`
Implements alignment algorithms.

```python
aligner = SequenceAligner(match_score=2, mismatch_score=-1, gap_penalty=-1)
result, seq1, seq2 = aligner.smith_waterman(ref_seq, sample_seq)
```

#### `ClinVarIntegrator`
Handles variant annotation.

```python
clinvar = ClinVarIntegrator()
variant_info = clinvar.search_variant("BRCA1", "c.185A>C")
```

### Visualization Functions

```python
from brca1_analyzer.visuals import (
    create_alignment_visualization,
    create_variant_distribution_chart,
    create_sequence_quality_metrics
)

# Generate interactive plots
alignment_fig = create_alignment_visualization(alignment_data)
variant_fig = create_variant_distribution_chart(variants)
quality_fig = create_sequence_quality_metrics(results)
```

## 🧪 Sample Data

The project includes sample data in `sample/patient_sample.fasta`:
- BRCA1 sequence variants for testing
- Demonstrates typical analysis workflow
- Shows pathogenic variant detection

## ⚙️ Configuration

### Alignment Parameters
```python
# Customize scoring scheme
aligner = SequenceAligner(
    match_score=2,      # Reward for matches
    mismatch_score=-1,  # Penalty for mismatches
    gap_penalty=-1      # Penalty for gaps
)
```

### Visualization Options
```python
# Control display parameters
alignment_viz = create_alignment_visualization(
    alignment_data,
    max_display_length=100  # Limit displayed sequence length
)
```

## 📋 Output Formats

### JSON Report
```json
{
  "analysis_date": "2024-01-15 10:30:00",
  "sample_info": {
    "length": 300,
    "reference_length": 300
  },
  "variants": [...],
  "summary": {
    "total_variants": 5,
    "pathogenic_variants": 2
  }
}
```

### CSV Export
Tabular format with columns:
- Position, Reference, Sample, Type, Variant Name, Clinical Significance, ClinVar ID

## 🚨 Clinical Interpretation

### Significance Categories
- **Pathogenic**: Disease-causing variants
- **Likely Pathogenic**: Probably disease-causing
- **Uncertain Significance**: Unknown clinical impact
- **Likely Benign**: Probably not disease-causing
- **Benign**: Not disease-causing

### Risk Assessment
The tool provides automated alerts for:
- Pathogenic variants detected
- Multiple variants of uncertain significance
- Novel variants not in ClinVar database

## 🔬 Technical Notes

### Performance Considerations
- Sequences are truncated to 100 nucleotides for visualization
- Large sequences may require increased processing time
- ClinVar queries are cached to improve performance

### Limitations
- Reference sequence is simplified for demonstration
- ClinVar integration uses mock data in demo mode
- Real deployment requires NCBI API configuration

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-analysis`)
3. Commit changes (`git commit -am 'Add new analysis feature'`)
4. Push to branch (`git push origin feature/new-analysis`)
5. Create Pull Request

### Development Setup
```bash
# Install development dependencies
pip install -e .
pip install pytest black flake8

# Run tests
pytest tests/

# Format code
black brca1_analyzer/
flake8 brca1_analyzer/
```

## 🆘 Support

For questions, issues, or contributions:
- Open an issue on GitHub
- Check existing documentation
- Review sample implementations

## 🙏 Acknowledgments

- NCBI ClinVar database for variant annotations
- BioPython community for sequence analysis tools
- Streamlit team for the web framework
- Plotly for interactive visualizations

---

**Note**: This tool is for research and educational purposes. Clinical decisions should always involve qualified healthcare professionals and validated diagnostic procedures.