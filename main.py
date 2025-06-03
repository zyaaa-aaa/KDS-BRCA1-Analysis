import time
import json
from io import StringIO
from Bio import SeqIO
import pandas as pd
import streamlit as st
from brca1_analyzer import BRCA1Analyzer
from brca1_analyzer.visuals import (
    create_alignment_visualization, 
    create_variant_distribution_chart,
    create_sequence_quality_metrics,
    create_nucleotide_composition_chart,
    create_variant_summary_table,
    display_alignment_details
)

def safe_get_alignment_value(align_data, key, default=0.0):
    """
    Safely extract numeric values from alignment data
    """
    try:
        value = align_data.get(key, default)
        if isinstance(value, (int, float)):
            return float(value)
        elif hasattr(value, key):  # AlignmentResult object
            return float(getattr(value, key, default))
        else:
            return float(default)
    except (ValueError, TypeError, AttributeError):
        return float(default)

def create_streamlit_app():
    """
    Create Streamlit web interface with enhanced visualizations
    """
    st.set_page_config(
        page_title="BRCA1 Sequence Analyzer",
        page_icon="🧬",
        layout="wide"
    )
    
    st.title("🧬 BRCA1 Sequence Analyzer")
    st.markdown("*Analisis Komparatif Varian Sekuens Gen BRCA1 dengan Integrasi ClinVar*")
    
    # Add information banner
    st.info("🔬 **Real-time integration with NCBI GenBank and ClinVar databases for accurate variant analysis**")
    
    # Initialize analyzer with error handling
    if 'analyzer' not in st.session_state:
        try:
            with st.spinner("Initializing BRCA1 Analyzer..."):
                # Initialize with email (required by NCBI)
                st.session_state.analyzer = BRCA1Analyzer(
                    email="abc@gmail.com",  # Replace with actual email
                    api_key=None  # Add NCBI API key if available
                )
                
                # Check if initialization was successful
                if not hasattr(st.session_state.analyzer, 'reference_sequence') or not st.session_state.analyzer.reference_sequence:
                    st.error("❌ Failed to load BRCA1 reference sequence from NCBI")
                    st.stop()
                else:
                    st.success("✅ BRCA1 Analyzer initialized successfully")
                    
        except Exception as e:
            st.error(f"Failed to initialize analyzer: {e}")
            st.exception(e)
            st.stop()
    
    # Sidebar for parameters
    with st.sidebar:
        st.header("⚙️ Analysis Parameters")
        
        st.subheader("🧬 Gene Selection")
        gene_name = st.selectbox(
            "Choose gene to analyze",
            ["BRCA1"],
            index=0,
            help="Select the gene for sequence analysis"
        )
        
        # Gene switching functionality
        if gene_name != "BRCA1":
            if st.button("🔄 Switch Gene", help=f"Load {gene_name} reference sequence from NCBI"):
                with st.spinner(f"Loading {gene_name} reference sequence from NCBI..."):
                    try:
                        success = st.session_state.analyzer.set_custom_gene(gene_name)
                        if success:
                            st.success(f"✅ Successfully switched to {gene_name}")
                            st.rerun()
                        else:
                            st.error(f"❌ Failed to load {gene_name} reference sequence")
                    except Exception as e:
                        st.error(f"❌ Error switching to {gene_name}: {e}")
        
        st.subheader("🔧 Alignment Parameters")
        alignment_type = st.selectbox(
            "Choose alignment method",
            ["Global (Needleman-Wunsch)", "Local (Smith-Waterman)", "Both"],
            index=2,
            help="Global: End-to-end alignment. Local: Best local similarity."
        )
        
        # Advanced alignment parameters
        with st.expander("⚙️ Advanced Alignment Settings"):
            match_score = st.slider("Match Score", 1, 5, 2, help="Score for matching nucleotides")
            mismatch_score = st.slider("Mismatch Penalty", -5, -1, -1, help="Penalty for mismatched nucleotides")
            gap_penalty = st.slider("Gap Penalty", -5, -1, -1, help="Penalty for insertions/deletions")
            
            if st.button("Update Alignment Parameters"):
                try:
                    # Update aligner parameters
                    st.session_state.analyzer.aligner.match_score = match_score
                    st.session_state.analyzer.aligner.mismatch_score = mismatch_score
                    st.session_state.analyzer.aligner.gap_penalty = gap_penalty
                    st.success("✅ Alignment parameters updated")
                except Exception as e:
                    st.error(f"❌ Error updating parameters: {e}")
        
        st.subheader("🏥 Variant Annotation")
        annotation_source = st.selectbox(
            "Database source",
            ["ClinVar"],
            disabled=True,
            help="Currently using NCBI ClinVar for variant classification"
        )
        
        st.subheader("📊 Visualization Options")
        show_alignment_viz = st.checkbox("Show Alignment Visualization", value=True)
        show_variant_charts = st.checkbox("Show Variant Distribution", value=True)
        show_quality_metrics = st.checkbox("Show Quality Metrics", value=True)
        show_composition = st.checkbox("Show Nucleotide Composition", value=True)
        max_display_length = st.slider("Max Display Length", 50, 500, 100, help="Maximum sequence length for visualization")
        
        # Reference sequence information
        if hasattr(st.session_state.analyzer, 'reference_sequence') and st.session_state.analyzer.reference_sequence:
            try:
                ref_info = st.session_state.analyzer.get_reference_info()
                if ref_info['length'] > 0:
                    st.subheader("📋 Reference Sequence Info")
                    st.write(f"**Gene:** {gene_name}")
                    st.write(f"**Length:** {ref_info['length']:,} bp")
                    st.write(f"**GC Content:** {ref_info['gc_content']:.1f}%")
                    
                    # Nucleotide composition
                    composition = ref_info.get('nucleotide_composition', {})
                    if composition:
                        st.write("**Composition:**")
                        for nt in ['A', 'T', 'G', 'C']:
                            st.write(f"  {nt}: {composition.get(nt, 0):.1f}%")
            except Exception as e:
                st.warning(f"Could not load reference info: {e}")
    
    # Main content area
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("📁 Upload FASTA File")
        uploaded_file = st.file_uploader(
            "Choose a FASTA file",
            type=['fasta', 'fa', 'txt', 'seq'],
            help="Upload your DNA sequence in FASTA format"
        )
        
        st.subheader("✏️ Paste Sequence")
        sequence_input = st.text_area(
            "Or paste your DNA sequence here:",
            height=150,
            placeholder="ATGGATTTATCTGCTCTTCGCGTT...",
            help="Enter DNA sequence using A, T, C, G nucleotides (case insensitive)"
        )
        
        # Example sequences
        st.subheader("🧪 Example Sequences")
        col_ex1, col_ex2 = st.columns(2)
        
        with col_ex1:
            if st.button("📋 Load Wild Type Example", help="Load a sample BRCA1 sequence segment"):
                # BRCA1 exon 11 partial sequence (wild type)
                example_seq = """ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAG
                TGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTT
                GCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAA
                AAGGAGCCTACAAGAAAGTACGAGATTTAGTCAACTTGTTGAAGAGCTATTGAAAATCATTTGTGCTTTT
                CAGCTTGACACAGGTTTGGAGTATGCAAACAGCTATAATTTTGCAAAAAAGGAAAATAACTGCCCTAAAC
                ATGAAGTGAAGATGTTCCGTTAAACAGATTTAACATTTCTCAGTACCCCTCCATTGTGTAGATTGCATCT
                TCTTCTCACCACCTAAATAAAGGGCCAGAAGAAGTACGAGATTTAGTCAACTTGTTGAAGAGCTATTGAA"""
                
                # Clean the sequence
                example_seq = ''.join(c.upper() for c in example_seq if c.upper() in 'ATCG')
                st.session_state.example_sequence = example_seq
                st.success(f"✅ Wild type example loaded! ({len(example_seq)} nucleotides)")
        
        with col_ex2:
            if st.button("🧪 Load Variant Example", help="Load a sequence with known BRCA1 variants"):
                # BRCA1 sequence with known pathogenic variant (5382insC)
                example_seq = """ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAG
                TGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTT
                GCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAA
                AAGGAGCCTACAAGAAAGTACGAGATTTAGTCAACTTGTTGAAGAGCTATTGAAAATCATTTGTGCTTTT
                CAGCTTGACACAGGTTTGGAGTATGCAAACAGCTATAATTTTGCAAAAAAGGAAAATAACTGCCCTAAAC
                ATGAAGTGAAGATGTTCCGTTAAACAGATTTAACATTTCTCAGTACCCCTCCATTGTGTAGATTGCATCT
                TCTTCTCACCACCTAAATAAAGGGCCAGAAGAAGTACGAGATTTAGTCAACTTGTTGAAGAGCTATTGAA"""
                
                # Clean and introduce a known variant
                example_seq = ''.join(c.upper() for c in example_seq if c.upper() in 'ATCG')
                # Introduce c.185A>G variant
                if len(example_seq) > 184:
                    example_seq = example_seq[:184] + "G" + example_seq[185:]
                st.session_state.example_sequence = example_seq
                st.success(f"✅ Variant example loaded! ({len(example_seq)} nucleotides)")
        
        if st.button("🗑️ Clear Example", help="Clear loaded example sequence"):
            if 'example_sequence' in st.session_state:
                del st.session_state.example_sequence
                st.info("Example sequence cleared")
    
    with col2:
        st.subheader("🔍 Analysis Status")
        
        # Process input with priority: example > file > manual input
        sample_sequence = ""
        sequence_source = ""
        sequence_info = {}
        
        # Check for example sequence first
        if 'example_sequence' in st.session_state:
            sample_sequence = st.session_state.example_sequence
            sequence_source = "Example"
            st.info(f"📝 Using example sequence: {len(sample_sequence):,} nucleotides")
            
        elif uploaded_file:
            try:
                content = StringIO(uploaded_file.getvalue().decode("utf-8"))
                records = list(SeqIO.parse(content, "fasta"))
                if records:
                    sample_sequence = str(records[0].seq).upper()
                    sequence_source = uploaded_file.name
                    sequence_info = {
                        'id': records[0].id,
                        'description': records[0].description,
                        'record_count': len(records)
                    }
                    st.success(f"✅ Loaded from file: {len(sample_sequence):,} nucleotides")
                    
                    # Show sequence metadata
                    with st.expander("📄 Sequence Information"):
                        st.write(f"**Sequence ID:** {sequence_info['id']}")
                        if sequence_info['description']:
                            st.write(f"**Description:** {sequence_info['description']}")
                        st.write(f"**Records in file:** {sequence_info['record_count']}")
                        if sequence_info['record_count'] > 1:
                            st.info("Multiple sequences found. Using the first sequence.")
                else:
                    st.error("❌ No valid sequences found in file")
            except Exception as e:
                st.error(f"❌ Error reading file: {e}")
                st.write("**Supported formats:** FASTA (.fasta, .fa, .txt)")
                
        elif sequence_input.strip():
            # Clean sequence input
            original_length = len(sequence_input.strip())
            sample_sequence = ''.join(c.upper() for c in sequence_input if c.upper() in 'ATCGN')
            sequence_source = "Manual input"
            
            if sample_sequence:
                st.success(f"✅ Input sequence: {len(sample_sequence):,} nucleotides")
                
                # Validation feedback
                invalid_chars = set(sequence_input.upper()) - set('ATCGN \n\t\r')
                if invalid_chars:
                    st.warning(f"⚠️ Removed {len(invalid_chars)} invalid character types: {', '.join(sorted(invalid_chars))}")
                
                removed_count = original_length - len(sample_sequence)
                if removed_count > 0:
                    st.info(f"ℹ️ Removed {removed_count} non-nucleotide characters")
            else:
                st.error("❌ No valid nucleotides found in input")
        
        # Sequence quality check
        if sample_sequence:
            quality_issues = []
            
            # Check sequence length
            if len(sample_sequence) < 20:
                quality_issues.append("Sequence very short (< 20 nucleotides)")
            elif len(sample_sequence) < 50:
                quality_issues.append("Sequence short (< 50 nucleotides)")
            
            # Check for excessive N's
            n_percentage = (sample_sequence.count('N') / len(sample_sequence)) * 100
            if n_percentage > 10:
                quality_issues.append(f"High ambiguous nucleotides ({n_percentage:.1f}% N's)")
            
            # Display quality warnings
            if quality_issues:
                st.warning("⚠️ **Sequence Quality Issues:**")
                for issue in quality_issues:
                    st.write(f"• {issue}")
        
        # Analysis button
        analyze_button = st.button(
            "🚀 Analyze Sequence", 
            type="primary",
            disabled=not sample_sequence or len(sample_sequence) < 10,
            help="Analyze sequence against BRCA1 reference" if sample_sequence else "Please provide a sequence first"
        )
        
        if analyze_button and sample_sequence and len(sample_sequence) >= 10:
            # Clear previous results
            if 'analysis_results' in st.session_state:
                del st.session_state.analysis_results
            
            with st.spinner("🔬 Analyzing sequence..."):
                try:
                    # Create progress tracking
                    progress_container = st.container()
                    with progress_container:
                        progress_bar = st.progress(0)
                        status_text = st.empty()
                        
                        # Step 1: Preprocessing
                        status_text.text("🔧 Preprocessing sequence...")
                        progress_bar.progress(10)
                        time.sleep(0.3)
                        
                        # Step 2: Sequence alignment
                        status_text.text("🧬 Performing sequence alignment...")
                        progress_bar.progress(30)
                        
                        # Perform analysis
                        alignment_map = {
                            "Global (Needleman-Wunsch)": "global",
                            "Local (Smith-Waterman)": "local",
                            "Both": "both"
                        }
                        
                        results = st.session_state.analyzer.analyze_sequence(
                            sample_sequence,
                            alignment_map[alignment_type]
                        )
                        
                        progress_bar.progress(60)
                        status_text.text("🔍 Identifying variants...")
                        time.sleep(0.3)
                        
                        progress_bar.progress(80)
                        status_text.text("🏥 Annotating with ClinVar...")
                        time.sleep(0.5)
                        
                        # Add metadata to results
                        results['sequence_source'] = sequence_source
                        results['gene_analyzed'] = gene_name
                        results['sequence_info'] = sequence_info
                        results['analysis_parameters'] = {
                            'alignment_type': alignment_type,
                            'match_score': match_score,
                            'mismatch_score': mismatch_score,
                            'gap_penalty': gap_penalty
                        }
                        
                        progress_bar.progress(100)
                        status_text.text("✅ Analysis completed!")
                        
                        # Store results
                        st.session_state.analysis_results = results
                        
                        # Clear progress indicators after a moment
                        time.sleep(0.5)
                        progress_bar.empty()
                        status_text.empty()
                        
                    st.success("✅ **Analysis completed successfully!**")
                    
                    # Quick summary
                    variants_count = len(results.get('variants', []))
                    if variants_count > 0:
                        st.info(f"🔍 Found {variants_count} variant(s) for further analysis")
                    else:
                        st.info("✅ No variants detected - sequence matches reference")
                        
                except Exception as e:
                    st.error(f"❌ **Analysis failed:** {str(e)}")
                    with st.expander("🐛 Error Details"):
                        st.exception(e)
                    
        elif analyze_button:
            if not sample_sequence:
                st.error("❌ Please provide a DNA sequence")
            elif len(sample_sequence) < 10:
                st.error("❌ Sequence too short (minimum 10 nucleotides required)")
    
    # Display results section
    if 'analysis_results' in st.session_state:
        results = st.session_state.analysis_results
        
        st.divider()
        st.header("📊 Analysis Results")
        
        # Summary metrics in columns
        metric_col1, metric_col2, metric_col3, metric_col4 = st.columns(4)
        
        with metric_col1:
            st.metric(
                "Sample Length",
                f"{results.get('sample_length', 0):,} bp",
                help="Length of analyzed sequence"
            )
        
        with metric_col2:
            st.metric(
                "Reference Length", 
                f"{results.get('reference_length', 0):,} bp",
                help=f"Length of {results.get('gene_analyzed', 'reference')} reference sequence"
            )
        
        with metric_col3:
            variants_found = len(results.get('variants', []))
            st.metric(
                "Variants Detected",
                variants_found,
                help="Total number of variants found"
            )
        
        with metric_col4:
            # Count pathogenic variants
            pathogenic_count = 0
            for v in results.get('variants', []):
                clinical_info = v.get('clinical_info')
                if clinical_info and hasattr(clinical_info, 'clinical_significance'):
                    if 'pathogenic' in clinical_info.clinical_significance.lower():
                        pathogenic_count += 1
            
            st.metric(
                "Pathogenic Variants",
                pathogenic_count,
                delta=pathogenic_count if pathogenic_count > 0 else None,
                delta_color="inverse" if pathogenic_count > 0 else "normal",
                help="Number of clinically significant pathogenic variants"
            )
        
        # Sequence metadata
        st.subheader("🧬 Sequence Information")
        info_col1, info_col2 = st.columns(2)
        
        with info_col1:
            st.write(f"**Source:** {results.get('sequence_source', 'Unknown')}")
            st.write(f"**Gene:** {results.get('gene_analyzed', 'BRCA1')}")
            
            # Analysis parameters
            params = results.get('analysis_parameters', {})
            if params:
                st.write(f"**Alignment:** {params.get('alignment_type', 'Unknown')}")
                
        with info_col2:
            ref_info = results.get('reference_info', {})
            sample_info = results.get('sample_info', {})
            st.write(f"**Reference GC:** {ref_info.get('gc_content', 0):.1f}%")
            st.write(f"**Sample GC:** {sample_info.get('gc_content', 0):.1f}%")
            
            # Show GC content difference
            if ref_info.get('gc_content') and sample_info.get('gc_content'):
                gc_diff = abs(ref_info['gc_content'] - sample_info['gc_content'])
                if gc_diff > 5:
                    st.write(f"⚠️ **GC difference:** {gc_diff:.1f}%")
        
        # Visualization sections
        tab1, tab2, tab3, tab4 = st.tabs(["📈 Quality Metrics", "🔗 Alignments", "🧬 Variants", "📊 Composition"])
        
        with tab1:
            if show_quality_metrics:
                try:
                    quality_fig = create_sequence_quality_metrics(results)
                    st.plotly_chart(quality_fig, use_container_width=True, key="quality_metrics_chart")
                except Exception as e:
                    st.error(f"Error creating quality metrics: {e}")
        
        with tab2:
            # Alignment results with enhanced visualization
            if 'alignment_results' in results and results['alignment_results']:
                alignment_tabs = st.tabs([
                    name.title() for name in results['alignment_results'].keys()
                ])
                
                for tab, (align_type, align_data) in zip(alignment_tabs, results['alignment_results'].items()):
                    with tab:
                        # Safe access to alignment data using helper function
                        score = safe_get_alignment_value(align_data, 'score', 0.0)
                        identity = safe_get_alignment_value(align_data, 'identity', 0.0)
                        gaps = safe_get_alignment_value(align_data, 'gaps', 0)
                        alignment_length = safe_get_alignment_value(align_data, 'alignment_length', 0)
                        
                        # Alignment metrics
                        align_col1, align_col2, align_col3, align_col4 = st.columns(4)
                        with align_col1:
                            st.metric("Alignment Score", f"{score:.1f}")
                        with align_col2:
                            st.metric("Identity", f"{identity:.1f}%")
                        with align_col3:
                            st.metric("Gaps", int(gaps))
                        with align_col4:
                            st.metric("Length", int(alignment_length))
                        
                        # Alignment visualization
                        if show_alignment_viz:
                            st.subheader(f"🎨 {align_type.title()} Alignment Visualization")
                            try:
                                align_fig = create_alignment_visualization(align_data, max_display_length)
                                st.plotly_chart(align_fig, use_container_width=True, key=f"alignment_{align_type}_viz")
                            except Exception as e:
                                st.error(f"Error creating alignment visualization: {e}")
                        
                        # Alignment details in expandable section
                        with st.expander("🔍 View Detailed Alignment"):
                            try:
                                display_alignment_details(align_data, max_length=max_display_length)
                            except Exception as e:
                                st.error(f"Error displaying alignment details: {e}")
            else:
                st.info("No alignment results available")
        
        with tab3:
            # Variant analysis with enhanced visualizations
            variants = results.get('variants', [])
            if variants:
                # Variant visualizations
                if show_variant_charts:
                    try:
                        variant_fig = create_variant_distribution_chart(variants)
                        st.plotly_chart(variant_fig, use_container_width=True, key="variant_distribution_chart")
                    except Exception as e:
                        st.error(f"Error creating variant chart: {e}")
                
                # Variant table
                st.subheader("📋 Variant Details")
                try:
                    variant_df = create_variant_summary_table(variants)
                    st.dataframe(variant_df, use_container_width=True)
                except Exception as e:
                    st.error(f"Error creating variant table: {e}")
                    
                    # Fallback simple table
                    st.subheader("📋 Variant Summary (Simplified)")
                    variant_data = []
                    for i, variant in enumerate(variants):
                        clinical_info = variant.get('clinical_info')
                        variant_data.append({
                            'Position': variant.get('position', i+1),
                            'Reference': variant.get('reference', 'N/A'),
                            'Sample': variant.get('sample', 'N/A'),
                            'Type': variant.get('type', 'Unknown'),
                            'Variant Name': variant.get('variant_name', 'Unknown'),
                            'Clinical Significance': getattr(clinical_info, 'clinical_significance', 'Unknown') if clinical_info else 'Unknown',
                            'ClinVar ID': getattr(clinical_info, 'clinvar_id', 'N/A') if clinical_info else 'N/A'
                        })
                    
                    df = pd.DataFrame(variant_data)
                    st.dataframe(df, use_container_width=True)
                
                # Clinical significance analysis
                st.subheader("🔍 Clinical Significance Analysis")
                
                pathogenic_variants = []
                uncertain_variants = []
                benign_variants = []
                
                for v in variants:
                    clinical_info = v.get('clinical_info')
                    if clinical_info and hasattr(clinical_info, 'clinical_significance'):
                        sig = clinical_info.clinical_significance.lower()
                        if 'pathogenic' in sig:
                            pathogenic_variants.append(v)
                        elif 'uncertain' in sig:
                            uncertain_variants.append(v)
                        elif 'benign' in sig:
                            benign_variants.append(v)
                
                # Display findings based on significance
                if pathogenic_variants:
                    st.error("⚠️ **Pathogenic Variants Detected**")
                    st.write("This sample contains pathogenic variants that may be clinically significant:")
                    
                    for variant in pathogenic_variants:
                        clinical_info = variant['clinical_info']
                        with st.container():
                            st.markdown(f"**• {variant.get('variant_name', 'Unknown')}** at position {variant.get('position', 'Unknown')}")
                            st.write(f"   **Significance:** {clinical_info.clinical_significance}")
                            if hasattr(clinical_info, 'description'):
                                st.write(f"   **Description:** {clinical_info.description}")
                            if hasattr(clinical_info, 'clinvar_id'):
                                st.write(f"   **ClinVar ID:** {clinical_info.clinvar_id}")
                            if hasattr(clinical_info, 'review_status'):
                                st.write(f"   **Review Status:** {clinical_info.review_status}")
                            st.write("")
                            
                elif uncertain_variants:
                    st.warning("❓ **Variants of Uncertain Significance Detected**")
                    st.write("Found variants with uncertain clinical significance:")
                    for variant in uncertain_variants[:3]:  # Show first 3
                        st.write(f"• **{variant.get('variant_name', 'Unknown')}**")
                    if len(uncertain_variants) > 3:
                        st.write(f"... and {len(uncertain_variants) - 3} more")
                        
                elif benign_variants:
                    st.success("✅ **Only Benign Variants Detected**")
                    st.write("All detected variants are classified as benign or likely benign.")
                    
                else:
                    st.success("✅ **No Pathogenic Variants Detected**")
                    st.write("No known pathogenic variants were found in this sample.")
                
            else:
                st.info("ℹ️ **No variants detected in this sequence.**")
                st.write("**Possible explanations:**")
                st.write("• The sequence perfectly matches the reference")
                st.write("• The sequence quality is very high")
                st.write("• The analyzed region contains no known variants")
                st.write("• The alignment parameters may need adjustment")
        
        with tab4:
            if show_composition:
                try:
                    composition_fig = create_nucleotide_composition_chart(results)
                    st.plotly_chart(composition_fig, use_container_width=True, key="composition_chart")
                except Exception as e:
                    st.error(f"Error creating composition chart: {e}")
        
        # Download reports section
        st.divider()
        st.subheader("💾 Download Analysis Reports")
        
        try:
            # Generate comprehensive report
            report_data = {
                'analysis_metadata': {
                    'analysis_date': time.strftime('%Y-%m-%d %H:%M:%S'),
                    'gene_analyzed': results.get('gene_analyzed', 'BRCA1'),
                    'sequence_source': results.get('sequence_source', 'Unknown'),
                    'analyzer_version': '1.0.0',
                    'databases_used': ['NCBI GenBank', 'ClinVar']
                },
                'sequence_information': {
                    'sample_length': results.get('sample_length', 0),
                    'reference_length': results.get('reference_length', 0),
                    'sample_gc_content': results.get('sample_info', {}).get('gc_content', 0),
                    'reference_gc_content': results.get('reference_info', {}).get('gc_content', 0),
                    'nucleotide_composition': results.get('sample_info', {}).get('nucleotide_composition', {})
                },
                'alignment_results': {
                    k: {
                        'score': safe_get_alignment_value(v, 'score', 0),
                        'identity': safe_get_alignment_value(v, 'identity', 0),
                        'gaps': safe_get_alignment_value(v, 'gaps', 0),
                        'alignment_length': safe_get_alignment_value(v, 'alignment_length', 0)
                    } for k, v in results.get('alignment_results', {}).items()
                },
                'variant_analysis': {
                    'total_variants': len(variants),
                    'pathogenic_variants': len(pathogenic_variants),
                    'uncertain_variants': len(uncertain_variants),
                    'benign_variants': len(benign_variants),
                    'variant_types': {
                        vtype: len([v for v in variants if v.get('type') == vtype]) 
                        for vtype in set(v.get('type', 'Unknown') for v in variants)
                    }
                },
                'clinical_summary': {
                    'risk_assessment': 'High' if pathogenic_variants else 'Low' if not variants else 'Uncertain',
                    'pathogenic_variant_names': [v.get('variant_name', 'Unknown') for v in pathogenic_variants],
                    'recommendations': [
                        'Consult with a genetic counselor' if pathogenic_variants else 'No immediate action required',
                        'Consider family screening' if pathogenic_variants else 'Continue routine screening'
                    ]
                }
            }
            
            # Download buttons
            download_col1, download_col2, download_col3 = st.columns(3)
            
            with download_col1:
                st.download_button(
                    label="📋 Download JSON Report",
                    data=json.dumps(report_data, indent=2),
                    file_name=f"{results.get('gene_analyzed', 'brca1').lower()}_analysis_{int(time.time())}.json",
                    mime="application/json",
                    help="Complete analysis results in JSON format"
                )
            
            with download_col2:
                # Generate CSV report for variants
                if variants:
                    try:
                        variant_df = create_variant_summary_table(variants)
                        if not variant_df.empty:
                            csv_data = variant_df.to_csv(index=False)
                            st.download_button(
                                label="📊 Download CSV Report",
                                data=csv_data,
                                file_name=f"{results.get('gene_analyzed', 'brca1').lower()}_variants_{int(time.time())}.csv",
                                mime="text/csv",
                                help="Variant details in CSV format"
                            )
                        else:
                            st.button("📊 CSV Report", disabled=True, help="No variants to export")
                    except Exception:
                        st.button("📊 CSV Report", disabled=True, help="Error generating CSV")
                else:
                    st.button("📊 CSV Report", disabled=True, help="No variants detected")
            
            with download_col3:
                # Generate text summary report
                summary_lines = [
                    f"BRCA1 Sequence Analysis Summary",
                    f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}",
                    f"",
                    f"Gene: {results.get('gene_analyzed', 'BRCA1')}",
                    f"Sample Length: {results.get('sample_length', 0):,} bp",
                    f"Variants Found: {len(variants)}",
                    f"Pathogenic Variants: {len(pathogenic_variants)}",
                    f"",
                    f"Clinical Assessment:",
                    f"Risk Level: {report_data['clinical_summary']['risk_assessment']}",
                    f"",
                ]
                
                if pathogenic_variants:
                    summary_lines.append("Pathogenic Variants Detected:")
                    for variant in pathogenic_variants:
                        summary_lines.append(f"- {variant.get('variant_name', 'Unknown')}")
                else:
                    summary_lines.append("No pathogenic variants detected.")
                
                summary_lines.extend([
                    f"",
                    f"Note: This analysis is for research purposes only.",
                    f"Consult with healthcare professionals for clinical decisions."
                ])
                
                summary_text = '\n'.join(summary_lines)
                
                st.download_button(
                    label="📄 Download Summary",
                    data=summary_text,
                    file_name=f"{results.get('gene_analyzed', 'brca1').lower()}_summary_{int(time.time())}.txt",
                    mime="text/plain",
                    help="Human-readable analysis summary"
                )
                
        except Exception as e:
            st.error(f"Error generating reports: {e}")
            st.write("Please try the analysis again or contact support.")
    
    # Help and information section
    st.divider()
    
    # Footer with information
    with st.expander("ℹ️ About This Tool"):
        st.markdown("""
        ### BRCA1 Sequence Analyzer
        
        This tool performs comprehensive analysis of DNA sequences against the BRCA1 reference sequence using:
        
        **🧬 Sequence Alignment:**
        - Needleman-Wunsch (Global): Best for full-length sequences
        - Smith-Waterman (Local): Best for partial sequences or finding regions of similarity
        
        **🏥 Variant Classification:**
        - Real-time integration with NCBI ClinVar database
        - Clinical significance assessment
        - Pathogenicity prediction
        
        **📊 Analysis Features:**
        - Nucleotide composition analysis
        - Quality metrics assessment
        - Interactive visualizations
        - Comprehensive reporting
        
        **⚠️ Important Disclaimers:**
        - This tool is for research and educational purposes only
        - Results should not be used for clinical diagnosis
        - Consult healthcare professionals for medical decisions
        - Variant classifications may change as new evidence emerges
        
        **🔗 Data Sources:**
        - Reference sequences: NCBI GenBank
        - Variant annotations: NCBI ClinVar
        - Gene information: NCBI Gene database
        """)
    
    # Technical information
    with st.expander("🔧 Technical Details"):
        st.markdown("""
        ### Algorithm Parameters
        
        **Alignment Scoring:**
        - Match Score: +2 (configurable)
        - Mismatch Penalty: -1 (configurable)
        - Gap Penalty: -1 (configurable)
        
        **Variant Detection:**
        - Minimum sequence length: 10 nucleotides
        - Supported variant types: SNV, insertion, deletion
        - HGVS nomenclature compliance
        
        **Quality Metrics:**
        - Sequence identity percentage
        - Alignment coverage
        - Gap distribution analysis
        - GC content comparison
        
        **Performance Notes:**
        - Analysis time depends on sequence length
        - ClinVar queries may take several seconds
        - Maximum recommended sequence length: 10,000 bp
        """)
    
    # Contact and support information
    st.markdown("""
    <div style='text-align: center; color: gray; margin-top: 2rem;'>
        <p>🧬 <strong>BRCA1 Sequence Analyzer v1.0</strong></p>
        <p>Powered by NCBI GenBank and ClinVar databases</p>
        <p><small>For research and educational purposes only • Not for clinical diagnosis</small></p>
    </div>
    """, unsafe_allow_html=True)

if __name__ == "__main__":
    create_streamlit_app()