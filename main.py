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
    create_sequence_quality_metrics
)

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
    st.markdown("*Analisis Komparatif Varian Sekuens Gen BRCA1*")
    
    # Initialize analyzer
    if 'analyzer' not in st.session_state:
        st.session_state.analyzer = BRCA1Analyzer()
    
    # Sidebar for parameters
    with st.sidebar:
        st.header("Analysis Parameters")
        
        st.subheader("Reference Sequence Accession")
        ref_accession = st.text_input("NCBI Accession", value="NM_007294.4", disabled=True)
        
        st.subheader("Alignment Type")
        alignment_type = st.selectbox(
            "Choose alignment method",
            ["Global (Needleman-Wunsch)", "Local (Smith-Waterman)", "Both"],
            index=2
        )
        
        st.subheader("Variant Annotation Source")
        annotation_source = st.selectbox(
            "Database source",
            ["ClinVar"],
            disabled=True
        )
        
        st.subheader("Visualization Options")
        show_alignment_viz = st.checkbox("Show Alignment Visualization", value=True)
        show_variant_charts = st.checkbox("Show Variant Distribution", value=True)
        show_quality_metrics = st.checkbox("Show Quality Metrics", value=True)
    
    # Main content area
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("📁 Upload FASTA File")
        uploaded_file = st.file_uploader(
            "Choose a FASTA file",
            type=['fasta', 'fa', 'txt'],
            help="Upload your DNA sequence in FASTA format"
        )
        
        st.subheader("✏️ Paste Sequence")
        sequence_input = st.text_area(
            "Or paste your DNA sequence here:",
            height=150,
            placeholder="ATGGATTTATCTGCTCTTCGCGTT..."
        )
    
    with col2:
        st.subheader("🔍 Analysis Results")
        
        # Process input
        sample_sequence = ""
        if uploaded_file:
            try:
                content = StringIO(uploaded_file.getvalue().decode("utf-8"))
                records = list(SeqIO.parse(content, "fasta"))
                if records:
                    sample_sequence = str(records[0].seq)
                    st.success(f"✅ Loaded sequence: {len(sample_sequence)} nucleotides")
            except Exception as e:
                st.error(f"❌ Error reading file: {e}")
        
        elif sequence_input.strip():
            # Clean sequence input
            sample_sequence = ''.join(c.upper() for c in sequence_input if c.upper() in 'ATCG')
            if sample_sequence:
                st.success(f"✅ Input sequence: {len(sample_sequence)} nucleotides")
        
        # Analysis button
        if st.button("🚀 Analyze Sequence", type="primary"):
            if sample_sequence:
                with st.spinner("Analyzing sequence..."):
                    # Simulate processing time
                    progress_bar = st.progress(0)
                    for i in range(100):
                        time.sleep(0.01)
                        progress_bar.progress(i + 1)
                    
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
                    
                    st.session_state.analysis_results = results
                    st.success("✅ Analysis completed!")
            else:
                st.error("❌ Please provide a DNA sequence")
    
    # Display results
    if 'analysis_results' in st.session_state:
        results = st.session_state.analysis_results
        
        st.divider()
        st.header("📊 Analysis Results")
        
        # Summary metrics
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric(
                "Sample Length",
                f"{results['sample_length']:,} bp"
            )
        
        with col2:
            st.metric(
                "Reference Length", 
                f"{results['reference_length']:,} bp"
            )
        
        with col3:
            variants_found = len(results['variants'])
            st.metric(
                "Variants Detected",
                variants_found
            )
        
        with col4:
            pathogenic_count = sum(
                1 for v in results['variants'] 
                if v.get('clinical_info') and 'pathogenic' in v['clinical_info'].clinical_significance.lower()
            )
            st.metric(
                "Pathogenic Variants",
                pathogenic_count,
                delta=pathogenic_count - (variants_found - pathogenic_count) if variants_found > 0 else 0
            )
        
        # Visualizations
        if show_quality_metrics:
            st.subheader("📈 Quality Metrics")
            quality_fig = create_sequence_quality_metrics(results)
            st.plotly_chart(quality_fig, use_container_width=True, key="quality_metrics_chart")
        
        # Alignment results with visualization
        if 'alignment_results' in results:
            st.subheader("🔗 Alignment Results")
            
            alignment_tabs = st.tabs([
                name.title() for name in results['alignment_results'].keys()
            ])
            
            for tab, (align_type, align_data) in zip(alignment_tabs, results['alignment_results'].items()):
                with tab:
                    align_result = align_data['result']
                    
                    # Alignment metrics
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Alignment Score", f"{align_result.score:.1f}")
                    with col2:
                        st.metric("Identity", f"{align_result.identity:.1f}%")
                    with col3:
                        st.metric("Gaps", align_result.gaps)
                    
                    # Alignment visualization
                    if show_alignment_viz:
                        st.subheader(f"🎨 {align_type.title()} Alignment Visualization")
                        align_fig = create_alignment_visualization(align_data)
                        st.plotly_chart(align_fig, use_container_width=True, key=f"alignment_{align_type}_viz")
                    
                    # Alignment preview (text)
                    if len(align_data['aligned_reference']) < 200:
                        st.text("Alignment Preview:")
                        st.code(f"Ref: {align_data['aligned_reference'][:100]}...")
                        st.code(f"Sam: {align_data['aligned_sample'][:100]}...")
        
        # Variant analysis with visualizations
        if results['variants']:
            st.subheader("🧬 Detected Variants")
            
            # Variant visualizations
            if show_variant_charts:
                variant_fig = create_variant_distribution_chart(results['variants'])
                st.plotly_chart(variant_fig, use_container_width=True, key="variant_distribution_chart")
            
            # Variant table
            variant_data = []
            for variant in results['variants']:
                clinical_info = variant.get('clinical_info')
                variant_data.append({
                    'Position': variant['position'],
                    'Reference': variant['reference'],
                    'Sample': variant['sample'],
                    'Type': variant['type'],
                    'Variant Name': variant['variant_name'],
                    'Clinical Significance': clinical_info.clinical_significance if clinical_info else 'Unknown',
                    'ClinVar ID': clinical_info.clinvar_id if clinical_info else 'N/A'
                })
            
            df = pd.DataFrame(variant_data)
            st.dataframe(df, use_container_width=True)
            
            # Key findings
            st.subheader("🔍 Key Findings")
            
            pathogenic_variants = [
                v for v in results['variants']
                if v.get('clinical_info') and 'pathogenic' in v['clinical_info'].clinical_significance.lower()
            ]
            
            if pathogenic_variants:
                st.error("⚠️ **Pathogenic Variants Detected**")
                st.write("This sample contains pathogenic variants that may be clinically significant:")
                for variant in pathogenic_variants:
                    clinical_info = variant['clinical_info']
                    st.write(f"• **{clinical_info.variant_name}** - {clinical_info.clinical_significance}")
                    st.write(f"  *{clinical_info.description}*")
            else:
                st.success("✅ **No Pathogenic Variants Detected**")
                st.write("No known pathogenic variants were found in this sample.")
            
            # Download results
            st.subheader("💾 Download Report")
            
            # Generate report
            report_data = {
                'analysis_date': time.strftime('%Y-%m-%d %H:%M:%S'),
                'sample_info': {
                    'length': results['sample_length'],
                    'reference_length': results['reference_length']
                },
                'variants': variant_data,
                'summary': {
                    'total_variants': len(results['variants']),
                    'pathogenic_variants': len(pathogenic_variants)
                }
            }
            
            col1, col2 = st.columns(2)
            with col1:
                st.download_button(
                    label="📋 Download JSON Report",
                    data=json.dumps(report_data, indent=2),
                    file_name=f"brca1_analysis_{int(time.time())}.json",
                    mime="application/json"
                )
            
            with col2:
                # Generate CSV report
                csv_data = pd.DataFrame(variant_data).to_csv(index=False)
                st.download_button(
                    label="📊 Download CSV Report",
                    data=csv_data,
                    file_name=f"brca1_variants_{int(time.time())}.csv",
                    mime="text/csv"
                )
        
        else:
            st.info("ℹ️ No variants detected in this sequence.")
            
            # Still show quality metrics even without variants
            if show_quality_metrics and not st.session_state.get('quality_shown', False):
                st.subheader("📈 Sequence Quality")
                quality_fig = create_sequence_quality_metrics(results)
                st.plotly_chart(quality_fig, use_container_width=True, key="quality_metrics_no_variants")
                st.session_state.quality_shown = True

if __name__ == "__main__":
    create_streamlit_app()