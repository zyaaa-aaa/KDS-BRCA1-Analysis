import plotly.graph_objects as go
import pandas as pd
from typing import Dict, List, Optional
from plotly.subplots import make_subplots
import streamlit as st

def create_alignment_visualization(alignment_data: Dict, max_display_length: int = 100) -> go.Figure:
    """
    Create interactive alignment visualization
    """
    ref_seq = alignment_data.get('aligned_reference', '')
    sample_seq = alignment_data.get('aligned_sample', '')
    
    if not ref_seq or not sample_seq:
        fig = go.Figure()
        fig.add_annotation(
            text="No alignment data available",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=16)
        )
        fig.update_layout(
            title="Sequence Alignment Visualization",
            height=300
        )
        return fig
    
    # Truncate for display if too long
    if len(ref_seq) > max_display_length:
        ref_seq = ref_seq[:max_display_length] + "..."
        sample_seq = sample_seq[:max_display_length] + "..."
    
    # Create alignment visualization
    fig = go.Figure()
    
    # Color mapping for nucleotides
    color_map = {'A': '#FF6B6B', 'T': '#4ECDC4', 'C': '#45B7D1', 'G': '#96CEB4', '-': '#DDD'}
    
    # Reference sequence
    ref_colors = [color_map.get(nt, '#DDD') for nt in ref_seq]
    sample_colors = [color_map.get(nt, '#DDD') for nt in sample_seq]
    
    # Create heatmap-like visualization
    positions = list(range(len(ref_seq)))
    
    # Reference track
    fig.add_trace(go.Scatter(
        x=positions,
        y=[1]*len(positions),
        mode='markers',
        marker=dict(
            color=ref_colors,
            size=15,
            symbol='square',
            line=dict(width=1, color='black')
        ),
        text=[f"Pos {i+1}: {nt}" for i, nt in enumerate(ref_seq)],
        hovertemplate="Reference<br>%{text}<extra></extra>",
        name="Reference",
        showlegend=True
    ))
    
    # Sample track
    fig.add_trace(go.Scatter(
        x=positions,
        y=[0]*len(positions),
        mode='markers',
        marker=dict(
            color=sample_colors,
            size=15,
            symbol='square',
            line=dict(width=1, color='black')
        ),
        text=[f"Pos {i+1}: {nt}" for i, nt in enumerate(sample_seq)],
        hovertemplate="Sample<br>%{text}<extra></extra>",
        name="Sample",
        showlegend=True
    ))
    
    # Highlight mismatches
    mismatches = []
    for i, (r, s) in enumerate(zip(ref_seq, sample_seq)):
        if r != s:
            mismatches.append(i)
    
    if mismatches:
        fig.add_trace(go.Scatter(
            x=mismatches,
            y=[0.5]*len(mismatches),
            mode='markers',
            marker=dict(
                color='red',
                size=20,
                symbol='x',
                line=dict(width=2, color='darkred')
            ),
            name="Variants",
            hovertemplate="Variant at position %{x}<extra></extra>",
            showlegend=True
        ))
    
    fig.update_layout(
        title="Sequence Alignment Visualization",
        xaxis_title="Position",
        yaxis=dict(
            tickvals=[0, 1],
            ticktext=["Sample", "Reference"],
            range=[-0.3, 1.3]
        ),
        height=300,
        showlegend=True,
        hovermode='closest'
    )
    
    return fig

def create_variant_distribution_chart(variants: List[Dict]) -> go.Figure:
    """
    Create variant distribution visualization
    """
    if not variants:
        fig = go.Figure()
        fig.add_annotation(
            text="No variants detected",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=16)
        )
        fig.update_layout(
            title="Variant Distribution",
            height=300
        )
        return fig
    
    # Group variants by type
    variant_types = {}
    for variant in variants:
        v_type = variant.get('type', 'Unknown')
        if v_type not in variant_types:
            variant_types[v_type] = []
        variant_types[v_type].append(variant)
    
    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=("Variant Types", "Clinical Significance", 
                       "Variant Positions", "Significance by Position"),
        specs=[[{"type": "pie"}, {"type": "pie"}],
               [{"type": "scatter"}, {"type": "scatter"}]]
    )
    
    # Variant types pie chart
    type_counts = {k: len(v) for k, v in variant_types.items()}
    fig.add_trace(
        go.Pie(
            labels=list(type_counts.keys()),
            values=list(type_counts.values()),
            name="Types"
        ),
        row=1, col=1
    )
    
    # Clinical significance pie chart
    significance_counts = {}
    for variant in variants:
        clinical_info = variant.get('clinical_info')
        if clinical_info and hasattr(clinical_info, 'clinical_significance'):
            sig = clinical_info.clinical_significance
            significance_counts[sig] = significance_counts.get(sig, 0) + 1
        else:
            significance_counts['Unknown'] = significance_counts.get('Unknown', 0) + 1
    
    if significance_counts:
        colors = {
            'Pathogenic': '#FF4444',
            'Likely pathogenic': '#FF8800',
            'Uncertain significance': '#FFDD00',
            'Likely benign': '#88DD00',
            'Benign': '#44DD44',
            'Unknown': '#CCCCCC'
        }
        
        sig_colors = [colors.get(sig, '#CCCCCC') for sig in significance_counts.keys()]
        
        fig.add_trace(
            go.Pie(
                labels=list(significance_counts.keys()),
                values=list(significance_counts.values()),
                marker_colors=sig_colors,
                name="Significance"
            ),
            row=1, col=2
        )
    
    # Variant positions scatter plot
    positions = [v.get('position', 0) for v in variants]
    types = [v.get('type', 'Unknown') for v in variants]
    
    # Create color mapping for types
    unique_types = list(set(types))
    type_colors = {t: i for i, t in enumerate(unique_types)}
    colors_list = [type_colors[t] for t in types]
    
    fig.add_trace(
        go.Scatter(
            x=positions,
            y=[1]*len(positions),
            mode='markers',
            marker=dict(
                size=10,
                color=colors_list,
                colorscale='viridis',
                showscale=False
            ),
            text=[f"{v.get('type', 'Unknown')}: {v.get('variant_name', 'Unknown')}" for v in variants],
            hovertemplate="%{text}<extra></extra>",
            name="Positions"
        ),
        row=2, col=1
    )
    
    # Clinical significance by position
    if significance_counts and len(significance_counts) > 1:
        sig_positions = []
        sig_values = []
        sig_colors_scatter = []
        sig_labels = []
        
        for variant in variants:
            clinical_info = variant.get('clinical_info')
            if clinical_info and hasattr(clinical_info, 'clinical_significance'):
                sig_positions.append(variant.get('position', 0))
                # Convert significance to numeric value for plotting
                sig_score = {
                    'Pathogenic': 5,
                    'Likely pathogenic': 4,
                    'Uncertain significance': 3,
                    'Likely benign': 2,
                    'Benign': 1
                }.get(clinical_info.clinical_significance, 3)
                sig_values.append(sig_score)
                sig_colors_scatter.append(colors.get(clinical_info.clinical_significance, '#CCCCCC'))
                sig_labels.append(clinical_info.clinical_significance)
        
        if sig_positions:
            fig.add_trace(
                go.Scatter(
                    x=sig_positions,
                    y=sig_values,
                    mode='markers',
                    marker=dict(
                        size=12,
                        color=sig_colors_scatter,
                        line=dict(width=1, color='black')
                    ),
                    text=sig_labels,
                    hovertemplate="Position: %{x}<br>Significance: %{text}<extra></extra>",
                    name="Clinical Significance"
                ),
                row=2, col=2
            )
    
    fig.update_layout(
        title="Variant Analysis Dashboard",
        height=600,
        showlegend=False
    )
    
    # Update y-axis for significance plot
    fig.update_yaxes(
        title_text="Clinical Significance Score",
        tickvals=[1, 2, 3, 4, 5],
        ticktext=["Benign", "Likely Benign", "Uncertain", "Likely Pathogenic", "Pathogenic"],
        row=2, col=2
    )
    
    fig.update_yaxes(title_text="", showticklabels=False, row=2, col=1)
    fig.update_xaxes(title_text="Position", row=2, col=1)
    fig.update_xaxes(title_text="Position", row=2, col=2)
    
    return fig

def create_sequence_quality_metrics(results: Dict) -> go.Figure:
    """
    Create sequence quality and alignment metrics visualization
    """
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=("Alignment Scores", "Sequence Identity", 
                       "Length Comparison", "Quality Metrics"),
        specs=[[{"type": "bar"}, {"type": "bar"}],
               [{"type": "bar"}, {"type": "indicator"}]]
    )
    
    # Extract alignment results
    alignment_results = results.get('alignment_results', {})
    
    if not alignment_results:
        # Show empty state
        fig.add_annotation(
            text="No alignment data available",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=16)
        )
        fig.update_layout(
            title="Sequence Analysis Quality Metrics",
            height=500
        )
        return fig
    
    # Get available algorithms
    algorithms = list(alignment_results.keys())
    
    # Extract scores - handle different possible data structures
    scores = []
    identities = []
    
    for alg in algorithms:
        alg_data = alignment_results[alg]
        
        # Handle different score formats
        if isinstance(alg_data.get('score'), (int, float)):
            scores.append(alg_data['score'])
        elif hasattr(alg_data.get('score'), 'score'):
            scores.append(alg_data['score'].score)
        else:
            scores.append(0)
        
        # Handle identity
        if 'identity' in alg_data:
            identities.append(alg_data['identity'])
        else:
            # Calculate identity if not provided
            ref_seq = alg_data.get('aligned_reference', '')
            sample_seq = alg_data.get('aligned_sample', '')
            if ref_seq and sample_seq:
                identity = calculate_sequence_identity(ref_seq, sample_seq)
                identities.append(identity)
            else:
                identities.append(0)
    
    # Alignment scores
    if scores:
        fig.add_trace(
            go.Bar(
                x=algorithms,
                y=scores,
                name="Alignment Score",
                marker_color=['#FF6B6B', '#4ECDC4'][:len(algorithms)],
                text=[f"{score:.1f}" for score in scores],
                textposition='auto'
            ),
            row=1, col=1
        )
    
    # Identity percentages
    if identities:
        fig.add_trace(
            go.Bar(
                x=algorithms,
                y=identities,
                name="Identity %",
                marker_color=['#45B7D1', '#96CEB4'][:len(algorithms)],
                text=[f"{identity:.1f}%" for identity in identities],
                textposition='auto'
            ),
            row=1, col=2
        )
    
    # Length comparison
    ref_length = results.get('reference_length', 0)
    sample_length = results.get('sample_length', 0)
    
    fig.add_trace(
        go.Bar(
            x=['Reference', 'Sample'],
            y=[ref_length, sample_length],
            name="Length",
            marker_color=['#FFA07A', '#98D8C8'],
            text=[f"{ref_length} bp", f"{sample_length} bp"],
            textposition='auto'
        ),
        row=2, col=1
    )
    
    # Overall quality indicator
    quality_score = max(identities) if identities else 0
    
    # Determine quality color
    if quality_score >= 90:
        quality_color = "green"
    elif quality_score >= 70:
        quality_color = "yellow"
    else:
        quality_color = "red"
    
    fig.add_trace(
        go.Indicator(
            mode="gauge+number+delta",
            value=quality_score,
            domain={'x': [0, 1], 'y': [0, 1]},
            title={'text': "Best Identity %"},
            gauge={
                'axis': {'range': [None, 100]},
                'bar': {'color': quality_color},
                'steps': [
                    {'range': [0, 50], 'color': "lightgray"},
                    {'range': [50, 80], 'color': "yellow"},
                    {'range': [80, 100], 'color': "lightgreen"}
                ],
                'threshold': {
                    'line': {'color': "red", 'width': 4},
                    'thickness': 0.75,
                    'value': 90
                }
            }
        ),
        row=2, col=2
    )
    
    fig.update_layout(
        title="Sequence Analysis Quality Metrics",
        height=500,
        showlegend=False
    )
    
    # Update axis labels
    fig.update_yaxes(title_text="Score", row=1, col=1)
    fig.update_yaxes(title_text="Identity (%)", row=1, col=2)
    fig.update_yaxes(title_text="Length (bp)", row=2, col=1)
    
    return fig

def calculate_sequence_identity(seq1: str, seq2: str) -> float:
    """
    Calculate sequence identity percentage between two aligned sequences
    """
    if not seq1 or not seq2 or len(seq1) != len(seq2):
        return 0.0
    
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != '-' and b != '-')
    valid_positions = sum(1 for a, b in zip(seq1, seq2) if a != '-' and b != '-')
    
    if valid_positions == 0:
        return 0.0
    
    return (matches / valid_positions) * 100

def create_nucleotide_composition_chart(results: Dict) -> go.Figure:
    """
    Create nucleotide composition comparison chart
    """
    fig = go.Figure()
    
    # Get reference and sample info
    ref_info = results.get('reference_info', {})
    sample_info = results.get('sample_info', {})
    
    ref_composition = ref_info.get('nucleotide_composition', {})
    sample_composition = sample_info.get('nucleotide_composition', {})
    
    nucleotides = ['A', 'T', 'G', 'C']
    
    # Reference composition
    ref_values = [ref_composition.get(nt, 0) for nt in nucleotides]
    sample_values = [sample_composition.get(nt, 0) for nt in nucleotides]
    
    fig.add_trace(go.Bar(
        x=nucleotides,
        y=ref_values,
        name='Reference',
        marker_color='lightblue',
        text=[f"{val:.1f}%" for val in ref_values],
        textposition='auto'
    ))
    
    fig.add_trace(go.Bar(
        x=nucleotides,
        y=sample_values,
        name='Sample',
        marker_color='lightcoral',
        text=[f"{val:.1f}%" for val in sample_values],
        textposition='auto'
    ))
    
    fig.update_layout(
        title="Nucleotide Composition Comparison",
        xaxis_title="Nucleotide",
        yaxis_title="Percentage (%)",
        barmode='group',
        height=400
    )
    
    return fig

def create_variant_summary_table(variants: List[Dict]) -> pd.DataFrame:
    """
    Create a summary table of variants
    """
    if not variants:
        return pd.DataFrame(columns=['Position', 'Type', 'Variant Name', 'Clinical Significance', 'ClinVar ID'])
    
    data = []
    for variant in variants:
        clinical_info = variant.get('clinical_info')
        
        row = {
            'Position': variant.get('position', 'Unknown'),
            'Type': variant.get('type', 'Unknown'),
            'Reference': variant.get('reference', 'Unknown'),
            'Sample': variant.get('sample', 'Unknown'),
            'Variant Name': variant.get('variant_name', 'Unknown'),
            'Clinical Significance': 'Unknown',
            'ClinVar ID': 'Unknown',
            'Description': 'No information available'
        }
        
        if clinical_info:
            if hasattr(clinical_info, 'clinical_significance'):
                row['Clinical Significance'] = clinical_info.clinical_significance
            if hasattr(clinical_info, 'clinvar_id'):
                row['ClinVar ID'] = clinical_info.clinvar_id
            if hasattr(clinical_info, 'description'):
                row['Description'] = clinical_info.description
        
        data.append(row)
    
    return pd.DataFrame(data)

def display_alignment_details(alignment_data: Dict, max_length: int = 200) -> None:
    """
    Display detailed alignment information in Streamlit
    """
    if not alignment_data:
        st.warning("No alignment data available")
        return
    
    ref_seq = alignment_data.get('aligned_reference', '')
    sample_seq = alignment_data.get('aligned_sample', '')
    score = alignment_data.get('score', 0)
    identity = alignment_data.get('identity', 0)
    
    st.write(f"**Alignment Score:** {score}")
    st.write(f"**Sequence Identity:** {identity:.2f}%")
    st.write(f"**Alignment Length:** {len(ref_seq)} bp")
    
    if len(ref_seq) > max_length:
        st.warning(f"Sequence too long for display. Showing first {max_length} characters.")
        ref_seq = ref_seq[:max_length] + "..."
        sample_seq = sample_seq[:max_length] + "..."
    
    # Create alignment view
    st.write("**Alignment Preview:**")
    
    # Split into chunks for better readability
    chunk_size = 60
    for i in range(0, len(ref_seq), chunk_size):
        ref_chunk = ref_seq[i:i+chunk_size]
        sample_chunk = sample_seq[i:i+chunk_size]
        
        # Create match string
        match_string = ""
        for r, s in zip(ref_chunk, sample_chunk):
            if r == s:
                match_string += "|"
            else:
                match_string += " "
        
        st.text(f"Ref:    {ref_chunk}")
        st.text(f"        {match_string}")
        st.text(f"Sample: {sample_chunk}")
        st.text("")  # Empty line for separation