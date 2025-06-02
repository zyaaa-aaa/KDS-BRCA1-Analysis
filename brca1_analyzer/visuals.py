import plotly.graph_objects as go
import pandas as pd
from typing import Dict, List
from plotly.subplots import make_subplots

def create_alignment_visualization(alignment_data: Dict, max_display_length: int = 100) -> go.Figure:
    """
    Create interactive alignment visualization
    """
    ref_seq = alignment_data['aligned_reference']
    sample_seq = alignment_data['aligned_sample']
    
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
        v_type = variant['type']
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
        if clinical_info:
            sig = clinical_info.clinical_significance
            significance_counts[sig] = significance_counts.get(sig, 0) + 1
    
    if significance_counts:
        colors = {
            'Pathogenic': '#FF4444',
            'Likely pathogenic': '#FF8800',
            'Uncertain significance': '#FFDD00',
            'Likely benign': '#88DD00',
            'Benign': '#44DD44'
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
    positions = [v['position'] for v in variants]
    types = [v['type'] for v in variants]
    
    fig.add_trace(
        go.Scatter(
            x=positions,
            y=[1]*len(positions),
            mode='markers',
            marker=dict(
                size=10,
                color=[hash(t) for t in types],
                colorscale='viridis',
                showscale=False
            ),
            text=[f"{v['type']}: {v['variant_name']}" for v in variants],
            hovertemplate="%{text}<extra></extra>",
            name="Positions"
        ),
        row=2, col=1
    )
    
    # Clinical significance by position
    if significance_counts:
        sig_positions = []
        sig_values = []
        sig_colors_scatter = []
        
        for variant in variants:
            clinical_info = variant.get('clinical_info')
            if clinical_info:
                sig_positions.append(variant['position'])
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
                text=[v.get('clinical_info').clinical_significance if v.get('clinical_info') else 'Unknown' 
                      for v in variants if v.get('clinical_info')],
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
                       "Gap Distribution", "Quality Metrics"),
        specs=[[{"type": "bar"}, {"type": "bar"}],
               [{"type": "bar"}, {"type": "indicator"}]]
    )
    
    # Extract alignment results
    alignment_results = results.get('alignment_results', {})
    
    # Alignment scores
    if alignment_results:
        algorithms = list(alignment_results.keys())
        scores = [alignment_results[alg]['result'].score for alg in algorithms]
        
        fig.add_trace(
            go.Bar(
                x=algorithms,
                y=scores,
                name="Alignment Score",
                marker_color=['#FF6B6B', '#4ECDC4'][:len(algorithms)]
            ),
            row=1, col=1
        )
        
        # Identity percentages
        identities = [alignment_results[alg]['result'].identity for alg in algorithms]
        
        fig.add_trace(
            go.Bar(
                x=algorithms,
                y=identities,
                name="Identity %",
                marker_color=['#45B7D1', '#96CEB4'][:len(algorithms)]
            ),
            row=1, col=2
        )
        
        # Gap distribution
        gaps = [alignment_results[alg]['result'].gaps for alg in algorithms]
        
        fig.add_trace(
            go.Bar(
                x=algorithms,
                y=gaps,
                name="Gaps",
                marker_color=['#FFA07A', '#98D8C8'][:len(algorithms)]
            ),
            row=2, col=1
        )
        
        # Overall quality indicator (using global alignment if available)
        quality_score = identities[0] if identities else 0
        
        fig.add_trace(
            go.Indicator(
                mode="gauge+number+delta",
                value=quality_score,
                domain={'x': [0, 1], 'y': [0, 1]},
                title={'text': "Sequence Quality"},
                gauge={
                    'axis': {'range': [None, 100]},
                    'bar': {'color': "darkblue"},
                    'steps': [
                        {'range': [0, 50], 'color': "lightgray"},
                        {'range': [50, 80], 'color': "yellow"},
                        {'range': [80, 100], 'color': "green"}
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
    
    return fig
