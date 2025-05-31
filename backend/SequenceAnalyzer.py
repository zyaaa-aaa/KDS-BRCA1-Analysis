#!/usr/bin/env python3
# BRCA1 Sequence Analyzer
# A tool for analyzing BRCA1 gene sequences and identifying variants

import os
import argparse
import time
import json
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align import PairwiseAligner
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
# Removed Bio.Align.Applications import to avoid deprecation warning
import subprocess
import tempfile
import requests

class SequenceAnalyzer:
    """
    A class for analyzing BRCA1 gene sequences and identifying variants.
    """
    
    def __init__(self, reference="NM_007294.4", email=None):
        """
        Initialize the SequenceAnalyzer with a reference sequence.
        
        Parameters:
        -----------
        reference : str
            The accession number of the reference BRCA1 sequence.
        email : str
            Email to use for NCBI Entrez queries. Required by NCBI.
        """
        self.reference_id = reference
        self.reference_seq = None
        self.variants_db = None
        self.clinvar_data = None
        
        # Set email for NCBI Entrez
        if email:
            Entrez.email = email
        else:
            Entrez.email = "your_email@example.com"  # Replace with actual email
        
        # Load reference sequence
        self._load_reference()
        
    def _load_reference(self):
        """
        Load the reference BRCA1 sequence from NCBI.
        """
        try:
            handle = Entrez.efetch(db="nucleotide", id=self.reference_id, rettype="fasta", retmode="text")
            self.reference_seq = SeqIO.read(handle, "fasta")
            handle.close()
            print(f"Successfully loaded reference sequence: {self.reference_id}")
            print(f"Sequence length: {len(self.reference_seq.seq)} bp")
        except Exception as e:
            print(f"Error loading reference sequence: {e}")
            raise
    
    def load_clinvar_data(self, gene="BRCA1"):
        """
        Load ClinVar data for BRCA1 variants.
        
        Parameters:
        -----------
        gene : str
            Gene symbol to search for in ClinVar.
        """
        try:
            # ClinVar API endpoint
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={gene}[gene]&retmax=1000&retmode=json"
            response = requests.get(url)
            data = response.json()
            
            # Get IDs of ClinVar entries
            ids = data['esearchresult']['idlist']
            
            # Fetch detailed records for each ID
            clinvar_records = []
            for batch in [ids[i:i + 50] for i in range(0, len(ids), 50)]:
                id_str = ",".join(batch)
                fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={id_str}&retmode=json"
                batch_response = requests.get(fetch_url)
                batch_data = batch_response.json()
                
                # The structure of batch_data can vary; handle both common formats
                if 'result' in batch_data:
                    result = batch_data['result']
                    # Some APIs return a 'result' key with an object containing IDs as keys
                    if isinstance(result, dict) and 'uids' not in result:
                        for id in batch:
                            if id in result:
                                clinvar_records.append(result[id])
                    # Others might have a 'uids' list within the result
                    elif isinstance(result, dict) and 'uids' in result:
                        for uid in result['uids']:
                            if uid in result:
                                clinvar_records.append(result[uid])
                    # Or directly a list of records
                    elif isinstance(result, list):
                        clinvar_records.extend(result)
            
            # Convert to DataFrame and handle potential empty results
            if clinvar_records:
                self.clinvar_data = pd.DataFrame(clinvar_records)
                print(f"Loaded {len(clinvar_records)} ClinVar records for {gene}")
                # Print a sample record to see structure
                if len(clinvar_records) > 0:
                    print("Sample ClinVar record keys:", list(clinvar_records[0].keys()))
            else:
                print("No ClinVar records found for query")
                self.clinvar_data = pd.DataFrame()
            
        except Exception as e:
            print(f"Error loading ClinVar data: {e}")
            self.clinvar_data = pd.DataFrame()
    
    def fetch_population_variants(self, population=None, max_variants=100):
        """
        Fetch BRCA1 variants from a specific population.
        
        Parameters:
        -----------
        population : str
            Population identifier (e.g., "European", "African", "Asian").
        max_variants : int
            Maximum number of variants to retrieve.
            
        Returns:
        --------
        list
            List of SeqRecord objects containing variant sequences.
        """
        query = "BRCA1[Gene]"
        if population:
            query += f" AND {population}[Population]"
        
        try:
            # Search for variants
            handle = Entrez.esearch(db="nuccore", term=query, retmax=max_variants)
            search_results = Entrez.read(handle)
            handle.close()
            
            # Get the IDs of variants
            id_list = search_results["IdList"]
            
            if not id_list:
                print(f"No variants found for query: {query}")
                return []
            
            # Fetch the sequences
            handle = Entrez.efetch(db="nuccore", id=",".join(id_list), rettype="fasta", retmode="text")
            variants = list(SeqIO.parse(handle, "fasta"))
            handle.close()
            
            print(f"Fetched {len(variants)} BRCA1 variants")
            return variants
            
        except Exception as e:
            print(f"Error fetching population variants: {e}")
            return []
    
    def align_sequences(self, seq_record, alignment_type="global"):
        """
        Align a sequence to the reference BRCA1 sequence.
        
        Parameters:
        -----------
        seq_record : SeqRecord
            The sequence to align to the reference.
        alignment_type : str
            Type of alignment: "global" or "local".
            
        Returns:
        --------
        alignment
            The alignment object.
        """
        if self.reference_seq is None:
            raise ValueError("Reference sequence not loaded")
        
        ref_seq = str(self.reference_seq.seq)
        query_seq = str(seq_record.seq)
        
        print(f"Aligning sequences: reference length={len(ref_seq)}, query length={len(query_seq)}")
        
        # For long sequences, we'll use a simpler approach that's less likely to overflow
        # Get first 100 amino acids for debugging
        if len(ref_seq) > 5000 or len(query_seq) > 5000:
            print("Sequences are very long, using simplified alignment")
            
            # Create a simple alignment by comparing sequences position by position
            ref_aligned = ""
            query_aligned = ""
            score = 0
            
            # Determine length for alignment (use the shorter sequence)
            max_length = min(len(ref_seq), len(query_seq))
            
            # Compare sequences position by position
            for i in range(max_length):
                ref_base = ref_seq[i]
                query_base = query_seq[i]
                
                ref_aligned += ref_base
                query_aligned += query_base
                
                # Simple scoring: +1 for match, -1 for mismatch
                if ref_base == query_base:
                    score += 1
                else:
                    score -= 1
            
            # Handle remaining bases in the longer sequence
            if len(ref_seq) > max_length:
                ref_aligned += ref_seq[max_length:]
                query_aligned += '-' * (len(ref_seq) - max_length)
            elif len(query_seq) > max_length:
                query_aligned += query_seq[max_length:]
                ref_aligned += '-' * (len(query_seq) - max_length)
            
            return (ref_aligned, query_aligned, score)
        
        # For shorter sequences, use the PairwiseAligner
        try:
            # Create a pairwise aligner
            aligner = PairwiseAligner()
            
            # Limit the number of alignments
            aligner.match_score = 2
            aligner.mismatch_score = -1
            aligner.open_gap_score = -0.5
            aligner.extend_gap_score = -0.1
            aligner.target_internal_gap_score = -0.5
            aligner.query_internal_gap_score = -0.5
            aligner.target_end_gap_score = -0.1
            aligner.query_end_gap_score = -0.1
            
            # Set alignment mode
            if alignment_type == "global":
                aligner.mode = 'global'
            else:
                aligner.mode = 'local'
            
            # Get only the best alignment
            print("Getting one optimal alignment...")
            alignment = aligner.align(ref_seq, query_seq)
            best_alignment = next(alignment)
            
            # Convert to format similar to old pairwise2 for compatibility
            aligned_ref, aligned_query = best_alignment.aligned
            
            # Convert to strings with gaps
            ref_aligned = ""
            query_aligned = ""
            ref_pos = 0
            query_pos = 0
            
            for ref_block, query_block in zip(aligned_ref, aligned_query):
                # Add any gaps before this block
                while ref_pos < ref_block[0]:
                    ref_aligned += ref_seq[ref_pos]
                    query_aligned += '-'
                    ref_pos += 1
                
                while query_pos < query_block[0]:
                    ref_aligned += '-'
                    query_aligned += query_seq[query_pos]
                    query_pos += 1
                
                # Add the aligned block
                ref_part = ref_seq[ref_block[0]:ref_block[1]]
                query_part = query_seq[query_block[0]:query_block[1]]
                ref_aligned += ref_part
                query_aligned += query_part
                
                ref_pos = ref_block[1]
                query_pos = query_block[1]
            
            # Add any remaining sequence
            while ref_pos < len(ref_seq):
                ref_aligned += ref_seq[ref_pos]
                query_aligned += '-'
                ref_pos += 1
            
            while query_pos < len(query_seq):
                ref_aligned += '-'
                query_aligned += query_seq[query_pos]
                query_pos += 1
            
            # Return in a format similar to pairwise2 for compatibility
            return (ref_aligned, query_aligned, best_alignment.score)
        
        except Exception as e:
            print(f"Error in sequence alignment: {e}")
            print("Falling back to simplified alignment method")
            
            # Use the same simplified alignment approach as above if PairwiseAligner fails
            ref_aligned = ""
            query_aligned = ""
            score = 0
            
            max_length = min(len(ref_seq), len(query_seq))
            
            for i in range(max_length):
                ref_base = ref_seq[i]
                query_base = query_seq[i]
                
                ref_aligned += ref_base
                query_aligned += query_base
                
                if ref_base == query_base:
                    score += 1
                else:
                    score -= 1
            
            if len(ref_seq) > max_length:
                ref_aligned += ref_seq[max_length:]
                query_aligned += '-' * (len(ref_seq) - max_length)
            elif len(query_seq) > max_length:
                query_aligned += query_seq[max_length:]
                ref_aligned += '-' * (len(query_seq) - max_length)
            
            return (ref_aligned, query_aligned, score)
    
    def identify_variants(self, alignment):
        """
        Identify variants in an aligned sequence.
        
        Parameters:
        -----------
        alignment : alignment
            Pairwise alignment between reference and query sequence.
            
        Returns:
        --------
        list
            List of dictionaries containing variant information.
        """
        if not alignment:
            return []
        
        ref_aligned = alignment[0]
        query_aligned = alignment[1]
        
        variants = []
        position = 0
        
        for i in range(len(ref_aligned)):
            ref_base = ref_aligned[i]
            query_base = query_aligned[i]
            
            # Skip gaps in reference sequence
            if ref_base != '-':
                position += 1  # Fixed a typo: position += a1 -> position += 1
                
                # Check for mismatch (SNP) or indel
                if ref_base != query_base:
                    variant_type = "SNP" if query_base != '-' and ref_base != '-' else "indel"
                    variant = {
                        "position": position,
                        "reference": ref_base,
                        "variant": query_base,
                        "type": variant_type,
                        "notation": f"c.{position}{ref_base}>{query_base}" if variant_type == "SNP" else f"c.{position}del" if query_base == '-' else f"c.{position}ins{query_base}"
                    }
                    variants.append(variant)
        
        return variants
    
    def annotate_variants(self, variants):
        """
        Annotate variants with information from ClinVar.
        
        Parameters:
        -----------
        variants : list
            List of variant dictionaries.
            
        Returns:
        --------
        list
            Annotated variant dictionaries.
        """
        if self.clinvar_data is None:
            print("ClinVar data not loaded. Cannot annotate variants.")
            return variants
        
        # First, let's inspect the actual structure of the ClinVar data
        if len(variants) > 0 and not self.clinvar_data.empty:
            print("ClinVar data columns:", self.clinvar_data.columns.tolist())
            
            for variant in variants:
                # Initialize with default values
                variant["clinvar_id"] = "Not found"
                variant["clinical_significance"] = "Unknown"
                variant["review_status"] = "N/A"
                
                # Try a more general search approach - look for position in title
                position_str = str(variant["position"])
                
                # Search through available columns - adjust based on actual data structure
                if 'title' in self.clinvar_data.columns:
                    matches = self.clinvar_data[self.clinvar_data['title'].astype(str).str.contains(position_str, na=False)]
                    
                    if not matches.empty:
                        variant["clinvar_id"] = str(matches.iloc[0].get('uid', 'Unknown'))
                        
                        # Get clinical significance if available
                        if 'clinical_significance' in self.clinvar_data.columns:
                            variant["clinical_significance"] = str(matches.iloc[0]['clinical_significance'])
                        
                        # Get review status if available
                        if 'review_status' in self.clinvar_data.columns:
                            variant["review_status"] = str(matches.iloc[0]['review_status'])
        else:
            print("Warning: No variants found or ClinVar data is empty")
        
        return variants
    
    def analyze_sample(self, sample_file, output_file=None):
        """
        Analyze a sample sequence file and identify variants.
        
        Parameters:
        -----------
        sample_file : str
            Path to the FASTA file containing the sample sequence.
        output_file : str
            Path to the output file to save results.
            
        Returns:
        --------
        dict
            Analysis results including variants.
        """
        start_time = time.time()
        
        # Load sample sequence
        try:
            sample_seq = SeqIO.read(sample_file, "fasta")
        except Exception as e:
            print(f"Error reading sample file: {e}")
            return None
        
        # Align sample to reference
        alignment = self.align_sequences(sample_seq, "global")
        
        # Identify variants
        variants = self.identify_variants(alignment)
        
        # Annotate variants if ClinVar data is available
        if self.clinvar_data is not None:
            variants = self.annotate_variants(variants)
        
        # Calculate elapsed time
        elapsed_time = time.time() - start_time
        
        # Prepare results
        results = {
            "sample_id": sample_seq.id,
            "sample_length": len(sample_seq.seq),
            "reference_id": self.reference_seq.id,
            "alignment_score": alignment[2],
            "variants": variants,
            "processing_time": elapsed_time,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
        }
        
        # Save results to file if specified
        if output_file:
            with open(output_file, 'w') as f:
                json.dump(results, f, indent=4)
            print(f"Results saved to {output_file}")
        
        return results
    
    def perform_multiple_alignment(self, sequences, output_file=None):
        """
        Perform multiple sequence alignment on a set of BRCA1 sequences.
        
        Parameters:
        -----------
        sequences : list
            List of SeqRecord objects.
        output_file : str
            Path to the output file to save the alignment.
            
        Returns:
        --------
        MultipleSeqAlignment
            The multiple sequence alignment object.
        """
        if not sequences:
            print("No sequences provided for alignment")
            return None
        
        # Include reference sequence in the alignment
        all_sequences = [self.reference_seq] + sequences
        
        # Create a temporary file for the input sequences
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as temp_in:
            SeqIO.write(all_sequences, temp_in, "fasta")
            temp_in_name = temp_in.name
        
        # Create a temporary file for the output alignment
        temp_out_name = tempfile.mktemp(suffix=".fasta")
        
        try:
            # Run MUSCLE using subprocess directly instead of MuscleCommandline
            # This avoids the Bio.Application deprecation warning
            muscle_cmd = ["muscle", "-in", temp_in_name, "-out", temp_out_name]
            subprocess.run(muscle_cmd, check=True)
            
            # Read the alignment
            alignment = AlignIO.read(temp_out_name, "fasta")
            
            # Save to output file if specified
            if output_file:
                AlignIO.write(alignment, output_file, "fasta")
                print(f"Multiple alignment saved to {output_file}")
            
            return alignment
            
        except Exception as e:
            print(f"Error performing multiple alignment: {e}")
            return None
            
        finally:
            # Clean up temporary files
            os.remove(temp_in_name)
            if os.path.exists(temp_out_name):
                os.remove(temp_out_name)
    
    def create_phylogenetic_tree(self, alignment, output_file=None):
        """
        Create a phylogenetic tree from a multiple sequence alignment.
        
        Parameters:
        -----------
        alignment : MultipleSeqAlignment
            The multiple sequence alignment.
        output_file : str
            Path to the output file to save the tree (Newick format).
            
        Returns:
        --------
        Tree
            The phylogenetic tree object.
        """
        if not alignment:
            print("No alignment provided")
            return None
        
        try:
            # Calculate distance matrix
            calculator = DistanceCalculator('identity')
            distance_matrix = calculator.get_distance(alignment)
            
            # Construct tree using UPGMA algorithm
            constructor = DistanceTreeConstructor()
            tree = constructor.upgma(distance_matrix)
            
            # Save to output file if specified
            if output_file:
                Phylo.write(tree, output_file, 'newick')
                print(f"Phylogenetic tree saved to {output_file}")
            
            return tree
            
        except Exception as e:
            print(f"Error creating phylogenetic tree: {e}")
            return None
    
    def visualize_alignment(self, alignment, output_file=None):
        """
        Visualize a pairwise alignment as a dot plot.
        
        Parameters:
        -----------
        alignment : alignment
            Pairwise alignment between reference and query sequence.
        output_file : str
            Path to the output file to save the visualization.
        """
        if not alignment:
            print("No alignment provided")
            return
        
        ref_aligned = alignment[0]
        query_aligned = alignment[1]
        
        # Remove gaps to get original sequences
        ref_seq = ''.join([b for b in ref_aligned if b != '-'])
        query_seq = ''.join([b for b in query_aligned if b != '-'])
        
        # Create a dot plot matrix
        dot_matrix = np.zeros((len(ref_seq), len(query_seq)))
        
        # Fill the matrix
        for i in range(len(ref_seq)):
            for j in range(len(query_seq)):
                if ref_seq[i] == query_seq[j]:
                    dot_matrix[i, j] = 1
        
        # Plot the matrix
        plt.figure(figsize=(10, 10))
        plt.imshow(dot_matrix, cmap='Blues', interpolation='nearest')
        plt.xlabel('Query Sequence')
        plt.ylabel('Reference Sequence')
        plt.title('Sequence Alignment Dot Plot')
        
        # Save to output file if specified
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Dot plot saved to {output_file}")
        
        plt.close()

    def visualize_variants_distribution(self, variants, output_file=None):
        """
        Visualize the distribution of variants along the reference sequence.
        
        Parameters:
        -----------
        variants : list
            List of variant dictionaries.
        output_file : str
            Path to the output file to save the visualization.
        """
        if not variants:
            print("No variants provided")
            return
        
        # Extract positions and types
        positions = [v['position'] for v in variants]
        types = [v['type'] for v in variants]
        
        # Create colors based on variant types
        colors = ['red' if t == 'SNP' else 'blue' for t in types]
        
        # Plot the distribution
        plt.figure(figsize=(12, 6))
        plt.scatter(positions, [1] * len(positions), c=colors, s=50, alpha=0.7)
        plt.xlabel('Position on Reference Sequence')
        plt.yticks([])
        plt.title('Variant Distribution')
        
        # Add legend
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='SNP'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label='Indel')
        ]
        plt.legend(handles=legend_elements)
        
        # Add reference sequence length
        plt.axvline(x=len(self.reference_seq.seq), color='gray', linestyle='--')
        plt.text(len(self.reference_seq.seq) + 100, 1, 'End of Reference', rotation=90)
        
        # Save to output file if specified
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Variant distribution plot saved to {output_file}")
        
        plt.close()


def main():
    """
    Main function to run the BRCA1 sequence analyzer.
    """
    parser = argparse.ArgumentParser(description='BRCA1 Sequence Analyzer')
    parser.add_argument('--reference', default='NM_007294.4', help='Reference sequence accession')
    parser.add_argument('--email', required=True, help='Email for NCBI Entrez queries')
    parser.add_argument('--input', required=True, help='Input FASTA file with sample sequence')
    parser.add_argument('--output', default='results.json', help='Output file for results')
    parser.add_argument('--align', choices=['global', 'local'], default='global', help='Alignment type')
    parser.add_argument('--variants', choices=['clinvar', 'none'], default='clinvar', help='Variant annotation source')
    parser.add_argument('--visualize', action='store_true', help='Generate visualizations')
    parser.add_argument('--population', help='Fetch variants from specific population')
    
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = SequenceAnalyzer(reference=args.reference, email=args.email)
    
    # Load ClinVar data if requested
    if args.variants == 'clinvar':
        analyzer.load_clinvar_data()
    
    # Analyze sample
    results = analyzer.analyze_sample(args.input, args.output)
    
    # Generate visualizations if requested
    if args.visualize and results:
        # Align sample for visualization
        sample_seq = SeqIO.read(args.input, "fasta")
        alignment = analyzer.align_sequences(sample_seq, args.align)
        
        # Generate dot plot
        dot_plot_file = args.output.replace('.json', '_dotplot.png')
        analyzer.visualize_alignment(alignment, dot_plot_file)
        
        # Generate variant distribution plot
        var_dist_file = args.output.replace('.json', '_variants.png')
        analyzer.visualize_variants_distribution(results['variants'], var_dist_file)
    
    # Fetch population variants if requested
    if args.population:
        pop_variants = analyzer.fetch_population_variants(args.population)
        if pop_variants:
            # Perform multiple alignment
            msa_file = args.output.replace('.json', '_msa.fasta')
            msa = analyzer.perform_multiple_alignment(pop_variants, msa_file)
            
            # Create phylogenetic tree
            tree_file = args.output.replace('.json', '_tree.newick')
            analyzer.create_phylogenetic_tree(msa, tree_file)

if __name__ == "__main__":
    main()