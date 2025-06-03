from .aligner import SequenceAligner
from .clinvar import ClinVarIntegrator
import requests
import xml.etree.ElementTree as ET
import time
from typing import Dict, List, Optional
import streamlit as st
import re
from Bio import Entrez, SeqIO
from io import StringIO

class NCBISequenceFetcher:
    """
    Real implementation for fetching reference sequences from NCBI GenBank
    """
    
    def __init__(self, email: str = None, api_key: str = None):
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        self.email = email or "your.email@example.com"  # Replace with actual email
        self.api_key = api_key
        self.sequence_cache = {}
        
        # Set Entrez email for Biopython if available
        try:
            if self.email:
                Entrez.email = self.email
            if self.api_key:
                Entrez.api_key = self.api_key
        except ImportError:
            # Biopython not available, use requests only
            pass
    
    def get_reference_sequence(self, gene_name: str = "BRCA1", 
                              sequence_type: str = "mRNA") -> str:
        """
        Get reference sequence from NCBI GenBank
        """
        cache_key = f"{gene_name}_{sequence_type}"
        
        if cache_key in self.sequence_cache:
            return self.sequence_cache[cache_key]
        
        try:
            # First, search for the gene to get accession numbers
            accession_id = self._search_gene_accession(gene_name, sequence_type)
            
            if not accession_id:
                st.warning(f"Could not find accession ID for {gene_name}")
                return ""
            
            # Fetch the actual sequence
            sequence = self._fetch_sequence_by_accession(accession_id)
            
            if sequence:
                self.sequence_cache[cache_key] = sequence
                st.success(f"Successfully fetched {gene_name} sequence from NCBI (Accession: {accession_id})")
                return sequence
            else:
                st.error(f"Could not fetch sequence for accession {accession_id}")
                return ""
                
        except Exception as e:
            st.error(f"Error fetching {gene_name} sequence: {str(e)}")
            return ""
    
    def _search_gene_accession(self, gene_name: str, sequence_type: str) -> Optional[str]:
        """
        Search for gene accession ID in NCBI
        """
        try:
            # Define search terms for BRCA1 with known accessions as fallback
            known_accessions = {
                "BRCA1": {
                    "mRNA": "NM_007294.4",  # BRCA1 mRNA reference
                    "genomic": "NC_000017.11",  # Chromosome 17 genomic
                    "cds": "NM_007294.4"
                },
                "BRCA2": {
                    "mRNA": "NM_000059.4",
                    "genomic": "NC_000013.11",
                    "cds": "NM_000059.4"
                }
            }
            
            # Try known accession first
            if gene_name.upper() in known_accessions and sequence_type in known_accessions[gene_name.upper()]:
                return known_accessions[gene_name.upper()][sequence_type]
            
            # Search terms for dynamic lookup
            search_terms = {
                "mRNA": f"{gene_name}[gene] AND homo sapiens[orgn] AND refseq[filter] AND mRNA[filter]",
                "genomic": f"{gene_name}[gene] AND homo sapiens[orgn] AND refseq[filter]",
                "cds": f"{gene_name}[gene] AND homo sapiens[orgn] AND refseq[filter] AND CDS[filter]"
            }
            
            search_term = search_terms.get(sequence_type, search_terms["mRNA"])
            
            # Use nucleotide database for DNA/RNA sequences
            database = "nucleotide"
            
            # Search using ESearch
            search_url = f"{self.base_url}esearch.fcgi"
            params = {
                'db': database,
                'term': search_term,
                'retmode': 'json',
                'retmax': '5',
                'sort': 'relevance'
            }
            
            if self.api_key:
                params['api_key'] = self.api_key
            
            response = requests.get(search_url, params=params, timeout=10)
            response.raise_for_status()
            
            data = response.json()
            
            if 'esearchresult' in data and data['esearchresult']['idlist']:
                # Get the first (most relevant) result
                gene_id = data['esearchresult']['idlist'][0]
                
                # Get accession number using ESummary
                accession = self._get_accession_from_id(gene_id, database)
                return accession
            
            return None
            
        except Exception as e:
            st.warning(f"Search failed for {gene_name}: {e}")
            return None
    
    def _get_accession_from_id(self, gene_id: str, database: str) -> Optional[str]:
        """
        Get accession number from NCBI ID using ESummary
        """
        try:
            summary_url = f"{self.base_url}esummary.fcgi"
            params = {
                'db': database,
                'id': gene_id,
                'retmode': 'json'
            }
            
            if self.api_key:
                params['api_key'] = self.api_key
            
            response = requests.get(summary_url, params=params, timeout=10)
            response.raise_for_status()
            
            data = response.json()
            
            if 'result' in data and gene_id in data['result']:
                result = data['result'][gene_id]
                
                # Try to get accession number from different fields
                accession = (result.get('accessionversion') or 
                           result.get('caption') or 
                           gene_id)
                
                return accession
                
            return gene_id  # Fallback to gene_id
            
        except Exception as e:
            st.warning(f"Failed to get accession for ID {gene_id}: {e}")
            return gene_id
    
    def _fetch_sequence_by_accession(self, accession_id: str) -> Optional[str]:
        """
        Fetch sequence using accession ID
        """
        try:
            # Method 1: Using Biopython Entrez (if available)
            try:
                handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
                record = SeqIO.read(handle, "fasta")
                handle.close()
                return str(record.seq)
            except (ImportError, NameError):
                # Biopython not available, use direct API
                pass
            except Exception:
                # Biopython failed, try direct API
                pass
            
            # Method 2: Using direct API call
            fetch_url = f"{self.base_url}efetch.fcgi"
            params = {
                'db': 'nucleotide',
                'id': accession_id,
                'rettype': 'fasta',
                'retmode': 'text'
            }
            
            if self.api_key:
                params['api_key'] = self.api_key
            
            response = requests.get(fetch_url, params=params, timeout=30)
            response.raise_for_status()
            
            # Parse FASTA format
            lines = response.text.strip().split('\n')
            if lines and lines[0].startswith('>'):
                # Remove header and join sequence lines
                sequence = ''.join(lines[1:])
                # Remove any whitespace or numbers, keep only nucleotides
                sequence = re.sub(r'[^ATCGNatcgn]', '', sequence)
                return sequence.upper()
            
            return None
            
        except Exception as e:
            st.warning(f"Failed to fetch sequence for {accession_id}: {e}")
            return None


class BRCA1Analyzer:
    """Main analyzer class combining all components with real NCBI integration"""
    
    def __init__(self, email: str = None, api_key: str = None):
        self.aligner = SequenceAligner()
        self.clinvar = ClinVarIntegrator()
        self.ncbi_fetcher = NCBISequenceFetcher(email, api_key)
        self.reference_sequence = ""
        self.reference_accession = ""
        
        # Initialize reference sequence
        self._initialize_reference_sequence()
    
    def _initialize_reference_sequence(self, gene_name: str = "BRCA1") -> None:
        """
        Initialize reference sequence from NCBI
        """
        try:
            with st.spinner(f"Fetching {gene_name} reference sequence from NCBI..."):
                self.reference_sequence = self.ncbi_fetcher.get_reference_sequence(
                    gene_name=gene_name, 
                    sequence_type="mRNA"
                )
                
                if self.reference_sequence:
                    st.info(f"✅ Reference sequence loaded: {len(self.reference_sequence)} nucleotides")
                else:
                    st.error(f"❌ Failed to load reference sequence for {gene_name}")
                    
        except Exception as e:
            st.error(f"Error initializing reference sequence: {e}")
    
    def set_custom_gene(self, gene_name: str) -> bool:
        """
        Set a different gene for analysis
        """
        try:
            with st.spinner(f"Fetching {gene_name} reference sequence..."):
                new_sequence = self.ncbi_fetcher.get_reference_sequence(
                    gene_name=gene_name, 
                    sequence_type="mRNA"
                )
                
                if new_sequence:
                    self.reference_sequence = new_sequence
                    st.success(f"✅ Successfully loaded {gene_name} reference sequence")
                    return True
                else:
                    st.error(f"❌ Failed to load reference sequence for {gene_name}")
                    return False
                    
        except Exception as e:
            st.error(f"Error setting custom gene: {e}")
            return False
    
    def get_reference_info(self) -> Dict:
        """
        Get information about the current reference sequence
        """
        return {
            'length': len(self.reference_sequence),
            'accession': self.reference_accession,
            'gc_content': self._calculate_gc_content(self.reference_sequence),
            'nucleotide_composition': self._get_nucleotide_composition(self.reference_sequence)
        }
    
    def _calculate_gc_content(self, sequence: str) -> float:
        """Calculate GC content of sequence"""
        if not sequence:
            return 0.0
        gc_count = sequence.count('G') + sequence.count('C')
        return (gc_count / len(sequence)) * 100
    
    def _get_nucleotide_composition(self, sequence: str) -> Dict:
        """Get nucleotide composition"""
        if not sequence:
            return {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0}
        
        total = len(sequence)
        return {
            'A': round((sequence.count('A') / total) * 100, 2),
            'T': round((sequence.count('T') / total) * 100, 2),
            'G': round((sequence.count('G') / total) * 100, 2),
            'C': round((sequence.count('C') / total) * 100, 2),
            'N': round((sequence.count('N') / total) * 100, 2)
        }
    
    def analyze_sequence(self, sample_sequence: str, alignment_type: str = "both") -> Dict:
        """
        Main analysis function using real reference sequence
        """
        if not self.reference_sequence:
            st.error("No reference sequence available. Please initialize first.")
            return {}
        
        # Clean input sequence
        sample_sequence = re.sub(r'[^ATCGNatcgn]', '', sample_sequence.upper())
        
        results = {
            'sample_length': len(sample_sequence),
            'reference_length': len(self.reference_sequence),
            'reference_info': self.get_reference_info(),
            'sample_info': {
                'gc_content': self._calculate_gc_content(sample_sequence),
                'nucleotide_composition': self._get_nucleotide_composition(sample_sequence)
            },
            'variants': [],
            'alignment_results': {}
        }
        
        # Perform alignments
        try:
            if alignment_type in ["local", "both"]:
                local_result, local_seq1, local_seq2 = self.aligner.smith_waterman(
                    self.reference_sequence, sample_sequence
                )
                results['alignment_results']['local'] = {
                    'score': local_result,
                    'aligned_reference': local_seq1,
                    'aligned_sample': local_seq2,
                    'identity': self._calculate_identity(local_seq1, local_seq2)
                }
            
            if alignment_type in ["global", "both"]:
                global_result, global_seq1, global_seq2 = self.aligner.needleman_wunsch(
                    self.reference_sequence, sample_sequence
                )
                results['alignment_results']['global'] = {
                    'score': global_result,
                    'aligned_reference': global_seq1,
                    'aligned_sample': global_seq2,
                    'identity': self._calculate_identity(global_seq1, global_seq2)
                }
                
        except Exception as e:
            st.error(f"Alignment failed: {e}")
            return results
        
        # Identify variants
        aligned_ref = results['alignment_results'].get('global', {}).get('aligned_reference', '')
        aligned_sample = results['alignment_results'].get('global', {}).get('aligned_sample', '')
        
        if aligned_ref and aligned_sample:
            variants = self._identify_variants(aligned_ref, aligned_sample)
            
            # Classify variants using ClinVar (with progress bar for multiple variants)
            if variants:
                progress_bar = st.progress(0)
                for i, variant in enumerate(variants):
                    try:
                        clinvar_info = self.clinvar.search_variant("BRCA1", variant['variant_name'])
                        variant['clinical_info'] = clinvar_info
                        progress_bar.progress((i + 1) / len(variants))
                        time.sleep(0.1)  # Rate limiting
                    except Exception as e:
                        st.warning(f"ClinVar lookup failed for {variant['variant_name']}: {e}")
                        variant['clinical_info'] = None
                
                progress_bar.empty()
            
            results['variants'] = variants
        
        return results
    
    def _calculate_identity(self, seq1: str, seq2: str) -> float:
        """Calculate sequence identity percentage"""
        if not seq1 or not seq2 or len(seq1) != len(seq2):
            return 0.0
        
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != '-' and b != '-')
        valid_positions = sum(1 for a, b in zip(seq1, seq2) if a != '-' and b != '-')
        
        if valid_positions == 0:
            return 0.0
        
        return (matches / valid_positions) * 100
    
    def _identify_variants(self, aligned_ref: str, aligned_sample: str) -> List[Dict]:
        """
        Identify variants from aligned sequences with enhanced detection
        """
        variants = []
        
        if not aligned_ref or not aligned_sample:
            return variants
        
        i = 0
        while i < min(len(aligned_ref), len(aligned_sample)):
            if aligned_ref[i] != aligned_sample[i]:
                # Detect complex variants (insertions/deletions)
                variant = self._analyze_variant_at_position(aligned_ref, aligned_sample, i)
                if variant:
                    variants.append(variant)
                    # Skip processed positions
                    i += variant.get('length', 1)
                else:
                    i += 1
            else:
                i += 1
        
        return variants
    
    def _analyze_variant_at_position(self, ref_seq: str, sample_seq: str, pos: int) -> Optional[Dict]:
        """
        Analyze variant at specific position with complex variant detection
        """
        try:
            # Simple substitution
            if ref_seq[pos] != '-' and sample_seq[pos] != '-':
                return {
                    'position': pos + 1,
                    'reference': ref_seq[pos],
                    'sample': sample_seq[pos],
                    'type': 'substitution',
                    'variant_name': f"c.{pos+1}{ref_seq[pos]}>{sample_seq[pos]}",
                    'length': 1
                }
            
            # Insertion in sample (gap in reference)
            elif ref_seq[pos] == '-' and sample_seq[pos] != '-':
                # Find extent of insertion
                end_pos = pos
                while (end_pos < len(ref_seq) and end_pos < len(sample_seq) and 
                       ref_seq[end_pos] == '-' and sample_seq[end_pos] != '-'):
                    end_pos += 1
                
                inserted_seq = sample_seq[pos:end_pos]
                return {
                    'position': pos + 1,
                    'reference': '-',
                    'sample': inserted_seq,
                    'type': 'insertion',
                    'variant_name': f"c.{pos+1}ins{inserted_seq}",
                    'length': end_pos - pos
                }
            
            # Deletion in sample (gap in sample)
            elif ref_seq[pos] != '-' and sample_seq[pos] == '-':
                # Find extent of deletion
                end_pos = pos
                while (end_pos < len(ref_seq) and end_pos < len(sample_seq) and 
                       ref_seq[end_pos] != '-' and sample_seq[end_pos] == '-'):
                    end_pos += 1
                
                deleted_seq = ref_seq[pos:end_pos]
                return {
                    'position': pos + 1,
                    'reference': deleted_seq,
                    'sample': '-',
                    'type': 'deletion',
                    'variant_name': f"c.{pos+1}del{deleted_seq}",
                    'length': end_pos - pos
                }
            
            return None
            
        except Exception as e:
            st.warning(f"Error analyzing variant at position {pos}: {e}")
            return None
    
    def export_results(self, results: Dict, filename: str = None) -> str:
        """
        Export analysis results to text format
        """
        if not filename:
            filename = f"brca_analysis_{int(time.time())}.txt"
        
        try:
            output = []
            output.append("BRCA1 Sequence Analysis Results")
            output.append("=" * 40)
            output.append(f"Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}")
            output.append(f"Reference Length: {results.get('reference_length', 0)} bp")
            output.append(f"Sample Length: {results.get('sample_length', 0)} bp")
            
            # Reference info
            ref_info = results.get('reference_info', {})
            output.append(f"\nReference Sequence Info:")
            output.append(f"  GC Content: {ref_info.get('gc_content', 0):.2f}%")
            
            # Alignment results
            if 'alignment_results' in results:
                output.append(f"\nAlignment Results:")
                for align_type, align_data in results['alignment_results'].items():
                    output.append(f"  {align_type.title()} Alignment:")
                    output.append(f"    Score: {align_data.get('score', 0)}")
                    output.append(f"    Identity: {align_data.get('identity', 0):.2f}%")
            
            # Variants
            variants = results.get('variants', [])
            output.append(f"\nVariants Found: {len(variants)}")
            
            for i, variant in enumerate(variants, 1):
                output.append(f"\n  Variant {i}:")
                output.append(f"    Position: {variant.get('position', 'Unknown')}")
                output.append(f"    Type: {variant.get('type', 'Unknown')}")
                output.append(f"    Variant Name: {variant.get('variant_name', 'Unknown')}")
                output.append(f"    Reference: {variant.get('reference', 'Unknown')}")
                output.append(f"    Sample: {variant.get('sample', 'Unknown')}")
                
                clinical_info = variant.get('clinical_info')
                if clinical_info:
                    output.append(f"    Clinical Significance: {clinical_info.clinical_significance}")
                    output.append(f"    ClinVar ID: {clinical_info.clinvar_id}")
            
            result_text = '\n'.join(output)
            
            # In a real application, you would save to file
            # For Streamlit, we'll return the text for download
            return result_text
            
        except Exception as e:
            st.error(f"Export failed: {e}")
            return f"Export failed: {e}"