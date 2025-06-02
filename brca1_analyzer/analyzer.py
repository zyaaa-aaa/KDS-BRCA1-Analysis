from .aligner import SequenceAligner
from .clinvar import ClinVarIntegrator
from typing import Dict, List

class BRCA1Analyzer:
    """Main analyzer class combining all components"""
    
    def __init__(self):
        self.aligner = SequenceAligner()
        self.clinvar = ClinVarIntegrator()
        self.reference_sequence = self._get_reference_sequence()
    
    def _get_reference_sequence(self) -> str:
        """
        Get BRCA1 reference sequence (mock for demonstration)
        In real implementation, this would fetch from NCBI GenBank
        """
        # This is a simplified BRCA1 sequence for demonstration
        return """ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAG
TGTCCCATCTGTCTGGAGTTGATCAAGGATAAGATGGAGGGAAGATGAAATCAGTGTATGATCCAGTAGA
TGTGGCTCAATCTGAAGAAGTTAAAGAAGCTGAAAGAGCTGCAGAAGCATTATTAGAAGCACAGAAGGCT
GAACCCAAATGTAAGGATTCTGAAGATAAGGCATATAATCAGAGGAGAAAAAAGTGGAGGAAGAAGAAAT
GGAGAAAGAAGAGAAGAAAGCAGAGAAAATGGAAAGAATGCAGATGGAGAAAGAAGAAGAAGATGAAGAA"""
    
    def analyze_sequence(self, sample_sequence: str, alignment_type: str = "both") -> Dict:
        """
        Main analysis function
        """
        results = {
            'sample_length': len(sample_sequence),
            'reference_length': len(self.reference_sequence),
            'variants': [],
            'alignment_results': {}
        }
        
        # Perform alignments
        if alignment_type in ["local", "both"]:
            local_result, local_seq1, local_seq2 = self.aligner.smith_waterman(
                self.reference_sequence, sample_sequence
            )
            results['alignment_results']['local'] = {
                'result': local_result,
                'aligned_reference': local_seq1,
                'aligned_sample': local_seq2
            }
        
        if alignment_type in ["global", "both"]:
            global_result, global_seq1, global_seq2 = self.aligner.needleman_wunsch(
                self.reference_sequence, sample_sequence
            )
            results['alignment_results']['global'] = {
                'result': global_result,
                'aligned_reference': global_seq1,
                'aligned_sample': global_seq2
            }
        
        # Identify variants
        variants = self._identify_variants(
            results['alignment_results'].get('global', {}).get('aligned_reference', ''),
            results['alignment_results'].get('global', {}).get('aligned_sample', '')
        )
        
        # Classify variants using ClinVar
        for variant in variants:
            clinvar_info = self.clinvar.search_variant("BRCA1", variant['variant_name'])
            variant['clinical_info'] = clinvar_info
        
        results['variants'] = variants
        
        return results
    
    def _identify_variants(self, aligned_ref: str, aligned_sample: str) -> List[Dict]:
        """
        Identify variants from aligned sequences
        """
        variants = []
        
        if not aligned_ref or not aligned_sample:
            return variants
        
        i = 0
        while i < min(len(aligned_ref), len(aligned_sample)):
            if aligned_ref[i] != aligned_sample[i]:
                # Found a difference
                variant_type = self._classify_variant_type(
                    aligned_ref, aligned_sample, i
                )
                
                variant = {
                    'position': i + 1,
                    'reference': aligned_ref[i],
                    'sample': aligned_sample[i],
                    'type': variant_type,
                    'variant_name': f"c.{i+1}{aligned_ref[i]}>{aligned_sample[i]}"
                }
                
                variants.append(variant)
            
            i += 1
        
        return variants
    
    def _classify_variant_type(self, ref_seq: str, sample_seq: str, pos: int) -> str:
        """
        Classify the type of variant (SNV, insertion, deletion)
        """
        if ref_seq[pos] == '-':
            return "insertion"
        elif sample_seq[pos] == '-':
            return "deletion"
        else:
            return "substitution"