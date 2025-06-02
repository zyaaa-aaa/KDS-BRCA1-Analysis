from dataclasses import dataclass
from typing import Tuple
import numpy as np

@dataclass
class AlignmentResult:
    """Data class for storing alignment results"""
    score: float
    alignment_length: int
    identity: float
    gaps: int
    query_start: int
    query_end: int
    target_start: int
    target_end: int

class SequenceAligner:
    """Implementation of Smith-Waterman and Needleman-Wunsch algorithms"""
    
    def __init__(self, match_score=2, mismatch_score=-1, gap_penalty=-1):
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_penalty = gap_penalty
    
    def smith_waterman(self, seq1: str, seq2: str) -> Tuple[AlignmentResult, str, str]:
        """
        Local alignment using Smith-Waterman algorithm
        """
        m, n = len(seq1), len(seq2)
        
        # Initialize scoring matrix
        score_matrix = np.zeros((m + 1, n + 1))
        
        # Fill scoring matrix
        max_score = 0
        max_pos = (0, 0)
        
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match = score_matrix[i-1][j-1] + (
                    self.match_score if seq1[i-1] == seq2[j-1] else self.mismatch_score
                )
                delete = score_matrix[i-1][j] + self.gap_penalty
                insert = score_matrix[i][j-1] + self.gap_penalty
                
                score_matrix[i][j] = max(0, match, delete, insert)
                
                if score_matrix[i][j] > max_score:
                    max_score = score_matrix[i][j]
                    max_pos = (i, j)
        
        # Traceback
        aligned_seq1, aligned_seq2 = self._traceback_sw(
            score_matrix, seq1, seq2, max_pos
        )
        
        # Calculate alignment statistics
        identity = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-')
        gaps = aligned_seq1.count('-') + aligned_seq2.count('-')
        
        result = AlignmentResult(
            score=max_score,
            alignment_length=len(aligned_seq1),
            identity=identity / len(aligned_seq1) * 100,
            gaps=gaps,
            query_start=0,
            query_end=len(aligned_seq1),
            target_start=0,
            target_end=len(aligned_seq2)
        )
        
        return result, aligned_seq1, aligned_seq2
    
    def needleman_wunsch(self, seq1: str, seq2: str) -> Tuple[AlignmentResult, str, str]:
        """
        Global alignment using Needleman-Wunsch algorithm
        """
        m, n = len(seq1), len(seq2)
        
        # Initialize scoring matrix
        score_matrix = np.zeros((m + 1, n + 1))
        
        # Initialize first row and column
        for i in range(1, m + 1):
            score_matrix[i][0] = i * self.gap_penalty
        for j in range(1, n + 1):
            score_matrix[0][j] = j * self.gap_penalty
        
        # Fill scoring matrix
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match = score_matrix[i-1][j-1] + (
                    self.match_score if seq1[i-1] == seq2[j-1] else self.mismatch_score
                )
                delete = score_matrix[i-1][j] + self.gap_penalty
                insert = score_matrix[i][j-1] + self.gap_penalty
                
                score_matrix[i][j] = max(match, delete, insert)
        
        # Traceback
        aligned_seq1, aligned_seq2 = self._traceback_nw(
            score_matrix, seq1, seq2
        )
        
        # Calculate alignment statistics
        identity = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-')
        gaps = aligned_seq1.count('-') + aligned_seq2.count('-')
        
        result = AlignmentResult(
            score=score_matrix[m][n],
            alignment_length=len(aligned_seq1),
            identity=identity / len(aligned_seq1) * 100,
            gaps=gaps,
            query_start=0,
            query_end=len(seq1),
            target_start=0,
            target_end=len(seq2)
        )
        
        return result, aligned_seq1, aligned_seq2
    
    def _traceback_sw(self, matrix, seq1, seq2, max_pos):
        """Traceback for Smith-Waterman"""
        aligned_seq1 = ""
        aligned_seq2 = ""
        i, j = max_pos
        
        while i > 0 and j > 0 and matrix[i][j] > 0:
            current_score = matrix[i][j]
            diagonal_score = matrix[i-1][j-1]
            up_score = matrix[i-1][j]
            left_score = matrix[i][j-1]
            
            if current_score == diagonal_score + (
                self.match_score if seq1[i-1] == seq2[j-1] else self.mismatch_score
            ):
                aligned_seq1 = seq1[i-1] + aligned_seq1
                aligned_seq2 = seq2[j-1] + aligned_seq2
                i -= 1
                j -= 1
            elif current_score == up_score + self.gap_penalty:
                aligned_seq1 = seq1[i-1] + aligned_seq1
                aligned_seq2 = "-" + aligned_seq2
                i -= 1
            else:
                aligned_seq1 = "-" + aligned_seq1
                aligned_seq2 = seq2[j-1] + aligned_seq2
                j -= 1
        
        return aligned_seq1, aligned_seq2
    
    def _traceback_nw(self, matrix, seq1, seq2):
        """Traceback for Needleman-Wunsch"""
        aligned_seq1 = ""
        aligned_seq2 = ""
        i, j = len(seq1), len(seq2)
        
        while i > 0 or j > 0:
            if i > 0 and j > 0:
                current_score = matrix[i][j]
                diagonal_score = matrix[i-1][j-1]
                up_score = matrix[i-1][j]
                left_score = matrix[i][j-1]
                
                if current_score == diagonal_score + (
                    self.match_score if seq1[i-1] == seq2[j-1] else self.mismatch_score
                ):
                    aligned_seq1 = seq1[i-1] + aligned_seq1
                    aligned_seq2 = seq2[j-1] + aligned_seq2
                    i -= 1
                    j -= 1
                elif current_score == up_score + self.gap_penalty:
                    aligned_seq1 = seq1[i-1] + aligned_seq1
                    aligned_seq2 = "-" + aligned_seq2
                    i -= 1
                else:
                    aligned_seq1 = "-" + aligned_seq1
                    aligned_seq2 = seq2[j-1] + aligned_seq2
                    j -= 1
            elif i > 0:
                aligned_seq1 = seq1[i-1] + aligned_seq1
                aligned_seq2 = "-" + aligned_seq2
                i -= 1
            else:
                aligned_seq1 = "-" + aligned_seq1
                aligned_seq2 = seq2[j-1] + aligned_seq2
                j -= 1
        
        return aligned_seq1, aligned_seq2
