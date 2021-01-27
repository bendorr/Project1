# Project 1 - Sequence Alignment
## Due 01/27/2021

![BuildStatus](https://github.com/bendorr2415/Project1/workflows/HW1/badge.svg?event=push)

In this assignment, you will implement two classical alignment algorithms and then evaluate each algorithm’s performance with a range of parameters. There are two parts to this assignment and Part 2 requires completion of Part 1. We recommend reading through both Part 1 and Part 2 before beginning this assignment. 

* Part 1 - API and implementation
* Part 2 - Evaluating alignments

### main
Runs all code in align/\_\_main\_\_.py, useful for part 2
```
python -m align
```

### testing
Testing is as simple as running
```
python -m pytest test/*
```
from the root directory of this project.

### API

```
# Class PairwiseAligner:

  - Contains methods and attributes that are shared between both Smith-Waterman and Needleman-Wunsch alignment algorithms.
    - Attributes:
      - gap opening and gap extension penalties
      - scoring matrix and the matrix header (order of amino acids in the rows and columns of the matrix)
      - seq1 and seq2: the two sequences which are being aligned
      - 4 different M x N matrices: the alignment score matrix, the traceback matrix, the gaps in sequence A matrix (IA), and the gaps in sequence B matrix (IB)
      - the aligned sequences as a tuple of two strings
      - the alignment score
      - swHighestVal and swHighestValIndex: the value and index of the highest value in the alignment score matrix, which is only set during S-W alignments
      
    - Methods:
      - readScoringMatrix:
      Reads the scoring matrix from the given file into two numpy arrays; one for the matrix header (1D) and
        one for the scores (2D)

            Param:
                scoringMatrixFile = the name of the scoring matrix file in the scoring_matrices folder

            Returns:
                a tuple containing two numpy arrays: one containing the 2D scoring matrix itself (containing the scores),
                    the other containing the 1D header of the scoring matrix (with the amino acids in order)
                    
      - readFafsaSeq:
      Reads the given fafsa (.fa) file and returns the sequence contained in it as a string

            Param:
                fafsa = the name of the .fa file containing a sequence

            Returns:
                the sequence contained in the .fa file as a string
      
      - getScoringMatrixCharIndex:
      Returns the index of the given amino acid one-character-code in the scoring matrix header

            Param:
                char = one character code for an amino acid, contained in the scoring matrix header

            Returns:
                the index of the given character in the scoringMatrixHeader array
                
      - setter functions: setGapOpening, setGapExtension, setM- T- IA- and IB- MatrixCell
      
      - makeMxNMatrix: returns a numpy matrix with M+1 ((length of sequence 1) + 1) columns and N+1 ((length of sequence 2) + 1) rows
      Creates and returns a 2D numpy array whose dimensions are the lengths of the two input
        sequences plus one

            Params:
                seq1, seq2 = the sequences which are being aligned by the PairwiseAligner object

            Returns:
                a 2D numpy matrix with len(seq2)+1 rows and len(seq1)+1 columns
                
      - pairwiseAlign:
      Iterates through each cell in the four M x N matrices (M - scoring matrix, T - traceback matrix, IA - gaps in sequence 1
        matrix, and IB - gaps in sequence 2 matrix) and fills them with their appropriate values.

            Param:
                isNW = a boolean that is True if the current alignment algorithm is Needleman-Wunsch, False if the current
                    algorithm is Smith-Waterman

            Returns:
                None
                
      - calcTracebackCell:
      Determines which of the possible values for the alignment's M Matrix was set at the current cell, and depending on which
        value was inserted into the M Matrix, inserts a 0, 1, 2, or 3 in the T Matrix's corresponding cell.

            Params:
                maxVal = the value that was inserted into the M Matrix at [row+1][col+1]
                i = the column in the scoring matrix that corresponds to the alignment score of the
                    current pair of amino acids
                j = the row in the scoring matrix that corresponds to the alignment score of the
                    current pair of amino acids
                row = the previous row to the T Matrix cell that is currently being considered. The cell at row (row+1)
                    will be set by this function
                col = the previous column to the T Matrix cell that is currently being considered. The cell at column (column+1)
                    will be set by this function
                isNW = a boolean that is True if the current alignment algorithm is Needleman-Wunsch, False if the current
                    algorithm is Smith-Waterman

            Returns:
                None
```
