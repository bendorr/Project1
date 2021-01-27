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
Class PairwiseAligner:
  - Contains methods and attributes that are shared between both Smith-Waterman and Needleman-Wunsch alignment algorithms.
    - Attributes:
      - gap opening and gap extension penalties
      - scoring matrix and the matrix header (order of amino acids in the rows and columns of the matrix)
      - the two sequences which are being aligned
      - 4 different M x N matrices: the alignment score matrix, the traceback matrix, the gaps in sequence A matrix (IA), and the gaps in sequence B matrix (IB)
      - the aligned sequences
      - the alignment score
      - (this class also contains the value and index of the highest value in the alignment score matrix, which is only set during S-W alignments)
    - Methods:
      - readScoringMatrix: read and store (as numpy arrays) the scoring matrix and header
      - readFafsaSeq: read and store the sequence in the given fafsa file
      - getScoringMatrixCharIndex: return the index of the given amino acid in the scoring matrix header
      - setter functions: setGapOpening, setGapExtension, setM- T- IA- and IB- MatrixCell
      - makeMxNMatrix: returns a numpy matrix with M+1 ((length of sequence 1) + 1) columns and N+1 ((length of sequence 2) + 1) rows
      - 
```
