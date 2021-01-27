import numpy as np


class PairwiseAligner:
    
    def __init__(self, scoringMatrixFile, fafsa1, fafsa2, gapOpening, gapExtension):
        self.setGapOpening(gapOpening)
        self.setGapExtension(gapExtension)
        self.scoringMatrix, self.scoringMatrixHeader = self.readScoringMatrix(scoringMatrixFile)
        self.seq1 = self.readFafsaSeq(fafsa1) #columns
        self.seq2 = self.readFafsaSeq(fafsa2) #rows
        self.mMatrix = self.makeMxNMatrix(self.seq1, self.seq2)
        self.tMatrix = self.makeMxNMatrix(self.seq1, self.seq2)

        self.swHighestVal = 0
        self.swHighestValIndex = (0, 0)

        self.IAMatrix = self.makeMxNMatrix(self.seq1, self.seq2)
        self.IBMatrix = self.makeMxNMatrix(self.seq1, self.seq2)

        self.alignedSequences = ('','')
        self.alignmentScore = -999


    def readScoringMatrix(self, scoringMatrixFile):
        """
        Reads the scoring matrix from the given file into two numpy arrays; one for the matrix header (1D) and
        one for the scores (2D)

            Param:
                scoringMatrixFile = the name of the scoring matrix file in the scoring_matrices folder

            Returns:
                a tuple containing two numpy arrays: one containing the 2D scoring matrix itself (containing the scores),
                    the other containing the 1D header of the scoring matrix (with the amino acids in order)
        """
        matrix = []
        header = []
        with open('scoring_matrices/'+scoringMatrixFile,'r') as f:
            for line in f:
                splitLine = line.split()
                if splitLine[0].isalpha():
                    header = line.split()
                elif splitLine[0].isnumeric() or splitLine[0][0] == '-':
                    matrix.append(line.split())
        f.close()
        npMatrix = np.array(matrix, np.int32)
        npHeader = np.array(header)
        return npMatrix, npHeader


    def readFafsaSeq(self, fafsa):
        """
        Reads the given fafsa (.fa) file and returns the sequence contained in it as a string

            Param:
                fafsa = the name of the .fa file containing a sequence

            Returns:
                the sequence contained in the .fa file as a string
        """
        seq = ''
        with open(fafsa, 'r') as f:
            for line in f:
                if '>' not in line: #this might not catch every exception
                    seq += line.upper().rstrip() #strip end of line char
        f.close()
        return seq


    def getScoringMatrixCharIndex(self, char): #could do this with a hash table (dict)
        """
        Returns the index of the given amino acid one-character-code in the scoring matrix header

            Param:
                char = one character code for an amino acid, contained in the scoring matrix header

            Returns:
                the index of the given character in the scoringMatrixHeader array
        """
        for i in range(len(self.scoringMatrixHeader)):
            if self.scoringMatrixHeader[i] == char:
                return i
        raise Exception("Character is not in the Scoring Matrix Header :( ")


    def setGapOpening(self, val):
        self.gapOpening = val * -1


    def setGapExtension(self, val):
        self.gapExtension = val * -1


    def makeMxNMatrix(self, seq1, seq2):
        """
        Creates and returns a 2D numpy array whose dimensions are the lengths of the two input
        sequences plus one

            Params:
                seq1, seq2 = the sequences which are being aligned by the PairwiseAligner object

            Returns:
                a 2D numpy matrix with len(seq2)+1 rows and len(seq1)+1 columns
        """
        matrix = []
        #maybe use numpy.zeros instead
        for r in range(len(self.seq2)+1):
            row = []
            for c in range(len(self.seq1)+1):
                row.append(0)
            matrix.append(row)
        return np.array(matrix)


    def setMMatrixCell(self, value, row, col):
        self.mMatrix[row][col] = value


    def setTMatrixCell(self, value, row, col):
        self.tMatrix[row][col] = value


    def setIAMatrixCell(self, value, row, col):
        self.IAMatrix[row][col] = value


    def setIBMatrixCell(self, value, row, col):
        self.IBMatrix[row][col] = value


    def calcTracebackCell(self, maxVal, i, j, row, col, isNW):
        """
        Determines which of the possible values for the alignment's M Matrix was set at the current cell, and depending on which
        value was inserted into the M Matrix, inserts a 0, 1, 2, or 3 in the T Matrix's corresponding cell.  If the alignment algorithm is
        Smith-Waterman and the value inserted into the M Matrix was 0, then a 0 is inserted into the T Matrix.  If the M Matrix's value came
        from adding the alignment score of the current amino acid pair to the value in the diagonal (upper left) cell, a 1 is inserted into the
        T Matrix. If the M Matrix's value came from adding a gap penalty to the cell to the left, then a 2 is inserted.  And finally, if the 
        M Matrix's value came from adding a gap penalty from the cell above, a 3 is inserted into the T Matrix.

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
        """
        if not isNW and maxVal == 0: #a zero in a SmithWaterman scoring matrix
            self.setTMatrixCell(0, row+1, col+1)
        elif maxVal == self.mMatrix[row][col] + self.scoringMatrix[j][i]: #diagonal
            self.setTMatrixCell(1, row+1, col+1) #plus one here because of offset. Matches mMatrix
        elif maxVal == self.IAMatrix[row+1][col+1]: #left
            self.setTMatrixCell(2, row+1, col+1)
        elif maxVal == self.IBMatrix[row+1][col+1]: #up
            self.setTMatrixCell(3, row+1, col+1)


    def pairwiseAlign(self, isNW = True):
        """
        Iterates through each cell in the four M x N matrices (M - scoring matrix, T - traceback matrix, IA - gaps in sequence 1
        matrix, and IB - gaps in sequence 2 matrix) and fills them with their appropriate values.
        IA Matrix cells are filled with the maximum between 1) the alignment score (M Matrix value) in the cell above  + the gap
        opening penalty and 2) the IA matrix's value in the cell above + the gap extension penalty). IB Matrix cells are filled in
        the same way as IA matrix cells but looking at the cells to the left instead of the cells above.
        M Matrix cells are set as the maximum value between 1) the upper-left (diagonal) M Matrix value + the scoring matrix's
        alignment score for the current amino acid pair, 2) the IA Matrix's corresponding cell's value, 3) the IB Matrix's
        corresponding cell's value, and, for Smith-Waterman alignments only, 4) 0.
        T Matrix cells are set depending on what value was set in the corresponding M Matrix cell, using the function calcTracebackCell


            Param:
                isNW = a boolean that is True if the current alignment algorithm is Needleman-Wunsch, False if the current
                    algorithm is Smith-Waterman

            Returns:
                None
        """
        for col in range(len(self.seq1)): #col
            for row in range(len(self.seq2)): #row
                char1 = self.seq1[col]
                char2 = self.seq2[row]
                i = self.getScoringMatrixCharIndex(char2) #row
                j = self.getScoringMatrixCharIndex(char1) #col

                #since we're setting r+1 and c+1 in mMatrix, we should look at r and c to get values from previous rows and columns
                self.setIAMatrixCell(max(self.mMatrix[row][col+1] + self.gapOpening, self.IAMatrix[row][col+1] + self.gapExtension), row+1, col+1)
                self.setIBMatrixCell(max(self.mMatrix[row+1][col] + self.gapOpening, self.IBMatrix[row+1][col] + self.gapExtension), row+1, col+1)

                if isNW:
                    maxVal = max(self.mMatrix[row][col] + self.scoringMatrix[j][i], self.IAMatrix[row+1][col+1], self.IBMatrix[row+1][col+1])
                    self.calcTracebackCell(maxVal, i, j, row, col, isNW)
                else:
                    maxVal = max(self.mMatrix[row][col] + self.scoringMatrix[j][i], self.IAMatrix[row+1][col+1], self.IBMatrix[row+1][col+1], 0)
                    if maxVal > self.swHighestVal:
                        self.swHighestVal = maxVal
                        self.swHighestValIndex = (row+1, col+1)
                    self.calcTracebackCell(maxVal, i, j, row, col, isNW = False)
                self.setMMatrixCell(maxVal, row+1, col+1) #plus ones are for the offset with the gaps are row and column zero


    def pairwiseTraceback(self, row, col):
        """
        Performs traceback algorithm over the tMatrix. Returns a tuple of the two aligned sequences.
        This method must be called after the align method, since align fills the tMatrix.

            Params:
                row = row from which to begin the traceback
                col = column from which to begin the traceback

            Returns:
                A tuple containing each of the two aligned sequences as strings

        """

        # Note that the matrix indeces are all +1 in relation to corresponding seq indeces, due to the offset caused 
        #   by inserting the first row and column in the matrix
        s1, s2 = '', ''
        while(row > -1 or col > -1):
            if(self.tMatrix[row][col] == 0):
                break
            elif(self.tMatrix[row][col] == 1):
                s1 = self.seq1[col-1] + s1
                s2 = self.seq2[row-1] + s2
                row -= 1
                col -= 1
            elif(self.tMatrix[row][col] == 2): #left
                s1 = self.seq1[col-1] + s1
                s2 = '-' + s2
                col -= 1
            elif(self.tMatrix[row][col] == 3): #up
                s1 = '-' + s1
                s2 = self.seq2[row-1] + s2
                row -= 1
        return (s1, s2)


class SmithWaterman(PairwiseAligner):

    def __init__(self, filename, fafsa1, fafsa2, gapOpening, gapExtension):
        PairwiseAligner.__init__(self, filename, fafsa1, fafsa2, gapOpening, gapExtension)


    def align(self):
        """
        Runs the pairwiseAlign method from the PairwiseAligner class with the parameter isNW = False.
        Sets its alignmentScore attribute to equal the "Smith-Waterman Highest Value in the M Matrix" value
        that is set during the pairwiseAlign method.

            Params:
                None

            Returns:
                None
        """
        PairwiseAligner.pairwiseAlign(self, isNW = False)
        self.alignmentScore = self.swHighestVal


    def traceback(self):
        """
        Runs the pairwiseTraceback method from the PairwiseAligner class with parameters that indicate
        the location of the cell containing the highest alignment score in the Smith-Waterman alignment score
        matrix (which is where the traceback for the Smith-Waterman alignment algorithm begins), and saves
        the tuple of aligned sequences that is returned in its alignedSequences attribute.

            Params:
                None

            Returns:
                None
        """
        self.alignedSequences = self.pairwiseTraceback(self.swHighestValIndex[0], self.swHighestValIndex[1])



class NeedlemanWunsch(PairwiseAligner):

    def __init__(self, filename, fafsa1, fafsa2, gapOpening, gapExtension):
        PairwiseAligner.__init__(self, filename, fafsa1, fafsa2, gapOpening, gapExtension)
        self.fillFirstRowAndColumn()


    def fillFirstRowAndColumn(self):
        """
        The M Matrix is initialized to be full of zeros, but the Needleman-Wunsch algorithm requires
        the first row and first column of its scoring matrix (M Matrix) to be filled with compounding
        gap opening penalties (for example, if the gap opening penalty = 5, the values in the first row
        should be: -5, -10, -15, -20, etc.). This function fills the first row and first column of the
        M Matrix with compounding gap opening penalties.

            Params:
                None

            Returns:
                None
        """
        for j in range(len(self.mMatrix)-1):
            self.setMMatrixCell(self.mMatrix[j][0] + self.gapOpening, j+1, 0)
        for i in range(len(self.mMatrix[0])-1):
            self.setMMatrixCell(self.mMatrix[0][i] + self.gapOpening, 0, i+1)


    def align(self):
        """
        Runs the pairwiseAlign method from the PairwiseAligner class with the parameter isNW = True.
        Sets its alignmentScore attribute to equal the value in the bottom right corner of the
        alignment score matrix.

            Params:
                None

            Returns:
                None
        """
        PairwiseAligner.pairwiseAlign(self, isNW = True)
        self.alignmentScore = self.mMatrix[-1][-1]


    def traceback(self):
        """
        Runs the pairwiseTraceback method from the PairwiseAligner class with parameters that indicate
        the location of the cell in the bottow right corner of the Needleman-Wunsch alignment score matrix,
        and saves the tuple of aligned sequences that is returned in its alignedSequences attribute.

            Params:
                None

            Returns:
                None
        """
        self.alignedSequences = self.pairwiseTraceback(len(self.mMatrix)-1, len(self.mMatrix[0])-1)




