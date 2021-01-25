import numpy as np

class PairwiseAligner:
    #still need to def init methods here, despite init.py file
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
        matrix = []
        header = []
        with open('../scoring_matrices/'+scoringMatrixFile,'r') as f:
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
        seq = ''
        with open(fafsa, 'r') as f:
            for line in f:
                if '>' not in line: #this might not catch every exception
                    seq += line.upper().rstrip() #strip end of line char
        f.close()
        return seq

    def getScoringMatrixCharIndex(self, char): #could do this with a hash table (dict)
        for i in range(len(self.scoringMatrixHeader)):
            if self.scoringMatrixHeader[i] == char:
                return i
        raise Exception("Character is not in the Scoring Matrix Header :( ")

    def setGapOpening(self, val):
        self.gapOpening = val * -1

    def setGapExtension(self, val):
        self.gapExtension = val * -1

    def makeMxNMatrix(self, seq1, seq2):
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

    #0 corresp to 0 in SW, 1 is diagonal, 2 is left, 3 is up
    def calcTracebackCell(self, maxVal, i, j, row, col, isNW):
        if not isNW and maxVal == 0:
            self.setTMatrixCell(0, row+1, col+1)
        elif maxVal == self.mMatrix[row][col] + self.scoringMatrix[j][i]: #diagonal
            self.setTMatrixCell(1, row+1, col+1) #plus one here because of offset. Matches mMatrix
        elif maxVal == self.IAMatrix[row+1][col+1]: #left
            self.setTMatrixCell(2, row+1, col+1)
        elif maxVal == self.IBMatrix[row+1][col+1]: #up
            self.setTMatrixCell(3, row+1, col+1)

    def pairwiseAlign(self, isNW = True):
        for col in range(len(self.seq1)): #col
            for row in range(len(self.seq2)): #row
                char1 = self.seq1[col]
                char2 = self.seq2[row]
                i = self.getScoringMatrixCharIndex(char2) #row
                j = self.getScoringMatrixCharIndex(char1) #col

                #since we're setting r+1 and c+1 in mMatrix, we should take max of r and c (instead of (r,c) and (r-1,c-1) like in the notes)

                #update the IA and IB matrices first
                self.setIAMatrixCell(max(self.mMatrix[row][col+1] + self.gapOpening, self.IAMatrix[row][col+1] + self.gapExtension), row+1, col+1)
                self.setIBMatrixCell(max(self.mMatrix[row+1][col] + self.gapOpening, self.IBMatrix[row+1][col] + self.gapExtension), row+1, col+1)

                #AFFINE just depends on whether prev align IN GIVEN SEQ had gap or not
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

                #note 1.18.21: I need to set row 0 and col 0 of the mMatrix before I can calc the other values in the matrix ([0][0] = 0), then
                #   I can start using neighboring values to calculate Sij one row (or col) at a time
                #I can have each subclass call the setMMatrixCell method with row=0 and col=0, and have them give their respective starting
                #   values that way
                #1.21.21: don't I already set all cells to 0 when I make the matrix, though?
                #Hopefully the np.int32 object is not iterable error goes away once I give the max f(x) multiple arguments - it does


    #traceback method must be called after align method
    #Important to note that the matrix indeces are all +1 in relation to corresponding seq indeces
    def pairwiseTraceback(self, row, col):
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
        PairwiseAligner.pairwiseAlign(self, isNW = False)
        self.alignmentScore = self.swHighestVal

    def traceback(self):
        self.alignedSequences = self.pairwiseTraceback(self.swHighestValIndex[0], self.swHighestValIndex[1])


class NeedlemanWunsch(PairwiseAligner):
    def __init__(self, filename, fafsa1, fafsa2, gapOpening, gapExtension):
        PairwiseAligner.__init__(self, filename, fafsa1, fafsa2, gapOpening, gapExtension)
        self.fillFirstRowAndColumn()

    def fillFirstRowAndColumn(self):
        for j in range(len(self.mMatrix)-1):
            self.setMMatrixCell(self.mMatrix[j][0] + self.gapOpening, j+1, 0)
        for i in range(len(self.mMatrix[0])-1):
            self.setMMatrixCell(self.mMatrix[0][i] + self.gapOpening, 0, i+1)

    def align(self):
        PairwiseAligner.pairwiseAlign(self, isNW = True)
        self.alignmentScore = self.mMatrix[-1][-1]

    def traceback(self):
        self.alignedSequences = self.pairwiseTraceback(len(self.mMatrix)-1, len(self.mMatrix[0])-1)



fa1 = '../sequences/prot-0004.fa'
fa2 = '../sequences/prot-0008.fa'
scoringMatrixFile = 'BLOSUM50.mat'
gapOpening = 11
gapExtension = 3

"""
nwa = NeedlemanWunsch(scoringMatrixFile, fa1, fa2, gapOpening, gapExtension)
nwa.align()
print(nwa.mMatrix)
print(nwa.mMatrix.shape)
print(nwa.tMatrix)
print(nwa.tMatrix.shape)
print(len(nwa.seq1))
print(len(nwa.seq2))
nwa.traceback()
print(nwa.alignedSequences)


swa = SmithWaterman(scoringMatrixFile, fa1, fa2, gapOpening, gapExtension)
swa.align()
print(swa.mMatrix)
print(swa.mMatrix.shape)
print(swa.tMatrix)
print(swa.tMatrix.shape)
print(len(swa.seq1))
print(len(swa.seq2))
swa.traceback()
print(swa.alignedSequences)
"""



"""
######
# Unit Tests
######
fa1 = '../sequences/prot-0004.fa'
fa2 = '../sequences/prot-0008.fa'
scoringMatrixFile = 'BLOSUM50.mat'
gapOpening = 11
gapExtension = 3

swaTest = SmithWaterman(scoringMatrixFile, fa1, fa2, gapOpening, gapExtension)
nwaTest = NeedlemanWunsch(scoringMatrixFile, fa1, fa2, gapOpening, gapExtension)


def ioFastaTest(swaTest):
    assert swaTest.readFafsaSeq('../sequences/prot-0004.fa') == 'SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM'
    print('"I/O of Protein FAFSA Sequence" Test Passed')

ioFastaTest(swaTest)

def ioScoringMatrixTest(swaTest):
    testScoringMatrix, testHeader = swaTest.readScoringMatrix('BLOSUM50.mat')
    assert testScoringMatrix[0][0] == 5 #how do I avoid truth value of array ValueError?
    assert testHeader[0] == 'A'
    print('"I/O of Scoring Matrix" Test Passed')

ioScoringMatrixTest(swaTest)

def identicalSequencesTest():
    swaTest2 = SmithWaterman(scoringMatrixFile, fa1, fa1, gapOpening, gapExtension)
    nwaTest2 = NeedlemanWunsch(scoringMatrixFile, fa1, fa1, gapOpening, gapExtension)
    swaTest2.align()
    swaTest2.traceback()
    nwaTest2.align()
    nwaTest2.traceback()
    assert swaTest2.alignedSequences[0] == swaTest2.alignedSequences[1]
    assert nwaTest2.alignedSequences[0] == nwaTest2.alignedSequences[1]
    print('"Identical Sequences" Test Passed')

identicalSequencesTest()

def alignmentScoreTest():
    pass

"""


"""

#####
#Part 2
#####

import matplotlib.pyplot as plt

#1)
scoringMatrixFile = 'BLOSUM50.mat'
gapOpening = 11
gapExtension = 3

allAlignmentScores = []

with open('../scoring_matrices/Pospairs.txt','r') as f:
    for line in f:
        fa1, fa2 = line.split()
        fa1 = '../' + fa1
        fa2 = '../' + fa2
        swa = SmithWaterman(scoringMatrixFile, fa1, fa2, gapOpening, gapExtension)
        swa.align()
        allAlignmentScores.append(swa.alignmentScore)
f.close()


with open('../scoring_matrices/Negpairs.txt','r') as f:
    for line in f:
        fa1, fa2 = line.split()
        fa1 = '../' + fa1
        fa2 = '../' + fa2
        swa = SmithWaterman(scoringMatrixFile, fa1, fa2, gapOpening, gapExtension)
        swa.align()
        allAlignmentScores.append(swa.alignmentScore)
f.close()

print(allAlignmentScores)

plt.scatter(list(range(0,len(allAlignmentScores))), allAlignmentScores)
plt.show()

"""




