import numpy as np

class PairwiseAligner:
    #still need to def init methods here, despite init.py file
    def __init__(self, filename, fafsa1, fafsa2):
        self.filename = filename
        self.scoringMatrix, self.scoringMatrixHeader = self.readScoringMatrix()
        self.seq1 = self.readFafsaSeq(fafsa1) #columns
        self.seq2 = self.readFafsaSeq(fafsa2) #rows
        self.mMatrix = self.makeMxNMatrix(self.seq1, self.seq2)
        self.tMatrix = self.makeMxNMatrix(self.seq1, self.seq2)
        self.gapPenalty = self.readGapPenalty()

        self.swHighestVal = 0
        self.swHighestValIndex = (0, 0)

    def readScoringMatrix(self):
        matrix = []
        header = []
        with open('../scoring_matrices/'+self.filename,'r') as f:
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

    def readGapPenalty(self):
        i = self.getScoringMatrixCharIndex('*')
        return self.scoringMatrix[i][0]

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

    def pairwiseAlign(self, isNW = True):
        for col in range(len(self.seq1)): #col
            for row in range(len(self.seq2)): #row
                char1 = self.seq1[col]
                #print(char1)
                char2 = self.seq2[row]
                #print(char2)
                #print('pair')
                i = self.getScoringMatrixCharIndex(char2) #row
                j = self.getScoringMatrixCharIndex(char1) #col
                #print(self.scoringMatrix[j][i])
                #since we're setting r+1 and c+1 in mMatrix, we should take max of r and c (instead of (r,c) and (r-1,c-1) like in the notes)
                if isNW:
                    maxVal = max(self.mMatrix[row][col] + self.scoringMatrix[j][i], self.mMatrix[row][col+1] + self.gapPenalty, \
                                self.mMatrix[row+1][col] + self.gapPenalty)
                else:
                    maxVal = max(self.mMatrix[row][col] + self.scoringMatrix[j][i], self.mMatrix[row][col+1] + self.gapPenalty, \
                                self.mMatrix[row+1][col] + self.gapPenalty, 0)
                    if maxVal > self.swHighestVal:
                        self.swHighestVal = maxVal
                        self.swHighestValIndex = (row+1, col+1)
                self.setMMatrixCell(maxVal, row+1, col+1) #plus ones are for the offset with the gaps are row and column zero


                #self.setTMatrixCell() based on what maxVal is equal to


                #note 1.18.21: I need to set row 0 and col 0 of the mMatrix before I can calc the other values in the matrix ([0][0] = 0), then
                #   I can start using neighboring values to calculate Sij one row (or col) at a time
                #I can have each subclass call the setMMatrixCell method with row=0 and col=0, and have them give their respective starting
                #   values that way
                #1.21.21: don't I already set all cells to 0 when I make the matrix, though?
                #Hopefully the np.int32 object is not iterable error goes away once I give the max f(x) multiple arguments - it does


class SmithWaterman(PairwiseAligner):
    def align(self):
        PairwiseAligner.pairwiseAlign(self, isNW = False)

class NeedlemanWunsch(PairwiseAligner):
    def __init__(self, filename, fafsa1, fafsa2):
        PairwiseAligner.__init__(self, filename, fafsa1, fafsa2)
        self.fillFirstRowAndColumn()

    def fillFirstRowAndColumn(self):
        for j in range(len(self.mMatrix)-1):
            self.setMMatrixCell(self.mMatrix[j][0] + self.gapPenalty, j+1, 0)
        for i in range(len(self.mMatrix[0])-1):
            self.setMMatrixCell(self.mMatrix[0][i] + self.gapPenalty, 0, i+1)

    def align(self):
        PairwiseAligner.pairwiseAlign(self, isNW = True)



fa1 = '../sequences/prot-0004.fa'
fa2 = '../sequences/prot-0008.fa'

"""
#Testing PWA class
pwa = PairwiseAligner('BLOSUM50.mat', fa1, fa2)
print(len(pwa.scoringMatrixHeader))
print(pwa.scoringMatrixHeader)
# print(pwa.mMatrix.shape)
# print(pwa.mMatrix)
print('gap penalty is ' + pwa.gapPenalty)
pwa.align()
#print(pwa.scoringMatrix)
print(pwa.mMatrix)
print(pwa.mMatrix.shape)
print(len(pwa.seq1))
print(len(pwa.seq2))
"""

nwa = NeedlemanWunsch('BLOSUM50.mat', fa1, fa2)
nwa.align()
print(nwa.scoringMatrix)
print(nwa.mMatrix)
print(nwa.mMatrix.shape)
print(len(nwa.seq1))
print(len(nwa.seq2))

swa = SmithWaterman('BLOSUM50.mat', fa1, fa2)
swa.align()
print(swa.mMatrix)
print(swa.mMatrix.shape)





#pwa.align('../sequences/prot-0004.fa', '../sequences/prot-0008.fa') #works

#pwa2 = PairwiseAligner('PAM250.mat')

#pwa.readScoreMatrix()
#pwa2.readScoreMatrix()