import numpy as np

class PairwiseAligner:
    #still need to def init methods here, despite init.py file
    def __init__(self, filename):
        self.filename = filename
        self.scoreMatrix = self.readScoreMatrix()

    def readScoreMatrix(self):
        matrix = []
        with open('../scoring_matrices/'+self.filename,'r') as f:
            for line in f:
                splitLine = line.split()
                if splitLine[0].isalpha() or splitLine[0].isnumeric() or splitLine[0][0] == '-':
                    matrix.append(line.split())
        f.close()
        npMatrix = np.array(matrix)
        print(npMatrix)
        print(npMatrix.shape) #(25,24) with extra row with AA labels
        #print(npmatrix[1][0]) #is 5, checks out
        return npMatrix

    def align(self, seq1, seq2):
        with open(seq1, 'r'):


    pass

class SmithWaterman(PairwiseAligner):
    pass

class NeedlemanWunsch(PairwiseAligner):
    pass


pwa = PairwiseAligner('BLOSUM50.mat')
#pwa2 = PairwiseAligner('PAM250.mat')

#pwa.readScoreMatrix()
#pwa2.readScoreMatrix()