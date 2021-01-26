import pytest
from align import algs

fa1 = 'sequences/prot-0004.fa'
fa2 = 'sequences/prot-0008.fa'
scoringMatrixFile = 'scoring_matrices/BLOSUM50.mat'
gapOpening = 11
gapExtension = 3

swaTest = algs.SmithWaterman(scoringMatrixFile, fa1, fa2, gapOpening, gapExtension)
nwaTest = algs.NeedlemanWunsch(scoringMatrixFile, fa1, fa2, gapOpening, gapExtension)

@pytest.fixture
def some_relevant_data():
	return np.ones(10)

def test_fasta_io(swaTest):
	assert swaTest.readFafsaSeq('sequences/prot-0004.fa') == 'SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM'
	print('"I/O of Protein FAFSA Sequence" Test Passed')
	
def test_scoring_matrix_io(swaTest):
	testScoringMatrix, testHeader = swaTest.readScoringMatrix('BLOSUM50.mat')
	assert testScoringMatrix[0][0] == 5 #how do I avoid truth value of array ValueError?
	assert testHeader[0] == 'A'
	print('"I/O of Scoring Matrix" Test Passed')
	
def test_identical():
	swaTest2 = algs.SmithWaterman(scoringMatrixFile, fa1, fa1, gapOpening, gapExtension)
	nwaTest2 = algs.NeedlemanWunsch(scoringMatrixFile, fa1, fa1, gapOpening, gapExtension)
	swaTest2.align()
	swaTest2.traceback()
	nwaTest2.align()
	nwaTest2.traceback()
	assert swaTest2.alignedSequences[0] == swaTest2.alignedSequences[1]
	assert nwaTest2.alignedSequences[0] == nwaTest2.alignedSequences[1]
	print('"Identical Sequences" Test Passed')
	
def test_alignment_score():
	swaTest.seq1 = 'ARN'
	swaTest.seq2 = 'ARN'
	swaTest.align()
	assert swaTest.alignmentScore == 19
	print('"Alignment Score" Test Passed')
