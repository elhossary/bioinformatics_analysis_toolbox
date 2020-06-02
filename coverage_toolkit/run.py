from wiggle_matrix import WiggleMatrix
import glob
import os
fasta_path = os.path.abspath("../test_data/cholerae/v_cholerae_O1_N16961.fa")
file_pathes = [os.path.abspath(i) for i in glob.glob("../test_data/cholerae/TERM-SEQ_LAST-BASE_*")]
x = WiggleMatrix(fasta_path, file_pathes).get_matrix_by_orientation("+")
print(x)
#a.normalize(a.get_percentile(95))

exit()