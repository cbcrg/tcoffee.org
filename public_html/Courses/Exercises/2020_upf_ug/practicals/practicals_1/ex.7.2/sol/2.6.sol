python3.7 2.1.msa2ungappedmsa.sol.py msa1.fasta 0 >msa1.ungapped.fasta
python3.7 2.2.fastafile2lom.sol.py  msa1.fasta >msa1.lom
python3.7 2.2.fastafile2lom.sol.py  msa1.ungapped.fasta >msa1.ungapped.lom
python3.7 2.5.compare_matrix.sol.py msa1.lom msa1.ungapped.lom 
