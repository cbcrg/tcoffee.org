2.1.msa2ungappedmsa.sol.py msa1.fasta 0 >msa1.ungapped.fasta
2.2.fastafile2lom.sol.py  msa1.fasta >msa1.lom
2.2.fastafile2lom.sol.py  msa1.ungapped.fasta >msa1.ungapped.lom
2.5.compare_matrix.sol.py msa1.lom msa1.ungapped.lom 
