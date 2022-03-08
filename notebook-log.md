Files retrieved, put into .fa form using SnapGene program
files made into 1 .fasta file using following command
cat *fa > alldata.fasta

Files aligned using ClustalW
$/c/'Program Files (x86)'/clustalW2/clustalw2/ -ALIGN -INFILE=alldata.fasta -OUTFILE=ALIGNEDalldata.fasta -OUTPUT=FASTA