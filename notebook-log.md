Files retrieved, put into .fa form using SnapGene program
files made into 1 .fasta file using following command
cat *fa > alldata.fasta

Files aligned using ClustalW
$/c/'Program Files (x86)'/clustalW2/clustalw2/ -ALIGN -INFILE=alldata.fasta -OUTFILE=ALIGNEDalldata.fasta -OUTPUT=FASTA

RGui using Distance based and Parsimony based methods
Distance Based
	install.packages("adegenet", dep=TRUE)
	install.packages("phangorn", dep=TRUE)
	library(ape)
	library(adegenet)
	library(phangorn)
	dna <- fasta2DNAbin(file="C:/Users/logan/Desktop/Phylo/Actual Data/ALIGNEDalldata.fasta")
	D <- dist.dna(dna, model="TN93")
	tre <- fastme.bal(D)
	tre <- ladderize(tre)
	plot(tre, cex=.6)
	title("Relationship between Gardnerella Strains and Samples")
Parsimony Based
	dna2 <- as.phyDat(dna)
	class(dna2)
	dna2
	tre.ini <- fastme.bal(dist.dna(dna,model="raw"))
	tre.ini
	parsimony(tre.ini, dna2)
	tre.pars <- optim.parsimony(tre.ini, dna2)
	tre.pars
	plot(tre.pars, cex=6)
	title("Maximum-parsimony tree")

	