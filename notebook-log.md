Files retrieved, put into .fa form using SnapGene program
files made into 1 .fasta file using following command
cat *fa > alldata.fasta

Files aligned using ClustalW
Slow process (option 1) chosen
$/c/'Program File (x86)'/clustalW2/clustalw2/ -ALIGN -INFILE=Gardnerella.fasta -OUTFILE=ALIGNEDGardnerella.fasta -OUTPUT=FASTA

RGui using Distance based and Parsimony based methods
Distance Based
	install.packages("adegenet", dep=TRUE)
	install.packages("phangorn", dep=TRUE)
	library(ape)
	library(adegenet)
	library(phangorn)
	dna <- fasta2DNAbin(file="C:/Users/logan/Desktop/Phylo/Actual/ALIGNEDGardnerella.fasta")
	D <- dist.dna(dna, model="TN93")
	tre <- fastme.bal(D)
	tre <- ladderize(tre)
	plot(tre, cex=.6)
	title("Distance-Based Gardnerella Strains and Samples")
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

IQ-Tree 
iqtree -s C:\Users\logan\Desktop\Phylo\Actual\ALIGNEDGardnerella.fasta -bb 1000 -nt AUTO
	Part of output is in a splits.nex file due to the inconclusive nature of the data. This is due to the PCR process done in lab which severs the ends of the DNA.

Mesquite File is used to get a .nex file that will work with MrBayes as IQ-Tree does not use a file that works well. This process involves simplifying the data a bit so MrBayes is able to read it. This was done by taking the ALIGNEDGardnerella.fasta file from ClustalW and exporting it as a .nex file designed to function with MrBayes. File saved as
	ALIGNEDGardnerellaBAYES.nex

begin mrbayes;
	set autoclose=yes nowarn=yes
	lset nst=6 rates=invgamma;
	unlink statefreq=(all) revmat=(all)shape=(all)pinvar=(all);
	preset applyto=(all) ratepr=variable;
	mcmcp ngen=10000 relburnin=yes burninfrao=0.25 printfreq=1000 samplefreq=1000 nchains=4 savebrlens=yes
	mcmc;
	sumt;
end;

Execute ALIGNEDGardnerellaBAYES.nex

ASTRAL

java -jar astral.5.7.8.jar -i c:\Users\logan\Desktop\Phylo\Actual\ALIGNEDGardnerella.fasta.treefile -o c:\Users\logan\Desktop\Phylo\Actual\GardnellaAstral.tre 2> C:\Users\logan\Desktop\Phylo\Actual\GardnerellaAstral.log