#Install edgeR
bash
cd 
git clone https://github.com/ctb/2017-ucdavis-igg201b.git
sudo Rscript --no-save ~/2017-ucdavis-igg201b/lab7/install-edgeR.R

#install salmon
cd
curl -L -O https://github.com/COMBINE-lab/salmon/releases/download/v0.8.0/Salmon-0.8.0_linux_x86_64.tar.gz
tar xzf Salmon-0.8.0_linux_x86_64.tar.gz
export PATH=$PATH:$HOME/Salmon-latest_linux_x86_64/bin

mkdir yeast
cd yeast

#download data (muts then wts)
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR458/ERR458500/ERR458500.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR458/ERR458501/ERR458501.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR458/ERR458502/ERR458502.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR458/ERR458689/ERR458689.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR458/ERR458701/ERR458701.fastq.gz

curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR458/ERR458493/ERR458493.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR458/ERR458494/ERR458494.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR458/ERR458495/ERR458495.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR459/ERR459145/ERR459145.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR459/ERR459155/ERR459155.fastq.gz

#Yeast reference
curl -O http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding.fasta.gz

#Index yeast transcriptome
salmon index --index yeast_orfs --type quasi --transcripts orf_coding.fasta.gz

#Run salmon
for i in *.fastq.gz
do
    salmon quant -i yeast_orfs --libType U -r $i -o $i.quant --seqBias --gcBias
done

#collect sample counts
curl -L -O https://github.com/ngs-docs/2016-aug-nonmodel-rnaseq/raw/master/files/gather-counts.py
python2 gather-counts.py

#edgeR
R
library("edgeR")

files <- c(
"ERR458493.fastq.gz.quant.counts",
"ERR458494.fastq.gz.quant.counts",
"ERR458495.fastq.gz.quant.counts",
"ERR459145.fastq.gz.quant.counts",
"ERR459155.fastq.gz.quant.counts",
"ERR458500.fastq.gz.quant.counts",
"ERR458501.fastq.gz.quant.counts",
"ERR458502.fastq.gz.quant.counts",
"ERR458689.fastq.gz.quant.counts",
"ERR458701.fastq.gz.quant.counts"
)

labels=c("WT1", "WT2", "WT3", "WT4", "WT5", "MUT1", "MUT2", "MUT3", "MUT4", "MUT5")

data <- readDGE(files)

print(data)

###

group <- c(rep("wt", 5), rep("mut", 5))

dge = DGEList(counts=data, group=group)
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)

# make an MA-plot

et <- exactTest(dge, pair=c("wt", "mut"))
etp <- topTags(et, n=100000)
etp$table$logFC = -etp$table$logFC
pdf("yeast-edgeR-MA-plot.pdf")
plot(
  etp$table$logCPM,
  etp$table$logFC,
  xlim=c(-3, 20), ylim=c(-12, 12), pch=20, cex=.3,
  col = ifelse( etp$table$FDR < .2, "red", "black" ) )
dev.off()

# plot MDS
pdf("yeast-edgeR-MDS.pdf")
plotMDS(dge, labels=labels)
dev.off()

# output CSV for 0-6 hr
write.csv(etp$table, "yeast-edgeR.csv")


#Commnents:
#The additional four samples yielded 3590 upregulated and 593 downregulated transcripts that were significantly DE
#Lab 8 gave us 548 upregulated and 3298 downregulated transcripts that were significantly DE
