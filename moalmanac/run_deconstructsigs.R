#!/usr/bin/Rscript
library("deconstructSigs")

args = commandArgs(trailingOnly=TRUE)
patient_id = args[1]
snv_handle = args[2]
sample = args[3]
ref = args[4]
alt = args[5]
chr = args[6]
pos = args[7]

maf = read.csv(snv_handle, sep = '\t', comment.char = '#')
cols = c(sample, ref, alt, chr, pos)
maf <- maf[colnames(maf) %in% cols]

maf$Tumor_Sample_Barcode <- sapply(maf$Tumor_Sample_Barcode, as.factor)
maf$Reference_Allele <- sapply(maf$Reference_Allele, as.factor)
maf$Tumor_Seq_Allele2 <- sapply(maf$Tumor_Seq_Allele2, as.factor)
maf$Chromosome <- sapply(maf$Chromosome, as.factor)

unique.samples = unique(maf$Tumor_Sample_Barcode)

sigs.input <- mut.to.sigs.input(mut.ref = maf,
    sample.id = sample, chr = chr,
    pos = pos, ref = ref,
    alt = alt)
    
temp.filename <- paste(patient_id, ".sigs.context.txt", sep = "")
write.table(sigs.input, file = temp.filename, sep = '\t', row.names = FALSE)

for (sample_ in unique.samples) {
    output.sigs <- whichSignatures(tumor.ref = sigs.input,
        signatures.ref = signatures.cosmic, sample.id = sample_, 
        context = TRUE, tri.counts.method = 'default')
    
    temp.filename = paste(patient_id, ".sigs.cosmic.txt", sep = "")
    write.table(output.sigs, file = temp.filename, sep = '\t', row.names = FALSE)
}
