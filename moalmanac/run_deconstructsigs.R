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
folder = args[8]

maf = read.csv(snv_handle, sep = '\t', comment.char = '#')
names(maf) <- tolower(names(maf))
cols = c(sample, ref, alt, chr, pos)
maf <- maf[colnames(maf) %in% cols]

maf$tumor_sample_barcode <- sapply(maf$tumor_sample_barcode, as.factor)
maf$reference_allele <- sapply(maf$reference_allele, as.factor)
maf$tumor_seq_allele2 <- sapply(maf$tumor_seq_allele2, as.factor)
maf$chromosome <- sapply(maf$chromosome, as.factor)

unique.samples = unique(maf$tumor_sample_barcode)

sigs.input <- mut.to.sigs.input(mut.ref = maf,
    sample.id = sample, chr = chr,
    pos = pos, ref = ref,
    alt = alt)
    
temp.filename <- paste(folder, patient_id, ".sigs.context.txt", sep = "")
write.table(sigs.input, file = temp.filename, sep = '\t', row.names = FALSE)

for (sample_ in unique.samples) {
    output.sigs <- whichSignatures(tumor.ref = sigs.input,
        signatures.ref = signatures.cosmic, sample.id = sample_, 
        context = TRUE, tri.counts.method = 'default')
    
    temp.filename = paste(folder, patient_id, ".sigs.cosmic.txt", sep = "")
    write.table(output.sigs, file = temp.filename, sep = '\t', row.names = FALSE)
}
