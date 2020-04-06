#!/usr/bin/env Rscript

library("optparse")
options <- list(
    make_option(c("-r", "--reference"), type="character",
                default=NULL, 
                help="Reference haplotypes directory",
                metavar="file"))
parser <- OptionParser(usage = "%prog [options] prefix", option_list=options)
arguments <- parse_args(parser, positional_arguments = 1)

prefix <- arguments$args
ref.dir <- arguments$options$reference

if (is.null(ref.dir)) {
    cat("--reference option needs to be specified!\n")
    q("no")
}

cat("Aligning genotypes in plink file\n", prefix,
    "\nto reference haplotypes in \n", ref.dir, "\n")

## load our allele frequencies
our.snps <- read.table(paste0(prefix, ".frq"), header=T, stringsAsFactors=F)
## and map positions
our.map <- read.table(paste0(prefix, ".map"), header=T, stringsAsFactors=F)
colnames(our.map) <- c("CHR", "SNP", "cM", "position")
our.snps <- merge(our.snps, our.map)

## go through them by chromosome
aligned <- NULL
ignore.chr <- c("0", "23", "24", "25", "26")
for (chr in setdiff(unique(our.snps$CHR), ignore.chr)) {
  cat("chr:", chr, "\n")
  our.snps.chr <- our.snps[our.snps$CHR == chr,]
  ## has alleles A1 and A2 and the frequency of A1: we call A1 the reference allele

  ref.snps <- read.table(paste0(ref.dir, "/1000GP_Phase3_chr", chr, ".legend.gz"), header=T, stringsAsFactors=F)
  ## has columns a0 and a1 for the alleles and allele freq of a1 in 1kg populations
  ref.snps <- data.frame(ref.snps, rsid=sapply(strsplit(ref.snps$id, ":"), "[", 1), stringsAsFactors=F)

  ## match the snps by position
  merged <- merge(our.snps.chr, ref.snps, by="position")
  merged <- data.frame(merged, nonref.AF=1-merged$MAF)

  ## compare the alleles
  same <- with(merged, A1 == a0 & A2 == a1)
  summary(with(merged[same,], nonref.AF - EUR))

  ## flag problemantic cases with allele freq deviation of more than 10%
  af.mismatch <- rep(NA, nrow(merged))
  af.mismatch[same] <- with(merged[same,], abs(nonref.AF - EUR) > 0.1)

  ## check if we need to swap alleles
  swap <- !same & with(merged, A2 == a0 & A1 == a1)
  af.mismatch[swap] <- with(merged[swap,], abs(MAF - EUR) > 0.1)

  ## check if we need to flip strands
  complement <- c(A="T", C="G", G="C", T="A")
  flip <- !same & !swap & with(merged, complement[A1] == a0 & complement[A2] == a1)
  flip[is.na(flip)] <- FALSE ## this can hapen for indels
  af.mismatch[flip] <- with(merged[flip,], abs(nonref.AF - EUR) > 0.1)

  ## check if we need to flip strands and swap alleles
  flip.swap <- !same & !swap & !flip & with(merged, complement[A2] == a0 & complement[A1] == a1)
  flip.swap[is.na(flip.swap)] <- FALSE ## this can hapen for indels
  af.mismatch[flip.swap] <- with(merged[flip.swap,], abs(MAF - EUR) > 0.1)
  
  aligned <- rbind(aligned, data.frame(merged, same.encoding=same, needs.swap=swap, needs.flip=flip, needs.flip.and.swap=flip.swap, af.mismatch))
}

write.table(aligned, paste0(prefix, "_aligned_map.txt"), sep="\t", quote=F, row.names=F)

## in the end we write out a list of snps that has AF mismatches to be removed
clean <- aligned[which(with(aligned, !is.na(af.mismatch) & !af.mismatch)),]
clean <- clean[with(clean, !(!same.encoding & !needs.swap & !needs.flip & !needs.flip.and.swap)),]

dup <- duplicated(clean$rsid)
dupped.snp <- clean$SNP[dup]
clean <- clean[!dup,]

dup <- duplicated(clean$SNP)
dupped <- clean[clean$SNP %in% clean$SNP[dup],]
clean <- clean[!dup,]

remove <- setdiff(our.snps$SNP, clean$SNP)
remove <- union(remove, dupped.snp)
cat(remove, sep="\n", file=paste0(prefix, "_snps_with_AF_mismatch.txt"))

## and a file with the proper reference alles to run plink recode
ref.table <- clean[, c("SNP", "a0")]
ref.table <- ref.table[!duplicated(ref.table$SNP),]
write.table(ref.table, sep="\t", quote=F, row.names=F, file=paste0(prefix, "_reference_alleles.txt"))

## we also need a flip list
flip.list <- with(clean, SNP[needs.flip | needs.flip.and.swap])
cat(flip.list, sep="\n", file=paste0(prefix, "_flip_strands.txt"))

## also rename markers to the names in the reference
write.table(clean[,c("SNP", "rsid")], file=paste0(prefix, "_rename_snps.txt"), sep="\t", quote=F, row.names=F, col.names=F)

## finally we also create some summary statistics
print(with(aligned, table(same.encoding, needs.swap)))

print(with(aligned, table(CHR, af.mismatch)))

cat("removed", length(remove), "of", nrow(our.snps), "snps\n")
