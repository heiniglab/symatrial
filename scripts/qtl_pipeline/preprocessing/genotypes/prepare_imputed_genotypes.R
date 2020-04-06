# Prepare genotypes and snp locations for imputed data file

path <- "/home/icb/ines.assum/projects/symAtrial_QTL/scripts/"
setwd(path)
source("helper/helper.R")
local=F

if(F){
  library(data.table)
  # imputed dosage file, not filtered for maf and hwe
  geno.file <- paste0(get.path("data", local), "genotype/AFHRI_B_imputed.txt")
  map.file <- paste0(get.path("data", local), "genotype/AFHRI_B_imputed_1.map")
  
  if(!file.exists(map.file)){
    geno.imp <- fread(geno.file, sep="\t", header=T)
    table(geno.imp$typed_marker_id)
    geno.ind <- colnames(geno.imp)[-c(1:6)]
    
    texpr <- readRDS(paste0(get.path("data", local), "AFHRI_B/AFHRI_B_transcriptomics_QC_symbol.RDS"))
    pexpr <- readRDS(paste0(get.path("data", local), "AFHRI_B/AFHRI_B_proteomics_QC_symbol.RDS"))
    
    gt.ind <- intersect(geno.ind, sub(".RAA", "", colnames(texpr)))
    gp.ind <- intersect(geno.ind, sub(".RAA", "", colnames(pexpr)))
    gres.ind <- intersect(geno.ind, intersect(sub(".RAA", "", colnames(texpr)),
                                              sub(".RAA", "", colnames(pexpr))))
    
    write.table(data.frame(IID=gt.ind, FID=gt.ind),
                paste0(get.path("data", local), "genotype/ind_genes_transcripts.txt"),
                quote = F, row.names = F)
    write.table(data.frame(IID=gp.ind, FID=gp.ind),
                paste0(get.path("data", local), "genotype/ind_genes_proteins.txt"),
                quote = F, row.names = F)
    write.table(data.frame(IID=gres.ind, FID=gres.ind),
                paste0(get.path("data", local), "genotype/ind_genes_residuals.txt"),
                quote = F, row.names = F)
    
    # ???rs2875098:15000001:T:C???, ???rs9478072:170000001:A:G???
    # which(duplicated(geno.imp$imputed_marker_id))
    # [1] 2441685 4876505
    # test <- geno.imp[which(geno.imp$imputed_marker_id %in% c("rs2875098:15000001:T:C", "rs9478072:170000001:A:G")), ]
    # test
    geno.imp <- geno.imp[!duplicated(geno.imp$imputed_marker_id), ]
    map <- geno.imp[, c("chrom", "imputed_marker_id", "position", "allele1", "allele2")]
    write.table(map,
                file=map.file,
                quote=F, sep="\t", row.names = F)
  }
  
  
  # system(paste0("cd /home/icb/ines.assum/projects/symAtrial_QTL/data/current/genotype \n",
  #               "/home/icb/ines.assum/ext_tools/plink/plink --bfile AFHRI-B_imputed --extract AFHRI_B_imputed_1.map --keep ind_genes_transcripts.txt --hwe 0.0001 --hardy --out stats_trans"))
  # 
  # system(paste0("cd /home/icb/ines.assum/projects/symAtrial_QTL/data/current/genotype \n",
  #               "/home/icb/ines.assum/ext_tools/plink/plink --bfile AFHRI-B_imputed --extract AFHRI_B_imputed_1.map --keep ind_genes_proteins.txt --hwe 0.0001 --hardy --out stats_prot"))
  # 
  # system(paste0("cd /home/icb/ines.assum/projects/symAtrial_QTL/data/current/genotype \n",
  #               "/home/icb/ines.assum/ext_tools/plink/plink --bfile AFHRI-B_imputed --extract AFHRI_B_imputed_1.map --keep ind_genes_residuals.txt --hwe 0.0001 --hardy --out stats_res"))
  
  # lmiss <- read.table("/home/icb/ines.assum/projects/symAtrial_QTL/data/current/genotype/stats_imputed.lmiss",
  #                     h=T,
  #                     stringsAsFactors = F)
  # tab <- table(lmiss$N_MISS)
  
  # hwe.trans <- read.table("/home/icb/ines.assum/projects/symAtrial_QTL/data/current/genotype/stats_trans.hwe",
  #                         h=T, stringsAsFactors = F)
  # hwe.prot <- read.table("/home/icb/ines.assum/projects/symAtrial_QTL/data/current/genotype/stats_prot.hwe",
  #                        h=T, stringsAsFactors = F)
  # hwe.res <- read.table("/home/icb/ines.assum/projects/symAtrial_QTL/data/current/genotype/stats_res.hwe",
  #                       h=T, stringsAsFactors = F)
  # hwe.all <- rbind(data.frame(hwe.trans[hwe.trans$TEST=="ALL", ], SET="trans", stringsAsFactors = F),
  #                  data.frame(hwe.prot[hwe.prot$TEST=="ALL", ], SET="prot", stringsAsFactors = F),
  #                  data.frame(hwe.res[hwe.res$TEST=="ALL", ], SET="residuals", stringsAsFactors = F))
  # 
  # colnames(hwe.all) <- c("CHR", "SNP", "TEST", "A1", "A2", "P11/P12/P22", "O.HET.", "E.HET.", "P", "SET")
  # head(hwe.all)
  # write.table(hwe.all,
  #             file ="/home/icb/ines.assum/projects/symAtrial_QTL/data/current/genotype/stats_imputed_splitslash.hwe",
  #             col.names = T, row.names = F, sep = "/", quote = F)
  # hwe.all <- read.table(file ="/home/icb/ines.assum/projects/symAtrial_QTL/data/current/genotype/stats_imputed_splitslash.hwe",
  #                       h = T, sep = "/", stringsAsFactors = F)
  # saveRDS(hwe.all,
  #         file = paste0(get.path("data", local), "genotype/stats_imputed.hwe.RDS"))
  
  if(!(file.exists(paste0(get.path("data", local), "genotype/recode_snps.RDS")) &
       file.exists(paste0(get.path("data", local), "genotype/filtered_snps.RDS")))){
      hwe.all <- readRDS(paste0(get.path("data", local), "genotype/stats_imputed.hwe.RDS"))
      filtered <- intersect(intersect(hwe.all[hwe.all$SET=="trans", "SNP"],
                                      hwe.all[hwe.all$SET=="prot", "SNP"]),
                            hwe.all[hwe.all$SET=="residuals", "SNP"])
      saveRDS(filtered,
              file=paste0(get.path("data", local), "genotype/filtered_snps.RDS"))
      
      # don't do anything for SNPs with no homozygote minor allele genotypes
      hwe <- hwe.all[hwe.all$P11>0, ]
      table(hwe$P11)
      # select SNPs with less or equal 3 individuals homozygote minor allele genotypes
      hwe <- hwe[hwe$P11<=3, ]
      table(hwe$P11)
      
      # hwe.0 <- hwe[hwe$P11==0,]
      # hwe.0 <- hwe.0[!duplicated(hwe.0$SNP), ]
      
      recode <- unique(hwe$SNP)
      length(recode)
      saveRDS(recode,
              file = paste0(get.path("data", local), "genotype/recode_snps_new.RDS"))
      saveRDS(hwe,
              file = paste0(get.path("data", local), "genotype/stats_imputed_recode_snps.hwe.RDS"))
      write.table(hwe[!duplicated(hwe$SNP), ],
                  file = paste0(get.path("data", local), "genotype/stats_imputed_recode_snps.hwe.txt"),
                  row.names=F, col.names=T, quote=F, sep="\t")
      dim(hwe)
      dim(hwe[!duplicated(hwe$SNP), ])
      colnames(hwe)
      
  }else{
    recode <- readRDS(paste0(get.path("data", local), "genotype/recode_snps_new.RDS"))
    filtered <- readRDS(paste0(get.path("data", local), "genotype/filtered_snps.RDS"))
  }
  
  print(geno.file)
  geno.imp <- fread(geno.file, sep="\t", header=T)
  dim(geno.imp)
  geno.imp <- data.frame(geno.imp[!duplicated(geno.imp$imputed_marker_id), ])
  rownames(geno.imp) <- geno.imp$imputed_marker_id
  dim(geno.imp)
  
  geno.imp2 <- geno.imp[geno.imp$imputed_marker_id %in% filtered, ]
  dim(geno.imp2)
  recode.rows <- geno.imp2$imputed_marker_id[geno.imp2$imputed_marker_id %in% recode]
  geno.imp3 <- as.matrix(geno.imp2[recode.rows, -(1:6)], row.names = recode.rows)
  
  geno.imp3.map <- as.matrix(geno.imp2[recode.rows, 1:6], row.names = recode.rows)
  coding.file <- paste0(get.path("data", local), "genotype/stats_imputed_recode_snps.hwe.txt")
  snp.file <- tempfile()
  write.table(recode, file=snp.file,
              row.names = F, col.names = F, sep="\n", quote=F)
  #coding = fread(paste("fgrep -f ", snp.file, " -w ", coding.file, sep=""), header=F)
  coding = fread(paste("awk 'FNR==NR { a[$1]; next } $2 in a'", snp.file, coding.file, sep=" "), header=F)
  setnames(coding, colnames(read.csv(coding.file, sep="\t", nrows=5)))
  coding = as.data.frame(coding)
  
  colnames(geno.imp3.map)
  # [1] "chrom"             "typed_marker_id"   "imputed_marker_id" "position"         
  # [5] "allele1"           "allele2" 
  colnames(coding)
  # [1] "CHR"    "SNP"    "TEST"   "A1"     "A2"     "P11"    "P12"    "P22"   
  # [9] "O.HET." "E.HET." "P"      "SET"
  
  match.alleles <- merge(geno.imp3.map, coding,
                         by.x=c("imputed_marker_id"), by.y="SNP",
                         all.x=T, sort=F)
  head(match.alleles)
  rownames(match.alleles) <- match.alleles$imputed_marker_id
  match.alleles <- match.alleles[rownames(geno.imp3), ]
  identical(rownames(geno.imp3), rownames(match.alleles))
  
  # split into two datasets: minor allele first / minor allele last
  match.alleles$same = toupper(match.alleles[,"allele1"]) == toupper(match.alleles[,"A1"])
  table(match.alleles$same)
  head(match.alleles)
  t(geno.imp3[rownames(match.alleles)[1:2], ])
  
  geno.imp2[1:10, 1:10]
  geno.imp3[1:10, 1:10]
  dim(geno.imp3)
  geno.imp3[match.alleles$same, ][geno.imp3[match.alleles$same, ]<0.5] <- 1
  geno.imp3[!match.alleles$same, ][geno.imp3[!match.alleles$same, ]>1.5] <- 1

  dim(geno.imp3)
  geno.imp2[recode.rows, -(1:6)] <- geno.imp3
  
  saveRDS(match.alleles,
          file=paste0(get.path("data", local), "genotype/AFHRI_B_imputed_preprocessed_genotypes_allele_coding.RDS"))
  saveRDS(geno.imp2,
          file="/home/icb/ines.assum/projects/symAtrial_QTL/data/current/genotype/AFHRI_B_imputed_recoded_genotypes_new.RDS")
  
  write.table(geno.imp2[, -c(1:2, 4:6)],
              file=paste0(get.path("data", local), "genotype/AFHRI_B_imputed_preprocessed_genotypes.txt"),
              quote=F, sep="\t", row.names = F)
  write.table(geno.imp2[, c("imputed_marker_id", "chrom", "position")],
              #file="/home/icb/ines.assum/projects/symAtrial_QTL/data/current//annotation/snplocs_imputed.txt"
              file=paste0(get.path("snplocs imputed", local)),
              quote=F, sep="\t", row.names = F)
  
  # map file for import into GenABEL
  write.table(geno.imp2[, c("chrom", "imputed_marker_id", "position")],
              file=paste0(get.path("data", local), "genotype/AFHRI_B_imputed_2.map"),
              quote=F, sep="\t", row.names = F)
  # map file to recode alleles
  rec.all <- cbind(geno.imp2[, c("imputed_marker_id", "allele1", "allele2")],
                   "A", "B")
  write.table(rec.all,
              file=paste0(get.path("data", local), "genotype/AFHRI_B_imputed_recode_alleles.map"),
              quote=F, sep="\t", row.names = F)
}


setwd("/home/icb/ines.assum/projects/symAtrial_QTL/scripts")
# setwd("/Users/ines/Documents/ICB/PhD/projects/symAtrial_QTL/scripts")

source("helper/helper.R")
local=F

data_path <- get.path("dataset", local)

expression_p <- readRDS(paste0(data_path, "AFHRI_B_proteomics_QC_symbol.RDS"))
colnames(expression_p) <- sub(".RAA","", colnames(expression_p))

expression_t <- readRDS(paste0(data_path, "AFHRI_B_transcriptomics_QC_symbol.RDS"))
colnames(expression_t) <- sub(".RAA", "", colnames(expression_t))

residuals_t <- readRDS(paste0(data_path, "residuals_trans.RDS"))
colnames(residuals_t) <- sub(".RAA","", colnames(residuals_t))

residuals_p <- readRDS(paste0(data_path, "residuals_prot.RDS"))
colnames(residuals_p) <- sub(".RAA", "", colnames(residuals_p))

ratios <- readRDS(paste0(data_path, "prot_trans_ratios.RDS"))
colnames(ratios) <- sub(".RAA","", colnames(ratios))

locations_p <- readRDS(paste0(get.path("locations", local), "protein_locations_ensembl_tss.RDS"))
locations_t <- readRDS(paste0(get.path("locations", local), "transcript_locations_ensembl_tss.RDS"))

source("preprocessing/data_parsing/cut_tables.R")

new_expression_p <- cut.locations(expression_p, locations_p)
new_expression_t <- cut.locations(expression_t, locations_t)
new_expression_res_t <- cut.locations(residuals_t, locations_t)
new_expression_res_p <- cut.locations(residuals_p, locations_p)
new_expression_ratios <- cut.locations(ratios, locations_t)

library(data.table)

# processed_genotype_file <- get.path("base genotype")
processed_genotype_file <- get.path("base genotype imputed", local)

genotype <- fread(processed_genotype_file, sep="\t", header=T)
rownames(genotype) <- genotype[[1]]
genotype[,1] <- NULL


# drop samples

results_p <- common.samples(new_expression_p, genotype)
new_gen_p <- results_p[[1]]
new_exp_p <- results_p[[2]]

results_t <- common.samples(new_expression_t, genotype)
new_gen_t <- results_t[[1]]
new_exp_t <- results_t[[2]]

results_res_t <- common.samples(new_expression_res_t, genotype)
new_gen_res_t <- results_res_t[[1]]
new_exp_res_t <- results_res_t[[2]]

results_res_p <- common.samples(new_expression_res_p, genotype)
new_gen_res_p <- results_res_p[[1]]
new_exp_res_p <- results_res_p[[2]]

results_ratios <- common.samples(new_expression_ratios, genotype)
new_gen_ratios <- results_ratios[[1]]
new_exp_ratios <- results_ratios[[2]]



geno_path <- get.path("genotype", local)

geno_t <- paste0(geno_path, "genotype_imputed_common_samples_t_raa.txt")
geno_p <- paste0(geno_path, "genotype_imputed_common_samples_p_raa.txt")
geno_res_t <- paste0(geno_path, "genotype_imputed_common_samples_res_t_raa.txt")
geno_res_p <- paste0(geno_path, "genotype_imputed_common_samples_res_p_raa.txt")
geno_ratios <- paste0(geno_path, "genotype_imputed_common_samples_ratios.txt")

write.table(new_gen_t, geno_t, sep="\t", quote=F, col.names=NA)
write.table(new_gen_p, geno_p, sep="\t", quote=F, col.names=NA)
write.table(new_gen_res_t, geno_res_t, sep="\t", quote=F, col.names=NA)
write.table(new_gen_res_p, geno_res_p, sep="\t", quote=F, col.names=NA)
write.table(new_gen_ratios, geno_ratios, sep="\t", quote=F, col.names=NA)

