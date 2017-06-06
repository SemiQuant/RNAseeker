#prepare data

#multiqc
# do a multi QC in folder
# file named "multiqc_report.html"

#rRNA
#luca's rRNA checker script
# folder namded "rRNA_check"

#dds
sampleFiles <- c("T004-B2SalBAL-PB-miRNA_S9.Homo_sapiens.HTSeq.counts", "T006-B2SalBAL-PBmiRNA_S12.Homo_sapiens.HTSeq.counts",
                 "T007-B2SalBAL-PBmiRNA_S29.Homo_sapiens.HTSeq.counts", "T009-B2SalBAL-PBmiRNA_S6.Homo_sapiens.HTSeq.counts",
                 "T013-B2SalBAL-PBmiRNA_S5.Homo_sapiens.HTSeq.counts", "T014-B2SalBAL-PBmiRNA_S29.Homo_sapiens.HTSeq.counts",
                 "T015-B2SalBAL-PBmiRNA1_S19.Homo_sapiens.HTSeq.counts", "T015-B2SalBAL-PLmiRNA2_S20.Homo_sapiens.HTSeq.counts",
                 "T023-B2SalBAL-PBmiRNA_S8.Homo_sapiens.HTSeq.counts", "T025-B2SalBAL-PBmiRNA_S11.Homo_sapiens.HTSeq.counts",
                 "T026-B2SalBAL-PBmiRNA_S14.Homo_sapiens.HTSeq.counts","T028-B2SalBAL-PBmiRNA_S24.Homo_sapiens.HTSeq.counts",
                 "T031-B2SalBAL-PBmiRNA_S7.Homo_sapiens.HTSeq.counts", "T032-B2SalBAL-PBmiRNA_S24.Homo_sapiens.HTSeq.counts",
                 "T033-B2SalBAL-PBmiRNA_S17.Homo_sapiens.HTSeq.counts", "T035-B2SalBAL-PBmiRNA_S1.Homo_sapiens.HTSeq.counts",
                 "T038-B2SalBAL-PBmiRNA_S26.Homo_sapiens.HTSeq.counts", "T039-B2SalBAL-PBmiRNA_S11.Homo_sapiens.HTSeq.counts",
                 "T041-B2SalBAL-PBmiRNA_S23.Homo_sapiens.HTSeq.counts", "T042-B2SalBAL-PBmiRNA1_S7.Homo_sapiens.HTSeq.counts",
                 "T042-B2SalBAL-PLmiRNA2_S8.Homo_sapiens.HTSeq.counts", "T047-B2SalBAL-PBmiRNA1_S3.Homo_sapiens.HTSeq.counts",
                 "T047-B2SalBAL-PLmiRNA2_S4.Homo_sapiens.HTSeq.counts", "T048-B2SalBAL-PBmiRNA_S5.Homo_sapiens.HTSeq.counts",
                 "T051-B2SalBAL-PBmiRNA_S13.Homo_sapiens.HTSeq.counts", "T054-B2SalBAL-PBmiRNA_S18.Homo_sapiens.HTSeq.counts")
# sub("(*.this).*","\\1",sampleFiles) # this deltes everything following 'this'
# sampleCondition <- data.frame(gsub(".Homo_sapiens.HTSeq.counts", "", sampleFiles)); colnames(sampleCondition) <- "File_name"


#first col sample name
#filename
#metadata columns
sampleTable <- data.frame(sampleName = gsub(".Homo_sapiens.HTSeq.counts", "", sampleFiles),
                          sampleFiles = sampleFiles,
           Group = c("LTBI", "LTBI", "prev", "prev", "prev", "prev", "prev", "prev", "TB", "prev", "prev", "TB", "prev", "prev", "prev", "sterilizing", "prev", "prev", "self", "LTBI", "LTBI", "prev", "prev", "prev", "recurrent", "recurrent"),
           rep = "1")


require(DESeq2)
# parallel = T
# library(BiocParallel)
# register(MulticoreParam(10))
# DESeq(dds, parallel=10)
# results(DESeq_dds, parallel = T)

dds_HTSeq_in <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                        directory = count_dir,
                                        design = ~ Group)

save(dds_HTSeq_in,
     file = "example_data/dds_HTSeq_in")

rlt_in <- rlog(dds_HTSeq_in, blind=T)
save(rlt_in,
     file = "example_data/rlt_in")

vst_in <- vst(dds_HTSeq_in, blind=T)
save(vst_in,
     file = "example_data/vst_in")

#dds end




bio <- read.table("bio.txt")
# require(NOISeq)
noi_dat <- readData(assay(dds_HTSeq_in), colData(dds_HTSeq_in), biotype = bio)
noi_dat_saturation <-  dat(Heart_data, k = 0, ndepth = 10, type = "saturation")











# require(refGenome)
# beg <- ensemblGenome()
# basedir(beg) <- "/Users/jdlim/Desktop/Temp/MA/"
# read.gtf(beg, "Homo_sapiens.GRCh38.86.gtf")
# # bio <- extractSeqids(beg, "^1$")
# bio <- data.frame(beg@ev$genes)[c(26,22)]

# bio_detect <- dat(Heart_data, type = "biodetection")
# mysaturation = dat(Heart_data, k = 0, ndepth = 7, type = "saturation")
# load("/Users/jdlim/Library/Mobile Documents/com~apple~CloudDocs/Bioinformatics/TB_HEART/TB_heart_anlaysis/mysaturation.dat")
# mysaturation_in <<- mysaturation
mycountsbio = dat(noi_dat, factor = NULL, type = "countsbio")
save(mycountsbio,
     file = "example_data/mycountsbio")
bio_detect <- dat(noi_dat, k = 0, type = "biodetection", factor = NULL)
save(bio_detect,
     file = "example_data/bio_detect")









