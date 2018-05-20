library(Biobase)

# Preparation -------------------------------------------------------------

fileNames = gsub(".gz", "", as.character(sample_tab$File.Name))
system("cat ./metadata/gdc_sample_sheet.2018-05-19.tsv | grep FPKM.txt.gz > ./metadata/gdc_sample_sheet.FPKM.tsv")
system("cat ./metadata/gdc_sample_sheet.2018-05-19.tsv | grep FPKM-UQ.txt.gz > ./metadata/gdc_sample_sheet.FPKMUQ.tsv")
system("cat ./metadata/gdc_sample_sheet.2018-05-19.tsv | grep htseq.counts.gz > ./metadata/gdc_sample_sheet.htseqcounts.tsv")
system("wc -l ./metadata/*")

# get counts table --------------------------------------------------------
    
sample_tab = read.delim("./metadata/gdc_sample_sheet.htseqcounts.tsv", header = F)
head(sample_tab)    
sample_map = as.character(sample_tab$V7)
names(sample_map) = gsub(".gz", "", as.character(sample_tab$V2))
head(sample_map)

files = list.files("./data_raw/", pattern = "htseq.counts", full.names = F)

data = read.table(paste0("./data_raw/", files[1]))
head(data)

for (file in files[2:length(files)]) {
    print (file)
    tmp = read.table(paste0("./data_raw/", file))
    
    if (all (tmp$V1 == data$V1)) {
        data = cbind (data, tmp$V2)
    } else {
        break ("File format wrong.")
    }
}

colnames(data) = c("Gene", sample_map[files])
head(data)

data[1:10, 1:10]
P1RCC_countsTable = data[, 2:322]
rownames(P1RCC_countsTable) = as.character(data$Gene)
P1RCC_countsTable[1:10, 1:10]

ESET_countsTable = new("ExpressionSet", exprs = as.matrix(P1RCC_countsTable))
rownames(sample_tab) = gsub(".gz", "", as.character(sample_tab$V2))
sample_tab_sort = sample_tab[files, ]
head(sample_tab)
head(sample_tab_sort)
pData(ESET_countsTable) = sample_tab_sort[, c("V7", "V8")]


# get FPKM table ----------------------------------------------------------

sample_tab = read.delim("./metadata/gdc_sample_sheet.FPKM.tsv", header = F)
head(sample_tab)    
sample_map = as.character(sample_tab$V7)
names(sample_map) = gsub(".gz", "", as.character(sample_tab$V2))
head(sample_map)

files = list.files("./data_raw/", pattern = "FPKM.txt", full.names = F)

data = read.table(paste0("./data_raw/", files[1]))
head(data)

for (file in files[2:length(files)]) {
    print (file)
    tmp = read.table(paste0("./data_raw/", file))
    
    if (all (tmp$V1 == data$V1)) {
        data = cbind (data, tmp$V2)
    } else {
        break ("File format wrong.")
    }
}

colnames(data) = c("Gene", sample_map[files])
head(data)

data[1:10, 1:10]
P1RCC_FPKM = data[, 2:322]
rownames(P1RCC_FPKM) = as.character(data$Gene)
P1RCC_FPKM[1:10, 1:10]

ESET_FPKM = new("ExpressionSet", exprs = as.matrix(P1RCC_FPKM))
rownames(sample_tab) = gsub(".gz", "", as.character(sample_tab$V2))
sample_tab_sort = sample_tab[files, ]
head(sample_tab)
head(sample_tab_sort)
pData(ESET_FPKM) = sample_tab_sort[, c("V7", "V8")]



# delete one sample -------------------------------------------------------

del = which (pData(ESET_countsTable)[, 2] == "Additional - New Primary")
ESET_countsTable_Sel = ESET_countsTable[, -del]
saveRDS(ESET_countsTable_Sel, "./data_proc/ESET_countsTable.rds")

del = which (pData(ESET_FPKM)[, 2] == "Additional - New Primary")
ESET_FPKM_Sel = ESET_FPKM[, -del]
saveRDS(ESET_FPKM_Sel, "./data_proc/ESET_FPKM.rds")

