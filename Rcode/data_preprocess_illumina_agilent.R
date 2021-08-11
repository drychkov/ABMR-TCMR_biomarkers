library("limma")
library("GEOquery")
library("biomaRt")
library("foreach")

# Helpful functions -----
fixGeoSampleNames = function(sampleNames) {
    sampleNames = paste0(toupper(sampleNames), '_')
    regexResult = regexpr('^GSM[0-9]+[^0-9]', sampleNames)
    sampleNamesNew = mapply(function(sampleName, matchLength) substr(sampleName, 1, matchLength-1),
                            sampleNames, attr(regexResult, 'match.length'))
    return(sampleNamesNew)}

getGeneProbeMappingDirect = function(featureDf, geneColname, probeColname='ID') {
    mapping = featureDf[,c(probeColname, geneColname)]
    mapping = mapping[apply(mapping, MARGIN=1, function(x) all(!is.na(x) & x!='')),]
    mapping = data.frame(lapply(mapping, as.character), stringsAsFactors=FALSE)
    colnames(mapping) = c('probeSet', 'geneId')
    return(mapping)}

calcExprsByGeneEmat = function(emat, mapping) {
    cat("Mapping probes to genes...", "\n")
    mapping = mapping[mapping$probeSet %in% rownames(emat),]
    geneIds = unique(mapping[,'geneId'])
    exprsByGene = matrix(nrow = length(geneIds), 
                         ncol = ncol(emat), 
                         dimnames = list(geneIds, colnames(emat)))
    for (geneId in geneIds) {
        probeSet = mapping[mapping[,'geneId'] == geneId, 'probeSet']
        if (length(probeSet) == 0) {
            next
        }
        exprsTmp = emat[probeSet,, drop = F]
        if (nrow(exprsTmp) == 1) {
            exprsByGene[geneId,] = exprsTmp
        } else {
            exprsByGene[geneId,] = rowMeans(t(exprsTmp), na.rm = TRUE)
        }
    }
    # exprsByGene = exprsByGene[complete.cases(exprsByGene),]
    return(exprsByGene)}


# Loading BiomRt database ----
ensembl = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# View(listDatasets(ensembl))
# View(listFilters(ensembl))
q1 = listAttributes(ensembl)
View(q1[grep("affy|agilent|illumina", q1$name),])


# Processing datasets ----
parentFolderPath = "../data_in"

## GSE43974 -----
studyName = "GSE43974"
q = read.ilmn(files = "GSE43974_non-normalized_data.txt.gz", 
              path = paste0(parentFolderPath),
              probeid = "PROBE_ID", expr = "AVG_Signal", other.columns = "Detection")
y = neqc(q)
mapping = getBM(attributes = c("illumina_humanht_12_v4", "entrezgene_id"), 
                mart = ensembl,
                verbose = F,
                uniqueRows = T)
mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "entrezgene_id", 
                                    probeColname = "illumina_humanht_12_v4")

exprsByGene = calcExprsByGeneEmat(y$E, mapping)
sum(is.na(exprsByGene))
eset <- new("ExpressionSet", exprs = exprsByGene)

esetGEO = getGEO("GSE43974")[[1]]

names = esetGEO$geo_accession[match(colnames(eset), esetGEO$description)]
colnames(eset) = names

saveRDS(eset, "../data_out/GSE43974.rds")

affy::plotDensity(log2(q$E))
affy::plotDensity(y$E)
affy::plotDensity(exprsByGene)


## GSE47097-----
studyName = "GSE47097"
q = read.ilmn(files="GSE47097_non_normalized.txt.gz", path = paste0(parentFolderPath, "/"),
              probeid = "ID_REF", expr = "SAMPLE", other.columns = "Detection")
y = neqc(q)

esetGEO = getGEO(filename="../../Data/Human/kidney/GEO/GSE47097_series_matrix.txt.gz")

mapping = getBM(attributes = c("illumina_humanref_8_v3", "entrezgene"), 
                mart = ensembl,
                verbose = F,
                uniqueRows = T)
mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "entrezgene", 
                                    probeColname = "illumina_humanref_8_v3")

exprsByGene = calcExprsByGeneEmat(y$E, mapping)
sum(is.na(exprsByGene))
eset <- new("ExpressionSet", exprs = exprsByGene)

names = esetGEO$geo_accession[match(paste0("SAMPLE", colnames(q)), esetGEO$description)]
colnames(eset) = names

saveRDS(eset, "GSE47097.rds")


## GSE52694 -------
studyName = "GSE52694"
x = read.ilmn(files = "GSE52694_non_normalized.txt.gz", path = paste0(parentFolderPath, "/"),
              probeid = "ID_REF", expr = "BORDERLINE_PRAGUE_")
y = backgroundCorrect.matrix(x$E, method = "normexp", normexp.method = "mle", offset = 16)
y = log2(y)
y = normalizeQuantiles(y)

esetGEO = getGEO(filename = "../../Data/Human/kidney/GEO/GSE52694_series_matrix.txt.gz")

mapping = getBM(attributes = c("illumina_humanht_12_v4", "entrezgene"), 
                mart = ensembl,
                verbose = F,
                uniqueRows = T)
mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "entrezgene", 
                                    probeColname = "illumina_humanht_12_v4")

exprsByGene = calcExprsByGeneEmat(y, mapping)
eset <- new("ExpressionSet", exprs = exprsByGene)

match(paste0("BORDERLINE_PRAGUE_", colnames(eset)), as.character(esetGEO$description))

names = esetGEO$geo_accession[match(paste0("BORDERLINE_PRAGUE_", colnames(eset)), 
                                    as.character(esetGEO$description))]
colnames(eset) = names

saveRDS(eset, "GSE52694.rds")


## GSE65326 ------
studyName = "GSE65326"
q = read.ilmn(files="GSE65326_non-normalized.txt.gz", path = paste0(parentFolderPath),
              probeid = "ID_REF", expr = "57",  other.columns = "Detection")
esetGEO = getGEO(filename = "../../Data/Human/kidney/GEO/GSE65326_series_matrix.txt.gz")

y = neqc(q)
mapping = getBM(attributes = c("illumina_humanht_12_v4", "entrezgene"), 
                mart = ensembl,
                verbose = F,
                uniqueRows = T)
mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "entrezgene", 
                                    probeColname = "illumina_humanht_12_v4")

exprsByGene = calcExprsByGeneEmat(y$E, mapping)
sum(is.na(exprsByGene))
exprsByGene = exprsByGene[complete.cases(exprsByGene),]

eset <- new("ExpressionSet", exprs = exprsByGene)

names = esetGEO$geo_accession[match(paste0("57", colnames(eset)), esetGEO$title)]
colnames(eset) = names

saveRDS(eset, "GSE65326.rds")


## GSE69677 ------
studyName = "GSE69677"
q = read.ilmn(files = "GSE69677_non-normalized.txt.gz", path = paste0(parentFolderPath, "/"),
              probeid = "PROBE_ID", expr = "AVG_Signal",  other.columns = "Detection")
esetGEO = getGEO(filename = "../../Data/Human/kidney/GEO/GSE69677_series_matrix.txt.gz")

y = neqc(q)

mapping = read.csv("../../Data/Human/kidney/GEO/HumanHT-12_V4_0_R2_15002873_B_WGDASL.txt.TXT", 
                     skip = 8, sep = "\t")

mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "Entrez_Gene_ID", 
                                    probeColname = "Probe_Id")

exprsByGene = calcExprsByGeneEmat(y$E, mapping)
sum(is.na(exprsByGene))
exprsByGene = exprsByGene[complete.cases(exprsByGene),]

eset <- new("ExpressionSet", exprs = exprsByGene)

eset = eset[,colnames(eset)[(colnames(eset) %in% esetGEO$description)]]
names = esetGEO$geo_accession[match(colnames(eset), esetGEO$description)]
colnames(eset) = names

saveRDS(eset, "GSE69677.rds")



## GSE60807 --------
studyName = "GSE60807"

x = read.maimages(files = dir(paste0(parentFolderPath, "/", studyName)), 
                  path = paste0(parentFolderPath, "/", studyName), source = "agilent", green.only = TRUE)
      
colnames(x) = unname(fixGeoSampleNames(colnames(x)))

y = neqc(x, status = x$genes$ControlType, negctrl = -1, regular = 0)
rownames(y$E) = y$genes$ProbeName

mapping = getBM(attributes = c("agilent_wholegenome_4x44k_v1", "entrezgene"), 
                mart = ensembl,
                verbose = F,
                uniqueRows = T)
mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "entrezgene", 
                                    probeColname = "agilent_wholegenome_4x44k_v1")

exprsByGene = calcExprsByGeneEmat(y$E, mapping)
sum(is.na(exprsByGene))
eset <- new("ExpressionSet", exprs = exprsByGene)

saveRDS(eset, "GSE60807.rds")


## GSE10419 -----
q = getGEO(filename = "../../Data/Human/kidney/GEO/GSE10419_family.soft.gz")
a = foreach(name = names(q@gsms)) %do% {emat = q@gsms[name]}[[1]]
b = lapply(a, function(x) q@gsms)[[1]]
c = lapply(b, function(x) x@dataTable@table[,6])
y = do.call(cbind, c)

meta = getGEO("GSE10419")[[1]]
rownames(y) = meta@featureData@data$SPOT_ID

emat = backgroundCorrect.matrix(y, method = "normexp", normexp.method = 'mle', offset = 16)
emat = log2(emat)
emat = normalizeQuantiles(emat)

mapping = meta@featureData@data[, c("SPOT_ID", "GENE")]
mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "GENE", 
                                    probeColname = "SPOT_ID")

exprsByGene = calcExprsByGeneEmat(emat, mapping)
sum(is.na(exprsByGene))
eset <- new("ExpressionSet", exprs = exprsByGene)

saveRDS(eset, "GSE10419.rds")



## GSE343 --------
studyName = "GSE343"
q = getGEO(filename = "../../Data/Human/kidney/GEO/GSE343_family.soft.gz")
pl = foreach(name = names(q@gsms)) %do% {emat = q@gsms[name]}[[1]]
pl = lapply(1:length(pl), function(x) pl[[x]]@header$platform_id)

a = sapply(names(q@gsms), function(x) q@gsms[[x]]@dataTable@table$CH2I_MEAN)
y = as.matrix(as.data.frame(a[-67])) # Drop the last sample from different platform

meta = q@gpls$GPL271@dataTable@table
rownames(y) = meta$CLONE_ID

emat = backgroundCorrect.matrix(y, method = "normexp", normexp.method = 'mle', offset = 16)
emat = log2(emat)
emat = normalizeQuantiles(emat)

mapping = meta[, c("CLONE_ID", "GENE")]
mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "GENE", 
                                    probeColname = "CLONE_ID")

exprsByGene = calcExprsByGeneEmat(emat, mapping)
sum(is.na(exprsByGene))
eset <- new("ExpressionSet", exprs = exprsByGene)

saveRDS(eset, "GSE343.rds")


