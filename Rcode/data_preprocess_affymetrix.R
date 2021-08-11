BiocManager::install(c("readxl", "affy", "GEOquery", "tidyverse", "foreach", "doParallel", "SCAN.UPC"))

library("GEOquery")
library("tidyverse")
library("foreach")
library("doParallel")
library("SCAN.UPC") # read vignette

# Custom packages for probe-gene mapping:
# http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/25.0.0/entrezg.asp
# Read vignette from SCAN.UPC


fixCustomCdfGeneIds = function(geneIds) {
    return(sub('_at', '', geneIds))
}

# Loading metadata ----
metadata = readxl::read_xlsx("../Metadata/Biopsy_Tx.xlsx", sheet = 1) %>% as.data.frame()
studies = metadata[metadata$Platform_manufacturer == "Affymetrix", "Study", drop = T]

# Downloading samples from GEO ----
# Use it only if no data is not downloaded already
# Downloading only relevant to current analysis samples
for (name in studies) {
    message("Downloading study ", name)
    samples2down = metadata[metadata$study == name, "sample"]
    dir.create(paste0("../Data_in/", name))
    foreach (sample = samples2down, .packages = "GEOquery") %dopar% {
        tmp = getGEOSuppFiles(GEO = sample,
                              makeDirectory = FALSE,
                              baseDir = paste0("../Data_in/", name),
                              fetch_files = TRUE,
                              filter_regex = NULL)
    }
}

# File name fix ----
datapath = '../Data_in'
dirz = list.dirs(datapath)[-1]
for (onedir in dirz) {
    cwd = setwd(onedir)
    filez = dir() #list.celfiles(onedir)
    for (onefile in filez) {
        parts = strsplit(onefile, split = "[._]")[[1]]
        newfile = paste0(c(parts[1], tail(parts, 2)), collapse = ".")
        file.rename(onefile, newfile)
    }
    setwd(cwd)
}

# Registering the multithreading
cl <-  makeCluster(3, type = 'PSOCK', outfile = "")
doParallel::registerDoParallel(cl)

# Data processing with SCAN.UPC -----
for (studyName in studies) {
    pkgName = NA
    celFilePath = file.path("../Data_in", studyName, "*.CEL.gz")
    
    message("Working on ", studyName, "...")
    files = affy::list.celfiles(file.path("../Data_in/", studyName))
    
    # Reading platform info from a CEL file
    platform = cleancdfname(read.celfile.header(file.path("../Data_in/", studyName, files[1]), info = "full")$cdfName)
    platform = sub("cdf", "", platform)
    platform = sub("stv1|stv2", "st", platform)
    pkgName = paste(platform, "hsentrezgprobe", sep = "")
    
    # Attempting to load Brain Array custom package for mapping probes to genes
    # If a package is not installed, download and install it from Brain Array
    if (!require(pkgName, character.only = TRUE)) {
        pkgName = InstallBrainArrayPackage(file.path("../Data_in/", studyName, files[1]), "25.0.0", "hs", "entrezg") # read vignette
        require(pkgName, character.only = TRUE)
    }
    
    # Process the dataset
    outFilePath = file.path("../Data_out/SCAN_txt", paste0(studyName, "_SCAN.txt"))
    normSet = SCAN(celFilePath, probeSummaryPackage = pkgName, outFilePath, verbose = TRUE)
    
    # Removing '_at' from gene ids
    message("Fixing Gene Ids in ", studyName)
    rownames(normSet) = fixCustomCdfGeneIds(rownames(normSet))
    
    # message("Fixing sample names in ", studyName)
    # colnames(normSet) = fixCelSampleNames(colnames(normSet))
    
    saveRDS(normSet, file = file.path("../Data_out", paste0(studyName, ".rds")))
}

stopCluster(cl)

