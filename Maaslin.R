## MAASLIN ANALYSIS
# https://huttenhower.sph.harvard.edu/maaslin/
# NOVEMBER 2020
rm(list = ls(all = TRUE))
gc()
# MaAsLin2 is used for determining multivariable association between phenotypes, environments, exposures, covariates and microbial meta’omic features. 
# it relies on general linear models to accommodate most modern epidemiological study designs, including cross-sectional and longitudinal, and offers a variety of data exploration, normalization, and transformation methods.

directory <- "/home/eugloh/Documents/DAVID/Maaslin2Test"
setwd(directory)
# check for installation and packages 
source("EnvironmentSetup.R")
library("Maaslin2")
option_list = list(
  make_option(c("-p", "--pattern"), type="character", default="tsv", help="pattern for file names [default= %default]", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

tsvFile  <- grep(opt$pattern,list.files(paste0(directory,"/tsv")),value=TRUE) # find names
metaFile  <- grep("Meta",list.files(paste0(directory,"/Metadata")),value=TRUE) # find names

DF_DATA<- list()
for (el in tsvFile){
  el_name <- sapply(strsplit(as.character(el), "\\."), "[[", 1)
  DF_DATA[[el_name]] <- read.delim(paste(directory,"tsv",el, sep = "/"), row.names=1)
}

DF_METADATA<- list()
for (el in metaFile){
  el_name <- sapply(strsplit(as.character(el), "\\."), "[[", 1)
  DF_METADATA[[el_name]] <- read.delim(paste(directory,"Metadata",el, sep = "/"), row.names=1)
}

names(DF_METADATA)
# [1] "Metadataall"   "MetadataF_all" "MetadataF"     "MetadataFOB"   "MetadataS_all" "MetadataS"    
# [7] "MetadataSOB"  
names(DF_DATA)
# [1] "ec_metagenome_F"       "ec_metagenome_FOB"     "ec_metagenome_S"       "ec_metagenome_SOB"    
# [5] "ko_metagenome_F"       "ko_metagenome_FOB"     "ko_metagenome_S"       "ko_metagenome_SOB"    
# [9] "pathway_abundance_F"   "pathway_abundance_FOB" "pathway_abundance_S"   "pathway_abundance_SOB"
# [13] "pathway_rel-freq_F"    "pathway_rel-freq_FOB"  "pathway_rel-freq_S"    "pathway_rel-freq_SOB"

## Running MaAsLin2
### test
for (el in names(DF_DATA)){
colnames(DF_DATA[[el]]) <- str_replace_all(colnames(DF_DATA[[el]]), pattern=fixed("."), replacement="-")
}

colnames(DF_METADATA[["MetadataFOB"]]) <-c("Run","Studies","Sample", "Age","PhyYEScalActivity","BMI","R.BMI.OverTime","Gender","Smoke","Weight",
                                           "MetS","WHratio","InterventionType","PCR","TNF.","IL6","BarcodeSequence")

colnames(DF_METADATA[["MetadataSOB"]]) <- c("Run","Studies","Age" ,"PhysicalActivity" ,"BMI","R.BMI.OverTime" ,             
"Gender"  , "Smoke","Weight","MetS","MetS.No.Yes.","WHratio","InterventionType","PCR" ,                        
"TNFα","IL6", "MultipleInfection..beta.HPV.","PersistentInfection","Sample" , "BarcodeSequence" )            

colnames(DF_METADATA[["MetadataF"]]) <- c("Run","Sample","TimePoint","Gender","TwinType",       
                                          "Smoke","BMI","BarcodeSequence")
colnames(DF_METADATA[["MetadataS"]]) <- c("Run","TimePoint","Sample","Gender",       
                                          "Smoke","TwinType","BarcodeSequence")


for (statut in c("F","S","FOB","SOB")){
  for (file_type in c("ec_metagenome","ko_metagenome","pathway_rel-freq")){
    Maaslin2(
      input_data = DF_DATA[[eval(paste0(file_type,"_",statut))]], 
      input_metadata = DF_METADATA[[eval(paste0("Metadata",statut))]], 
      output = eval(paste0("Studies","/",file_type,"/",file_type,"_",statut)),
      fixed_effects = c("TimePoint"))
  }
}



for (statut in c("S","F","FOB","SOB")){
  for (file_type in c("ec_metagenome","ko_metagenome","pathway_rel-freq")){
    Maaslin2(
      input_data = DF_DATA[[eval(paste0(file_type,"_",statut))]], 
      input_metadata = DF_METADATA[[eval(paste0("Metadata",statut))]], 
      output = eval(paste0("Gender","/",file_type,"/",file_type,"_",statut)),
      random_effects = c("Gender"))
  }
}

colnames(DF_METADATA[["MetadataFOB"]]) <-c("Run","Studies","Sample", "Age","PhyYEScalActivity","BMI","R.BMI.OverTime","Gender","Smoke","Weight",
"MetS","WHratio","InterventionType","PCR","TNF.","IL6","BarcodeSequence")

for (statut in c("S","F","FOB","SOB")){
  for (file_type in c("ec_metagenome","ko_metagenome","pathway_rel-freq")){
    Maaslin2(
      input_data = DF_DATA[[eval(paste0(file_type,"_",statut))]], 
      input_metadata = DF_METADATA[[eval(paste0("Metadata",statut))]], 
      output = eval(paste0("add_random","/",file_type,"/",file_type,"_",statut)),
      random_effects = c("Run"),
      fixed_effects = c("Studies"))
  }
}


