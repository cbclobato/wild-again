###Tax4Fun2
#Reference: https://github.com/bwemheu/Tax4Fun2


setwd("C:/Users/Bziuk/Desktop/Test/NB02")

###Install package, reference data and dependencies
#For buildDependencies, there was a probelm with installing blastn. Details in Test Skript.
#Use use_force = T

#install.packages(pkgs = "Tax4Fun2_1.1.5.tar.gz", repos = NULL, source = TRUE)
library(Tax4Fun2)

#buildReferenceData()

#testReferenceData(path_to_reference_data = "Tax4Fun2_ReferenceData_v2")

# buildDependencies(path_to_reference_data = "./Tax4Fun2_ReferenceData_v2",
#                   install_suggested_packages = TRUE,
#                   use_force = TRUE
# )

#getExampleData(path_to_working_directory = "C:/Users/Bziuk/Desktop/Test")



###Making functional predictions
#Dna-sequences and otu table formatted as example data -> easy
#Folder: path_to_temp_folder -> made by me
#All files created by RefBlast and MakeFunctPred are saved automatically in this folder

## For Ref99NR
# 1. Run the reference blast
runRefBlast(
  path_to_otus = "dna-sequences.fasta",
  path_to_reference_data = "Tax4Fun2_ReferenceData_v2",
  path_to_temp_folder = "NB02_Ref99NR",
  database_mode = "Ref99NR",
  use_force = T,
  num_threads = 6
)

# 2. Predicting functional profiles
makeFunctionalPrediction(
  path_to_otu_table = "OTUtable.txt",
  path_to_reference_data = "Tax4Fun2_ReferenceData_v2",
  path_to_temp_folder = "NB02_Ref99NR",
  database_mode = "Ref99NR",
  normalize_by_copy_number = TRUE,
  min_identity_to_reference = 0.97,
  normalize_pathways = FALSE
)

# or:
makeFunctionalPrediction(
  path_to_otu_table = "OTUtable.txt",
  path_to_reference_data = "Tax4Fun2_ReferenceData_v2",
  path_to_temp_folder = "NB02_Ref99NR",
  database_mode = "Ref99NR",
  normalize_by_copy_number = TRUE,
  min_identity_to_reference = 0.97,
  normalize_pathways = TRUE
)


## For Ref100NR
# 1. Run the reference blast
runRefBlast(
  path_to_otus = "dna-sequences.fasta",
  path_to_reference_data = "Tax4Fun2_ReferenceData_v2",
  path_to_temp_folder = "NB02_Ref100NR",
  database_mode = "Ref100NR",
  use_force = T,
  num_threads = 6
)

# 2. Predicting functional profiles
makeFunctionalPrediction(
  path_to_otu_table = "OTUtable.txt",
  path_to_reference_data = "Tax4Fun2_ReferenceData_v2",
  path_to_temp_folder = "NB02_Ref100NR",
  database_mode = "Ref100NR",
  normalize_by_copy_number = TRUE,
  min_identity_to_reference = 0.97,
  normalize_pathways = FALSE
)

# or:
fp2 <- makeFunctionalPrediction(
  path_to_otu_table = "OTUtable.txt",
  path_to_reference_data = "Tax4Fun2_ReferenceData_v2",
  path_to_temp_folder = "NB02_Ref100NR",
  database_mode = "Ref100NR",
  normalize_by_copy_number = TRUE,
  min_identity_to_reference = 0.97,
  normalize_pathways = TRUE
)



###Calculating functional redundancy indices
# 1. Run the reference blast
# redundant, as it is the same blast as above
runRefBlast(
  path_to_otus = "dna-sequences.fasta",
  path_to_reference_data = "Tax4Fun2_ReferenceData_v2",
  path_to_temp_folder = "NB02_Ref99NR",
  database_mode = "Ref99NR",
  use_force = T,
  num_threads = 6
)

# 2. Calculating FRIs
calculateFunctionalRedundancy(
  path_to_otu_table = "OTUtable.txt",
  path_to_reference_data = "Tax4Fun2_ReferenceData_v2",
  path_to_temp_folder = "NB02_Ref99NR",
  database_mode = "Ref99NR",
  min_identity_to_reference = 0.97
)

#New in the latest pre-release (v1.1.6): prevalence_cutoff (see comment on pre-release)
calculateFunctionalRedundancy(
  path_to_otu_table = "OTUtable.txt",
  path_to_reference_data = "Tax4Fun2_ReferenceData_v2",
  path_to_temp_folder = "NB02_Ref99NR",
  database_mode = "Ref99NR",
  min_identity_to_reference = 0.97,
  prevalence_cutoff = 1.0
)
