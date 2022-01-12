library(tidyverse)
library(Seurat)

################################################################################################################################################
#### Add the zonation (clustering-based way) to the metadata of the seurat object
################################################################################################################################################

#### LSECs ########
output_lsecs = readRDS("output/output_liver_zonation_lsecs.rds")
output_lsecs$zonation_metadata  %>% saveRDS("output/lsec_zonation_ids.rds")
zonation_metadata_lsecs = output_lsecs$zonation_metadata 

#### Hepatocytes ########
output_hepatocytes = readRDS("output/output_liver_zonation_Hepatocytes.rds")
output_hepatocytes$zonation_metadata  %>% saveRDS("output/hepatocyte_zonation_ids.rds")
zonation_metadata_hepatocytes = output_hepatocytes$zonation_metadata

#### Stellate cells ########
zonation_metadata_stellate = readRDS("output/stellate_zonation_ids.rds")

#### Combine zonation metadata with the other metadata ##########
zonation_metadata = list(zonation_metadata_lsecs, zonation_metadata_hepatocytes,  zonation_metadata_stellate) %>% bind_rows()
zonation_metadata  %>% saveRDS("output/zonation_metadata.rds")

