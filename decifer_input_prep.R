#loading all required libraries
library(readxl)
library(data.table)
library(dplyr)
library(stringr)

#load the exome seq
Twist_All_reads <- read_excel("Patient1_Twist_All_reads.xlsx")
#load the copy numbers
CNA_Input_Devolution_cleaned <- read_excel("Patient1_CNA_Input Devolution_cleaned.xlsx")
#load the purity
Patient1_TCFs <- read_excel("Patient1_decifer_purity.xlsx")

# long table conversion of exome data
long <- melt(setDT(Twist_All_reads[,-3]), id.vars = c("chr","pos","ref","alt"), variable.name = "sampleID")

# splitting ref/alt reads to columns
long <- long %>%
  mutate(
    variantID = paste0(chr,":",pos,"_",ref,alt),
    alt = str_split_fixed(value, "/", n = 2)[, 1],  # Extract first element
    ref = str_split_fixed(value, "/", n = 2)[, 2],   # Extract second element
    ref=as.numeric(ref)-as.numeric(alt)
  ) %>%
  select(-c(value))

#declare background ploidy
if('2n' %in% unique(CNA_Input_Devolution_cleaned$Ploidi)){
  base.ploidy<-"1+1" } else if('3n' %in% unique(CNA_Input_Devolution_cleaned$Ploidi)){
  base.ploidy<-"2+1"} else if('4n' %in% unique(CNA_Input_Devolution_cleaned$Ploidi)){
  base.ploidy<-"2+2" } else if('5n' %in% unique(CNA_Input_Devolution_cleaned$Ploidi)){
  base.ploidy<-"3+2" } else if('6n' %in% unique(CNA_Input_Devolution_cleaned$Ploidi)){
  base.ploidy<-"3+3" } else {
  base.ploidy<-NA
}

# create a clean copy of the copy number data
cna <- CNA_Input_Devolution_cleaned %>% 
  filter(Type!="amp") %>% 
  mutate(Clone.size = as.numeric(Clone.size))

# Define function to check conditions within each row of data
create.cna.table <- function(cna = cna) { 
  unique_sample_ids <- unique(cna$Sample.ID)
  unique_Cytoband <- unique(cna$Cytoband)
  cn.test<-tibble()
  
  for (cn in unique_Cytoband) {
    dt <- cna %>% filter(Cytoband == cn)
    cn.states <- unique(c(dt$Type, base.ploidy))
    cn.store<-tibble()
  
      for (sample_id in unique_sample_ids) {
        # Filter rows for the current sample_id
        sample_rows <- dt[dt$Sample.ID == sample_id, ]
        
        # Check if all possible values in cn.states exist
        if (!all(cn.states %in% sample_rows$Type)) {
          # Create new rows for missing Type values
          missing_values <- setdiff(cn.states, sample_rows$Type)
    
            for (missing_value in missing_values) {
            # Copy an existing row and modify the Type value
            new_row <- dt[1, ]  # Choose an existing row (e.g., the first row)
            new_row$Type <- missing_value
            new_row$Sample.ID <- sample_id
            new_row$Clone.size <- 0  # Set Clone.size to 0
            sample_rows <- rbind(sample_rows, new_row)  # Append the new row to the dataframe
           }
        } 
        sample_rows$Clone.size[nrow(sample_rows)] = 100 - sum(sample_rows$Clone.size[-nrow(sample_rows)])
        cn.store <- rbind(cn.store,sample_rows)
      }
    cn.test <- rbind(cn.test,cn.store)
  }
return(cn.test)
}
cna<-create.cna.table(cna)

# Append base ploidy for variants without CNA and create CN column for DeCiFer 

apply_logic <- function(row,cna=cna,base.ploidy=base.ploidy) {
  filtered_cna <- cna %>%
    filter(Sample.ID == row$sampleID,
           Chr == row$chr,
           as.numeric(Start) < as.numeric(row$pos),
           as.numeric(End)   > as.numeric(row$pos))
  nr.clone=nrow(filtered_cna)
  if (sum(as.numeric(filtered_cna$Clone.size)) == 0) {
    base <-paste(base.ploidy, (100-sum(as.numeric(filtered_cna$Clone.size)))/100,
                 sep = "\t", collapse = "\t")
          cn <- paste(base,filtered_cna$Type,as.numeric(filtered_cna$Clone.size)/100,
                     sep = "\t", collapse = "\t")
          state <- paste(base.ploidy,filtered_cna$Type,sep = "\t", collapse = "\t")} else{
                       cn <- paste(filtered_cna$Type,
                        as.numeric(filtered_cna$Clone.size)/100,
                        sep = "\t", collapse = "\t")
                       state <- paste(filtered_cna$Type,sep = "\t", collapse = "\t")
                     }
  return(list(cn=cn,state=state))
}


# Create sample and mutation indices for DeCiFer 
long$sample_index = as.numeric(factor(long$sampleID))-1
long$Mutation_index = as.numeric(factor(long$variantID))+100

#Storing the sample name and index mapping for later use with purity
sample_map <- long %>% select(Sample = sampleID,sample_index) %>%  distinct()

var <- long %>%
    rowwise() %>%
    mutate(value = apply_logic(pick(everything()),cna,base.ploidy)$cn,
           state = apply_logic(pick(everything()),cna,base.ploidy)$state,
           CNinfo = gsub("\\+", "\t", value),
           states = gsub("\\+", ",", state),
           states = gsub("\t", ";", states),
           states = gsub("^;+|;+$", "", states)) %>%
  arrange(chr,pos,sampleID) %>% 
  select('#sample_index'=	sample_index,
         sample_label= sampleID,
         character_index=	Mutation_index,
         character_label=	variantID,
         ref=	ref,
         var= alt,
         states,
         'CN' = CNinfo) 

# Spread the CN info in tab separated columns
max_values <- max(sapply(strsplit(var$CN, "\t"), length))
new_columns <- paste0("col", 1:max_values)
var[new_columns] <- t(sapply(strsplit(var$CN, "\t"), function(x) {
  c(x, rep("", max_values - length(x)))
}))

#Extract all the CNA states to build state trees if needed
CNstates <- var %>% select(states) %>% distinct()

var <- var %>%
  select(-CN,-states)

purity  <- Patient1_TCFs %>% 
  full_join(sample_map, by='Sample') %>% 
  arrange(sample_index) %>% 
  select('#INPUT' = sample_index,
         PURITY=contains('TCF')) 

# Extracting DeCiFer main input file
write.table(var,"Patient1_decifer_input.tsv",sep="\t",row.names=F, col.names=F, quote=F)

# Extracting DeCiFer input file for purity
write.table(purity,"Patient1_decifer_purity.tsv",sep="\t",row.names=F, col.names=F, quote=F)

# Extracting DeCiFer CNA states for state tree
write.table(CNstates,"Patient1_decifer_CNAstates.tsv",sep="\t",row.names=F, col.names=F, quote=F)


