## Original code acquired from Chris Murray 
## Modified by Myung Chang Lee (Noah Lee) on April 8th, 2019
## to extend functionality to multiple files using functions
## and clean up formula


##requires dplyr/tidyr
library(dplyr)
## I'd rather output to Excel than csv to prevent dataloss
library(openxlsx)

setwd("C:/Users/noahl/Box/Sage Lab Server/Noah/sgRNA")

rm(list = ls()) # Clear workspace

options(stringsAsFactors = FALSE) 


input_dir <- "./Designs/"
output_dir <- "./Filtered/"


##Enter restriction sites that are used for barcoding
BamHI <- 'GGATCC'
BspEI <- 'TCCGGA'
XmaI <- 'CCCGGG'
RE_Sites <- c(BamHI, BspEI, XmaI)


# Get a list of files in the directory
files <- list.files(path = input_dir)

## If you want, you can simply have files be a list of files that you want program to read
# files <- c("mKDR sgRNA All.csv")



# Function to filter sgRNA ------------------------------------------------

filtersgRNA <- function(sgRNA_Scoring, offtarget = 0.6, percentpepmax = 0.7, ontarget = 0.6){
  ##Remove sgRNAs containing restriction sites
  Filtered_sgRNAs <- sgRNA_Scoring %>%
    filter(!grepl(paste(RE_Sites, collapse = '|'),bases))
  
  ##Set minimimum threshold for off-target scoring
  Filtered_sgRNAs <- Filtered_sgRNAs %>%
    filter(scores.offtarget.hsu_2013.value > offtarget)
  
  ##Subset guides with maximum representation across transcript isoforms, representation values range between 0.0 and 1.0
  Filtered_sgRNAs <- Filtered_sgRNAs %>%
    filter(scores.site.representation.value == max(scores.site.representation.value))
  
  ##Filter guides based on their position within the CDS, values range from 0.0 to 1.0, select sgRNAs that are upstream of major features but avoid extreme ends of CDS
  Filtered_sgRNAs <- Filtered_sgRNAs %>%
    filter(scores.site.percent_peptide.value > 0.05) %>%
    filter(scores.site.percent_peptide.value < percentpepmax)
  
  if(nrow(Filtered_sgRNAs) < 2){
    print("No more sgRNA candidates remaining after % pep filter")
    return(Filtered_sgRNAs)
  }
  
  ##Filter out sgRNAs containing homopolymers
  Filtered_sgRNAs <- Filtered_sgRNAs %>%
    filter(scores.site.no_homopolymer.value == 1)
  
  ##Filter out 'UUU' to avoid transcription termination when expressing from U6
  Filtered_sgRNAs <- Filtered_sgRNAs %>%
    filter(scores.site.no_uuu.value == 1)
  
  ##Filter GC content to maximize activity
  Filtered_sgRNAs <- Filtered_sgRNAs %>%
    filter(0.4 < scores.site.gc_content.value & scores.site.gc_content.value < 0.75)
  
  if(nrow(Filtered_sgRNAs) < 2){
    print("No more sgRNA candidates remaining after GC filter")
    return(Filtered_sgRNAs)
  }
  
  ##Set minimum threshold for on-target activity with Doench 2016 guidelines
  Filtered_sgRNAs <- Filtered_sgRNAs %>%
    filter(scores.activity.doench_2016_full.value > ontarget)
  
  if(nrow(Filtered_sgRNAs) < 2){
    print("No more sgRNA candidates remaining after on-target filter")
    return(Filtered_sgRNAs)
  }
  
  ##Alternative scoring algorithm if filtering on Doench 2016 guidelines yield little to no candidate sgRNAs
  #Filtered_sgRNAs <- Filtered_sgRNAs %>%
  #filter(scores.activity.chari_2015.value > 0.6)
  
  ##Set threshold for predicted microhomology score to increase likelihood of introducing frameshift (can shift down to 0.5 minimum if needed)
  Filtered_sgRNAs <- Filtered_sgRNAs %>%
    filter(scores.activity.bae_2014.value > 0.6)
  
  if(nrow(Filtered_sgRNAs) < 2){
    print("No more sgRNA candidates remaining after microhomology filter")
    return(Filtered_sgRNAs)
  }
  
  ##Create aggregate score for on-target activity, off-target activity and microhomology (max score of 4)
  Filtered_sgRNAs <- Filtered_sgRNAs %>%
    mutate(Aggregate_Score = scores.offtarget.hsu_2013.value + scores.activity.doench_2016_full.value*2 + scores.activity.bae_2014.value)
  
  ##Order sgRNAs by aggregate score (decreasing)
  Filtered_sgRNAs <- Filtered_sgRNAs[order(-Filtered_sgRNAs$Aggregate_Score),]
  
  ##View
  head(Filtered_sgRNAs)
  
  return(Filtered_sgRNAs)
}

# Function to get filtered sgRNA in a given file --------------------------------

getsgRNA <- function(filename){
  ##Import downloaded .csv from Desktop Genetics, make sure to indicate that there is a heading
  sgRNA_Scoring <- read.csv(paste0(input_dir, filename))
  
  Filtered_sgRNAs <- filtersgRNA(sgRNA_Scoring)
  
  offtarget = 0.6
  percentpepmax = 0.7
  ontarget = 0.6
  pass = 1
  
  while(nrow(Filtered_sgRNAs) < 2){
    if(pass == 1){
      print(paste(filename,"yielded no sgRNAs"))
      print("Lowering standards...")
    }
    if(pass %% 2 == 1){
      ontarget <- ontarget - 0.02
    }
    else {
      offtarget <- offtarget - 0.02
    }
    Filtered_sgRNAs <- filtersgRNA(sgRNA_Scoring, offtarget, percentpepmax, ontarget)
    pass <- pass + 1
    
    if(pass > 5){
      print("Max iterations reached. Breaking")
      break
    }
  }
  
  if(pass > 1){
    print(paste0("Obtained sgRNA on Pass #", pass))
    print(paste("Ontarget:",ontarget,"; offtarget:", offtarget))
  }
  
  ##Export to .xlsx, sgRNA seed sequences are located in 'bases' column
  write.xlsx(Filtered_sgRNAs, paste0(output_dir, gsub(" sgRNA.+\\.csv", "\\1", filename), '_filtered_sgRNAs.xlsx'))
}

# Go through files in the input directory and filter ----------------------

for(file in files){
  getsgRNA(file)
}

