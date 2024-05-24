#!/usr/bin/env Rscript

# Set up ------------------------------------------------------------------

rm(list = ls())

#This script could be easily written in base R, but I am used to writing code using these other packages.
library(dplyr) #For %>% pipes.
library(stringr) #For str_sub function.
library(readr) #For write_lines function.

# Read command line arguments ---------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

#The locations of the MaxEntScan scripts.
mes_5ss <- args[1]
mes_3ss <- args[2]

#When you run MaxEntScan, you need to have all the splice site sequences saved as a file.
seq_dir <- args[3]

# target sequence
sequence_5ss <- args[4]  # viral seq
sequence_3ss <- args[4]  # viral seq

# me2x5 file path
me2x5 <- args[5]

# splice5seqs file path
splice5seqs <- args[6]

# splicemodels dir path
splicemodels <- args[7]

# Functions ---------------------------------------------------------------

#Basic functions for taking a long sequence and getting the 3ss or 5ss scores.

get_5ss_from_seq <- function(sequence){
  
  #Chop the sequence up in to 9-mers at each position.
  short_seqs <- str_sub(sequence, 
                        c(1:(nchar(sequence) - 8)),
                        c(9:nchar(sequence)))
  
  
  #Save the 9-mers.
  write_lines(short_seqs, file = paste0(seq_dir, "_ss5"))
  
  #Use the system command to execute the MaxEntScan Perl script.
  ss5 <- system(paste0("perl ", mes_5ss, " ", seq_dir, "_ss5 ", me2x5, " ", splice5seqs),
                intern = T) %>% 
    
    #The output from the script comes directly to R. It contains the sequence and the score, but we only want the score.
    #The score starts at position 11 in the output string, so I take everything from position 11 onwards to the end (somewhere before position 20).
    str_sub(11,20) %>% 
    
    #The splice site scores are still classed as a string instead of a number, so I convert them.
    as.numeric()
  
  return(ss5)

}

get_3ss_from_seq <- function(sequence){
  
  #Chop the sequence up in to 23-mers at each position.
  short_seqs <- str_sub(sequence, 
                        c(1:(nchar(sequence) - 22)),
                        c(23:nchar(sequence)))
  
  
  #Save the 23-mers.
  write_lines(short_seqs, file = paste0(seq_dir, "_ss3"))
  
  #Use the system command to execute the MaxEntScan Perl script.
  ss3 <- system(paste0("perl ", mes_3ss, " ", seq_dir, "_ss3 ", splicemodels),
                intern = T) %>% 
    
    #The output from the script comes directly to R. It contains the sequence and the score, but we only want the score.
    #The score starts at position 11 in the output string, so I take everything from position 25 onwards to the end.
    str_sub(25,50) %>% 
    
    #The splice site scores are still classed as a string instead of a number, so I convert them.
    as.numeric()
  
  return(ss3)
  
}

# Running the functions ---------------------------------------------------

#Take a sequence.

#Get splice site scores for the positions along the sequence.
ss5_scores <- get_5ss_from_seq(sequence_5ss)

df5 <- data.frame(ss5_scores)
df5$location2<-seq_along(df5[,1])

ss3_scores <- get_3ss_from_seq(sequence_3ss)

df3 <- data.frame(ss3_scores)
df3$location<-seq_along(df3[,1])
df3$location2<-(df3$location + 17)

# Output ------------------------------------------------------------------
write_delim(df5, file = "5ss_scores.tsv", delim = "\t")
write_delim(df3, file = "3ss_scores.tsv", delim = "\t")