####################################################################################################################################################################
# Creator: Ilian Torres
# Date: August 3,3016
# Works for Version 0.11.5 
####################################################################################################################################################################
#Libraries
library(plyr)
library(ggplot2)
library(tcltk)
library(asbio)
library(nortest)
library(tidyr)
####################################################################################################################################################################
#   Functions   ####################################################################################################################################################
####################################################################################################################################################################
#Calculate_Confidence_Intervals_of_Table calculate confidence intervals using the median for Per sequence quality scores and
#returns a list containing the samples that dont fall inside the confidence intervals, Median, Upper and Lower confidence intervals.
#As input it takes the quality score values (x_values) and Per sequence quality scores Merged into one dataframe (Merged_File) 
#that consisnt of 3 rows:Quality,Count,Name (of file). 
Calculate_Confidence_Intervals_of_Table <- function(Merge_File,x_value,npar=TRUE){
  colnames(Merge_File)<-c("point","Name","value")
  Bad_data=NULL
  Bad_data_Names=NULL
  BG=NULL
  Median=NULL
  Lower=NULL
  Upper=NULL
  Data=NULL
  Median_table=NULL
  Upper_table=NULL
  Lower_table=NULL
  Names=NULL
  # For loop that loops through all the quality scores (x).
  # 1) Quality_Subset subsets all the Quality scores that is equal to x.
  # 2) Con_Int calculates the confidence intervals of Quality_Subset.
  # 3) Then the values are stored in thier corresponding dataframes
  # 4) Names Store the Names in which the data does not fall in between the confident intervals.
  for (x in x_value){
    Quality_Subset=subset(Merge_File,Merge_File$point==x)
    Con_Int=ci.median(Quality_Subset$value,.975)
    Con_Int=setNames(data.frame(t(Con_Int$ci)),Con_Int$ends)
    Median$point=x
    Median$value=Con_Int[[1]]
    Lower$point=x
    Lower$value=Con_Int[[2]]
    Upper$point=x
    Upper$value=Con_Int[[3]]
    Median=data.frame(Median)
    Lower=data.frame(Lower)
    Upper=data.frame(Upper)
    Median_table=rbind(Median_table,Median)
    Lower_table=rbind(Lower_table,Lower)
    Upper_table=rbind(Upper_table,Upper)
    Names = ifelse(Quality_Subset$value <= Lower$value[1] | Quality_Subset$value >= Upper$value[1],Quality_Subset$Name,0)
    Names = subset(Names,Names!=0)
    Bad_data_Names=c(Bad_data,Names)
  }
  Bad_data_Names=unique(Bad_data_Names)
  #For loop that loops through the names and subsets the data in Bad_data.
  for ( i in Bad_data_Names ){
    BG=subset(Merge_File,Merge_File$Name==i)
    Bad_data=rbind(Bad_data,BG)
  }
  #Returns a list 
  Data=list(Median_table,Lower_table,Upper_table,Bad_data)
  return(Data)
}

####################################################################################################################################################################
# Verifying_Intervals takes the data given by Calculate_Confidence_Intervals_of_Table and eliminates
# all the data from Quality 35-40 that are above the confindence interval but that Quality 2-34 lies under the upper confidence interval.
# As input it takes the data created by Calculate_Confidence_Intervals_of_Table.
Verifying_Intervals <- function(Bad_samples,npar=TRUE){
  colnames(Bad_samples)<-c("Quality","Name","Count","Lower","Upper")
  Verifying_Bad_Sample=NULL
  Subset_of_35=subset(Bad_samples,Bad_samples$Quality>=35)
  Backup=Bad_samples
  # Takes Subset_of_35, which is a subset all the Quality above 35, and subsets all the data that falls above the confedence interval.
  # Also all the data that above the confidence interval it eliminates it from Bad_samples.
  for (i in 1:nrow(Subset_of_35)){
    if (Subset_of_35$Count[i] >= Subset_of_35$Upper[i]){
      Subset_Questionable_sample=subset(Bad_samples,Bad_samples$Name==Subset_of_35$Name[i])
      Bad_samples=Bad_samples[Bad_samples$Name != Subset_of_35$Name[i],]
      Subset_Questionable_sample=data.frame(Subset_Questionable_sample)
      Verifying_Bad_Sample=rbind(Verifying_Bad_Sample,Subset_Questionable_sample)
    }
  }
  Subset_of_35_Under=subset(Verifying_Bad_Sample,Verifying_Bad_Sample$Quality<35)
  # Takes Subset_of_35_under, which is a subset all the Quality below 35 of Verifying_Bad_Sample, which is all the data thats is above 35, 
  # goes line by line if the data is above the upper confidence interval it places it back in Bad_samples 
  for (i in 1:nrow(Subset_of_35_Under)){
    if (is.na(Subset_of_35_Under$Count[i])){
      break
    }
    if (Subset_of_35_Under$Count[i] >= Subset_of_35_Under$Upper[i]){
      Bad=subset(Backup,Backup$Name==Subset_of_35_Under$Name[i])
      Subset_of_35_Under=Subset_of_35_Under[Subset_of_35_Under$Name != Subset_of_35_Under$Name[i],]
      Bad=data.frame(Bad)
      Bad_samples=rbind(Bad_samples,Bad)
    }
  }
  #returns Bad_Samples
  return(Bad_samples)
}
####################################################################################################################################################################
# Finding_Continuous_Line takes as input the name of all the files and data that is composed of
#"Base","G","A","T","C","Name","G_C","A_T","Legend". The purpose of this function is to find when AT and GC have the lonegest continous line.
Finding_Continuous_Line<- function(PBSC_Subset,Name_of_files){
  colnames(PBSC_Subset)<-c("Base","G","A","T","C","Name","G_C","A_T","Legend")
  Longest_Lines=NULL
  # For loop that takes each name of the batch and subsets the data by name of file. TF_List finds where the base is 
  # NA and displays true where it is NA and false whenit is a value. Then Index_Trues takes the index of where TF_List is true.
  for (x in Name_of_files){
    Straightlines=NULL
    Verify_Mean=NULL
    Subset_Name=subset(PBSC_Subset,PBSC_Subset$Name==x)
    TF_List=is.na(Subset_Name$Base)
    TF_List=c(TF_List,TRUE)
    Index_Trues=which(TF_List %in% TRUE)
    # For loop that goes through Index_Trues and finds all the staight lines and puts them in a table.
    for ( i in 1:(length(Index_Trues)-1)){
      if (Index_Trues[i+1]-Index_Trues[i]-1>1){
        Straight=matrix(c(Index_Trues[i]+1,Index_Trues[i+1]-1,Index_Trues[i+1]-Index_Trues[i]-1),byrow = TRUE,ncol = 3)
        Straight=data.frame(Straight)
        Straightlines=rbind(Straightlines,Straight)
      }
    }
    # if there are staright lines goes in loop.
    if (!is.null(Straightlines)){
      Straightlines=data.frame(Straightlines)
      bool=TRUE
      # While bool= true means that the longest line has not been found. Leaves loop when finds
      # longest line and verifies the the line stays stright by verifying it with its mean.
      while(bool){
        # Finds longest line
        Longest_line=Straightlines[(which.max( Straightlines[,3])),]
        colnames(Longest_line)<-c("L","H","D")
        # Takes the data of the longest lines.
        Verify_Mean=Subset_Name[Subset_Name$Base>=Longest_line$L[1] & Subset_Name$Base<=Longest_line$H[1],]
        Verify_Mean=na.omit(Verify_Mean)
        rownames(Verify_Mean)<-c(1:nrow(Verify_Mean))
        Verify_Mean=data.frame(Verify_Mean)
        MeanAT=mean(Verify_Mean$A)
        MeanGC=mean(Verify_Mean$G)
        # Verify that the longest line falls between the mean+-2
        for (i in 1:nrow(Verify_Mean)){
          if (MeanAT+2>=Verify_Mean$A[i] | MeanAT-2<=Verify_Mean$A[i] | MeanAT+2>=Verify_Mean$T[i] | MeanAT-2<=Verify_Mean$T[i] | 
              MeanGC+2>=Verify_Mean$G[i] | MeanGC-2<=Verify_Mean$G[i] | MeanGC+2>=Verify_Mean$C[i] | MeanGC-2<=Verify_Mean$C[i]){
            bool=TRUE
            break
          }
          else
            bool=FALSE
        }
        if (bool){
          Longest_Lines=rbind(Longest_Lines,Verify_Mean)
          bool=FALSE
        }
      }
    }
    # Doesent have a staight line fills it with NA.
    else{
      NA_Data=c(NA,NA,NA,NA,NA,Subset_Name$Name[i],NA,NA,NA)
      Longest_Lines=rbind(Longest_Lines,NA_Data)
    }
  }
  Longest_Lines$Base <- as.numeric(Longest_Lines$Base)
  Longest_Lines$A <- as.numeric(Longest_Lines$A)
  Longest_Lines$T <- as.numeric(Longest_Lines$T)
  Longest_Lines$G <- as.numeric(Longest_Lines$G)
  Longest_Lines$C <- as.numeric(Longest_Lines$C)
  return(Longest_Lines)
}
####################################################################################################################################################################
# Function that slits data for more understandible graphs. Input the data, name of the files that are in the data, and what number of data per graph(div).
Splitting_data <- function (Good_data,Name_of_files,div,npar=TRUE){
  len = length(Name_of_files)
  rest=len %% div
  value = len-rest
  Datalist=NULL
  value_div=value/div
  Data_sep_div=NULL
  t=1
  for (i in (1:(t*div))){
    name=Name_of_files[i]
    Data_subset=subset(Good_data,Good_data$Name==name)
    Data_sep_div=rbind(Data_sep_div,Data_subset)
  }
  Datalist[[t]]=Data_sep_div
  Data_sep_div=NULL
  t=t+1
  while(t!=value_div){
    for (i in ((t*div+1):(t*div+div))){
      name=Name_of_files[i]
      Data_subset=subset(Good_data,Good_data$Name==name)
      Data_sep_div=rbind(Data_sep_div,Data_subset)
    }
    Datalist[[t]]=Data_sep_div
    Data_sep_div=NULL
    t=t+1
  }
  for (i in ((t*div+1):(t*div+rest))){
    name=Name_of_files[i]
    Data_subset=subset(Good_data,Good_data$Name==name)
    Data_sep_div=rbind(Data_sep_div,Data_subset)
  }
  Datalist[[t]]=Data_sep_div
  return(Datalist)
}
####################################################################################################################################################################
#Normality test for Per Sequence GC Content
Normality_Test <- function(Merge_File,Name_of_files,npar=TRUE){
  colnames(Merge_File)<-c("point","value","Name")
  Bad_data=NULL
  for (x in Name_of_files){
    Data_Subset=subset(Merge_File,Merge_File$Name==x)
    Data=Data_Subset$value
    Shapiro=ad.test(Data)
    Shapiro=Shapiro$p.value
    if (Shapiro<.05){
      Bad_data=rbind(Bad_data,Data_Subset)
    }
  }
  return(Bad_data)
}


####################################################################################################################################################################
# Function getReadPairIDFromFq takes as input the name of the file and creates and returns his ID.
# Function Created by Holly Beale
# rules
# it's ok to have the read end indicator at the end of the name, e.g. holly_2.fq
# if it's not at the end, it can't be followed by a letter or number, and it can only be 
# e.g. holly_2a.fq  <- 2 is not considered the read end indicator
# e.g. holly_2_clean.fq  <- 2 is considered the read end indicator only if "clean" is in wordsAllowedAfterEndNumber
# e.g. holly_2_KC23072.fq  <- 2 is considered the read end indicator only if "KC23072" is in wordsAllowedAfterEndNumber
# this should fail, but doesn't:  holly_2_KC23072_2_clean"
getReadPairIDFromFq<-function(fileNames= miserableTestData, allowedSuffix=c("_fastqc.zip", "\\.fastq", "\\.fq.gz", "\\.tgz", "\\.fastq.gz", "\\.tar", "_fastqc"), wordsAllowedAfterEndNumber=c("clean", "trimmed")){
  #  extractEndsWithLowConfidencePattern=TRUE
  # remove path info if any
  
  trimPunct= function (x) gsub("^[[:punct:]]+|[[:punct:]]+$", "", x)
  #trimPunct (c("holly_", "holly__", "holly_2_KC23072_", "holly."))
  df=data.frame(originalOrder=1:length(fileNames), inputName= fileNames , baseFileName=basename(fileNames))
  
  #  strip suffix
  df$noSuffix= df$baseFileName
  for (s in allowedSuffix) df$noSuffix =sub(paste0(s, "$"), "", df$noSuffix)
  
  # if relevant, strip words that are allowed after end indicator; strip outer punctuation
  hasEndIndicatorAtStringEnd=grepl("[\\._][R]?[12]$", df $noSuffix) | grepl("[\\._][R]?[12](_[0-9]{3})$",  df $noSuffix)
  
  df$allowedWordsRemoved= df $noSuffix
  for (w in wordsAllowedAfterEndNumber) df$allowedWordsRemoved[!hasEndIndicatorAtStringEnd] =sub(w, "", df$allowedWordsRemoved[!hasEndIndicatorAtStringEnd])
  df$allowedWordsRemoved= trimPunct(df$allowedWordsRemoved)
  
  # extract read end ids
  df$readPairID=sub("[\\._][R]?[12]$", "", df$allowedWordsRemoved)
  df$rawReadEnd =sub("^.*[\\._][R]?([12])$", "\\1", df$allowedWordsRemoved)
  
  remainingToParse= df$readPairID == df$allowedWordsRemoved
  
  df$readPairID[remainingToParse ]=sub("[\\._][R]?[12](_[0-9]{3})$", "\\1",  df$allowedWordsRemoved)[remainingToParse]
  df$rawReadEnd[remainingToParse]=sub("^.*[\\._][R]?([12])_[0-9]{3}$", "\\1", df$allowedWordsRemoved[remainingToParse])
  
  df$readEnd= df$rawReadEnd
  df$readEnd[! df$rawReadEnd %in% 1:2]=NA
  
  # check for remaining read end identifiers
  df$candidateEndIndentifierRemainsInReadPairID=grepl("[\\._][R]?[12]$", df$readPairID) |grepl("[\\._][R]?[12][\\._]", df$readPairID)
  
  return(df)
  
}

####################################################################################################################################################################
# Start of Program #################################################################################################################################################
####################################################################################################################################################################
#sets working directory
setwd('Desktop/fastQC_treehouse_results')
#variables
Directories=list.dirs(recursive = FALSE)
Name_of_folder='Sep_fastqc_data'
name_of_txt = c("Adapter_Content.txt","Basic_Statistics.txt","Kmer_Content.txt","Overrepresented_sequences.txt",
                "Per_base_N_content.txt","Per_base_sequence_content.txt","Per_base_sequence_quality.txt",
                "Per_sequence_GC_content.txt","Per_sequence_quality_scores.txt","Per_tile_sequence_quality.txt",
                "Sequence_Duplication_Levels.txt", "Sequence_Length_Distribution.txt" )
Name_of_files=list.files()
Name_of_files=Name_of_files[Name_of_files != "Test.R"]
data=list()
# Makes a data list that is composed of dataframes names after their folder and contains the datatxt.
for (File in Directories){
  for (folder in File)
    data[[folder]]=list()
  for (Info in name_of_txt){
    if (Info == "Basic_Statistics.txt")
      data[[folder]][[Info]]=read.table(paste(File,"/",folder,"/",Name_of_folder,"/",Info,sep = ''),header=TRUE,comment.char = "",sep="\t",check.names = FALSE,skip = 2)
    else if (Info == "Sequence_Duplication_Levels.txt"){
      data[[folder]][['Total Duplicate Percentage']]=substring(readLines(paste(File,"/",folder,"/",Name_of_folder,"/",Info,sep = ''),n=1),"29")
      data[[folder]][[Info]]=read.table(paste(File,"/",folder,"/",Name_of_folder,"/",Info,sep = ''),header=TRUE,comment.char = "",sep="\t",check.names = FALSE,skip = 1)
    }
    else if (Info == "Overrepresented_sequences.txt"){
      txt_file=file.info(paste(File,"/",folder,"/",Name_of_folder,"/",Info,sep = ''))
      if (txt_file$size!=0)
        data[[folder]][['Overrepresented_sequences.txt']]=data[[folder]][[Info]]=
          read.table(paste(File,"/",folder,"/",Name_of_folder,"/",Info,sep = ''),header=TRUE,comment.char = "",sep="\t",check.names = FALSE)
    }
    else
      data[[folder]][[Info]]=read.table(paste(File,"/",folder,"/",Name_of_folder,"/",Info,sep = ''),header=TRUE,comment.char = "",sep="\t",check.names = FALSE)
  }
}
name_of_txt=append(name_of_txt, "Total Duplicate Percentage")
####################################################################################################################################################################
#Name of of pdf file where all graphs will be saved
pdf(file = "/Users/Ilian/Desktop/fastQC_treehouse_results/Merged.pdf", title="Graphs of Merged Data",height = 13,width = 15)
#ggplot Per_sequence_quality_scores
PSQS_Merged=NULL
gg_plot1=NULL
gg_plot2=NULL
# For loop that imports all the Data of Per_sequence_quality_scores
for (Name in Name_of_files){
  PSQS=data [[paste('./',Name,sep='')]][['Per_sequence_quality_scores.txt']]
  PSQS$Name <- Name
  PSQS_Merged <- rbind(PSQS_Merged,PSQS)
}
PSQS_Merged=rename(PSQS_Merged, c("#Quality"="Quality"))
PSQS_Merged=complete(PSQS_Merged, Quality, Name, fill = list(Count = 0))
Quality_values_PSQS=c(2:40)
Conf_Int_Table=Calculate_Confidence_Intervals_of_Table(PSQS_Merged,Quality_values_PSQS)
Lower=Conf_Int_Table[[2]]
Upper=Conf_Int_Table[[3]]
Bad_samples=Conf_Int_Table[[4]]
Bad_samples$Upper = sapply(Bad_samples$point, function(x){(Upper[Upper$point==x[1], "value"])})
Bad_samples$Lower = sapply(Bad_samples$point, function(x){(Lower[Lower$point==x[1], "value"])})
Bad_samples=data.frame(Bad_samples)
colnames(Bad_samples)<-c("Quality","Name","Count","Upper","Lower")

Bad_Data=Verifying_Intervals(Bad_samples)
Bad_Data_Names=unique(Bad_Data$Name)
Bad_data_DIV=Splitting_data(Bad_Data,Bad_Data_Names,30)


for (i in 1:(length(Bad_data_DIV))){
  Graph=Bad_data_DIV[[i]]
  PSQS_PLOT_Conf_Int <- ggplot(Graph, aes(x = Quality, y = Count,group=Name))+
    geom_line(aes(color =Name))+ggtitle("Per Sequence Quality Scores 95% Confidence Intervals Excluding Upper Outliers Above Quality 35")+
    geom_ribbon(aes(ymin=Lower,ymax=Upper,x=Quality, group=Name), alpha=0.01, linetype=1,colour="grey60", size=.05,fill="grey40")
  PSQS_PLOT_Conf_Int_IND = ggplot(Graph, aes(x = Quality, y = Count)) +
    geom_line(aes(color =Name)) +facet_wrap(~Name) + geom_ribbon(aes(ymin=Lower,ymax=Upper), alpha=0.1, linetype=1,colour="grey20", size=.05,fill="grey20")+guides(color=FALSE)
  gg_plot1[[i]]=PSQS_PLOT_Conf_Int
  gg_plot2[[i]]=PSQS_PLOT_Conf_Int_IND
}
for (i in 1:(length(gg_plot))){
  print  (gg_plot1[[i]])
  print (gg_plot2[[i]])
}

dev.off()
####################################################################################################################################################################
#ggplot per base sequence content plot
PBSC_Merged=NULL
PBSC_Sub=NULL
PBSC_Subset=NULL
g_plot=NULL
for (Name in Name_of_files){
  PBSC=data [[paste('./',Name,sep='')]][['Per_base_sequence_content.txt']]
  PBSC$Name <- Name
  PBSC_Merged <- rbind(PBSC_Merged,PBSC)
}
PBSC_Merged=rename(PBSC_Merged, c("#Base"="Base"))
PBSC_Merged=data.frame(PBSC_Merged)
PBSC_Sub=PBSC_Merged
PBSC_Sub$G_C= (abs(PBSC_Merged$G-PBSC_Merged$C))
PBSC_Sub$A_T= (abs(PBSC_Merged$A-PBSC_Merged$T))
PBSC_Sub$Legend = ifelse(PBSC_Merged$C<PBSC_Merged$T | PBSC_Merged$C<PBSC_Merged$A
                         | PBSC_Merged$G<PBSC_Merged$T | PBSC_Merged$G<PBSC_Merged$A,"AT Greater than GC","GC Greater than AT")
PBSC_Sub=data.frame(PBSC_Sub)

PBSC_Subset=PBSC_Sub
PBSC_Subset$Base=ifelse(PBSC_Sub$G_C<1 & PBSC_Sub$A_T<1 , PBSC_Sub$Base, NA )
PBSC_Subset=data.frame(PBSC_Subset)

PBSC_DIV=Splitting_data(PBSC_Subset,Name_of_files,50)
Good_data=Finding_Continuous_Line(PBSC_Subset,Name_of_files)
Good_data_DIV=Splitting_data(Good_data,Name_of_files,50)
for (i in 1:(length(PBSC_DIV))){
  Graph=PBSC_DIV[[i]]
  PBSC_Subset_Plot<-ggplot(Graph,aes(x=Base,y=Name)) +
    geom_point(aes(colour = Legend),alpha = 0.5,size = 1.5,position = position_jitter(width = 0.25, height = 0))+
    ggtitle("Per Base Sequence Content Complete")
  g_plot[[i]]=PBSC_Subset_Plot
}

for (i in 1:(length(g_plot))){
  print  (g_plot[[i]])
}
g_plot=NULL
for (i in 1:(length(Good_data_DIV))){
  Graph=Good_data_DIV[[i]]
  Good_data_Plot=ggplot(Graph,aes(x=Base,y=Name)) +
    geom_point(aes(colour = Legend),alpha = 0.5,size = 1.5,position = position_jitter(width = 0.25, height = 0))+
    ggtitle("Per Base Sequence Content Longest Lines")
  g_plot[[i]]=Good_data_Plot
}
for (i in 1:(length(g_plot))){
  print  (g_plot[[i]])
}

####################################################################################################################################################################
# Basic Statistics
Basic_Stat_Merged=NULL
BasicS=NULL
list_of_Ids=NULL
for (Name in Name_of_files){
  Basic_Stat=data [[paste('./',Name,sep='')]][['Basic_Statistics.txt']]
  BasicS$Total_Sequences=Basic_Stat[[2]][4]
  BasicS$Name <- Name
  Basic_Stat_Merged <- rbind(Basic_Stat_Merged,BasicS)
}
len=nrow(Basic_Stat_Merged)
row.names(Basic_Stat_Merged)<-c(1:len)
Basic_Stat_Merged=data.frame(Basic_Stat_Merged)

for ( i in (1:len)){
  name=Basic_Stat_Merged$Name[[i]]
  info=getReadPairIDFromFq(name)
  list_of_Ids[i]=info$readPairID
  Basic_Stat_Merged$ID[i]=info$readPairID
}
list_of_Ids=unique(list_of_Ids)
Data=NULL
for ( i in list_of_Ids){
  Basic_Stat_Subset=subset(Basic_Stat_Merged,Basic_Stat_Merged$ID==i)
  Test=length(unique(Basic_Stat_Subset$Total_Sequences))
  if (Test == 1){
    Basic_Stat_Subset$Verification=1
  }
  else{
    Basic_Stat_Subset$Verification=0
  }
  Basic_Stat_Subset$Total_Sequences <- as.numeric(Basic_Stat_Subset$Total_Sequences)
  Data=rbind(Data,Basic_Stat_Subset)
}

PBSC_Subset_Plot<-ggplot(Data,aes(x=ID,y=Total_Sequences)) 
PBSC_Subset_Plot + geom_point(aes(color = Verification),alpha = 0.5,size = 1.5,position = position_jitter(width = 0.25, height = 0))+
  ggtitle("Per Base Sequence Content")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
####################################################################################################################################################################
# Basic Statistics
Basic_Stat_Merged=NULL
BasicS=NULL
list_of_Ids=NULL
i=1
for (Name in Name_of_files){
  Basic_Stat=data [[paste('./',Name,sep='')]][['Basic_Statistics.txt']]
  BasicS$Total_Sequences=Basic_Stat[[2]][4]
  BasicS$Name <- Name
  Basic_Stat_Merged <- rbind(Basic_Stat_Merged,BasicS)
}
len=nrow(Basic_Stat_Merged)
row.names(Basic_Stat_Merged)<-c(1:len)
Basic_Stat_Merged=data.frame(Basic_Stat_Merged)

for ( i in (1:len)){
  name=Basic_Stat_Merged$Name[[i]]
  info=getReadPairIDFromFq(name)
  list_of_Ids[i]=info$readPairID
  Basic_Stat_Merged$ID[i]=info$readPairID
}
list_of_Ids=unique(list_of_Ids)
Data=NULL
for ( i in list_of_Ids){
  Basic_Stat_Subset=subset(Basic_Stat_Merged,Basic_Stat_Merged$ID==i)
  Test=length(unique(Basic_Stat_Subset$Total_Sequences))
  if (Test == 1){
    Basic_Stat_Subset$Verification=1
  }
  else{
    Basic_Stat_Subset$Verification=0
  }
  Basic_Stat_Subset$Total_Sequences <- as.numeric(Basic_Stat_Subset$Total_Sequences)
  Data=rbind(Data,Basic_Stat_Subset)
}

PBSC_Subset_Plot<-ggplot(Data,aes(x=ID,y=Total_Sequences)) 
PBSC_Subset_Plot + geom_point(aes(color = Verification),alpha = 0.5,size = 1.5,position = position_jitter(width = 0.25, height = 0))+
  ggtitle("Per Base Sequence Content")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
####################################################################################################################################################################
#Per Sequence GC Content
PSGCC_Merged=NULL
for (Name in Name_of_files){
  PSGCC=data [[paste('./',Name,sep='')]][['Per_sequence_GC_content.txt']]
  PSGCC$Name <- Name
  PSGCC_Merged <- rbind(PSGCC_Merged,PSGCC)
}
PSGCC_Merged=rename(PSGCC_Merged, c("#GC Content"="GC_Content"))
Bad_Data=Normality_Test(PSGCC_Merged,Name_of_files,npar=TRUE)

#dev.off()
####################################################################################################################################################################