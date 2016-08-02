#libraries
library(plyr)
library(ggplot2)
library(xtable)
library(tcltk)
library(asbio)
library(nortest)
library(data.table)
#Calculates the Confidence Intervals
Calculate_Confidence_Intervals_of_Table <- function(Merge_File,x_value,npar=TRUE){
  colnames(Merge_File)<-c("point","value","Name")
  Bad_data=NULL
  BG=NULL
  Median=NULL
  Lower=NULL
  Upper=NULL
  Data=NULL
  Median_List=NULL
  Upper_List=NULL
  Lower_List=NULL
  Con_Interval_Table=NULL
  # For loop that loops through all the quality scores x being th quality score.
  # The for loop Subsets all the quality scores, takes the confidence interval using ci.median and then stores it.
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
    Median_List=rbind(Median_List,Median)
    Lower_List=rbind(Lower_List,Lower)
    Upper_List=rbind(Upper_List,Upper)
    for (i in nrow(Quality_Subset)){
      if (Quality_Subset$value[i] <= Lower$value[1] || Quality_Subset$value[i] >= Upper$value[1]){
        Subset_by_name=subset(Merge_File,Merge_File$Name==Quality_Subset$Name[i])
        Merge_File=Merge_File[Merge_File$Name != Quality_Subset$Name[i],]
        Subset_by_name=data.frame(Subset_by_name)
        Bad_data=rbind(Bad_data,Subset_by_name)
      }
    }
  }
  Data=list(Median_List,Lower_List,Upper_List,Bad_data)
  #Bad_data$value <- as.numeric(Bad_data$value)
  return(Data)
}

# Checking if its higher that 35
Verifying_Intervals <- function(Bad_samples,npar=TRUE){
  colnames(Bad_samples)<-c("Quality","Count","Name","Lower","Upper")
  Verifying_Bad_Sample=NULL
  Subset_of_35=subset(Bad_samples,Bad_samples$Quality>=35)
  Backup=Bad_samples
  for (i in 1:nrow(Subset_of_35)){
    if (Subset_of_35$Count[i] >= Subset_of_35$Upper[i]){
      Subset_Questionable_sample=subset(Bad_samples,Bad_samples$Name==Subset_of_35$Name[i])
      Bad_samples=Bad_samples[Bad_samples$Name != Subset_of_35$Name[i],]
      Subset_Questionable_sample=data.frame(Subset_Questionable_sample)
      Verifying_Bad_Sample=rbind(Verifying_Bad_Sample,Subset_Questionable_sample)
    }
  }
  Subset_of_35_Under=subset(Verifying_Bad_Sample,Verifying_Bad_Sample$Quality<35)
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
  return(Bad_samples)
}

#Finding continuous lines and verifiying they are straigh in y axis
Finding_Continuous_Line<- function(PBSC_Subset,Name_of_files){
  colnames(PBSC_Subset)<-c("Base","G","A","T","C","Name","G_C","A_T","GC_bot")
  Longest_Lines=NULL
  for (x in Name_of_files){
    Straightlines=NULL
    Verify_Mean=NULL
    Subset_Name=subset(PBSC_Subset,PBSC_Subset$Name==x)
    TF_List=is.na(Subset_Name$Base)
    TF_List=c(TF_List,TRUE)
    Index_Trues=which(TF_List %in% TRUE)
    for ( i in 1:(length(Index_Trues)-1)){
      if (Index_Trues[i+1]-Index_Trues[i]-1>1){
        Straight=matrix(c(Index_Trues[i]+1,Index_Trues[i+1]-1,Index_Trues[i+1]-Index_Trues[i]-1),byrow = TRUE,ncol = 3)
        Straight=data.frame(Straight)
        Straightlines=rbind(Straightlines,Straight)
      }
    }
    if (!is.null(Straightlines)){
      Straightlines=data.frame(Straightlines)
      bool=TRUE
      while(bool){
        Longest_line=Straightlines[(which.max( Straightlines[,3])),]
        colnames(Longest_line)<-c("L","H","D")
        Verify_Mean=Subset_Name[Subset_Name$Base>=Longest_line$L[1] & Subset_Name$Base<=Longest_line$H[1],]
        Verify_Mean=na.omit(Verify_Mean)
        rownames(Verify_Mean)<-c(1:nrow(Verify_Mean))
        Verify_Mean=data.frame(Verify_Mean)
        MeanAT=mean(Verify_Mean$A)
        MeanGC=mean(Verify_Mean$G)
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
    else{
      NA_Data=c(NA,NA,NA,NA,NA,Subset_Name$Name[i],NA,NA,NA)
      Longest_Lines=rbind(Longest_Lines,NA_Data)
    }
  }
  return(Longest_Lines)
}
# Spliting into 50 per set
Splitting_data <- function (Good_data,Name_of_files,npar=TRUE){
  len = length(Name_of_files)
  rest=len %% 10
  value = len-rest
  Datalist=NULL
  value_50=value/50
  Data_sep_50=NULL
  t=1
  for (i in (1:(t*50))){
    name=Name_of_files[i]
    Data_subset=subset(Good_data,Good_data$Name==name)
    Data_sep_50=rbind(Data_sep_50,Data_subset)
  }
  Datalist[[t]]=Data_sep_50
  Data_sep_50=NULL
  t=t+1
  while(t!=value_50){
    for (i in ((t*50+1):(t*50+50))){
      name=Name_of_files[i]
      Data_subset=subset(Good_data,Good_data$Name==name)
      Data_sep_50=rbind(Data_sep_50,Data_subset)
    }
    Datalist[[t]]=Data_sep_50
    Data_sep_50=NULL
    t=t+1
  }
  for (i in ((t*50+1):(t*50+rest))){
    name=Name_of_files[i]
    Data_subset=subset(Good_data,Good_data$Name==name)
    Data_sep_50=rbind(Data_sep_50,Data_subset)
  }
  Datalist[[t]]=Data_sep_50
  return(Datalist)
}
#Normality test for Per Sequence GC Content
Normality_Test <- function(Merge_File,Name_of_files,npar=TRUE){
  colnames(Merge_File)<-c("point","value","Name")
  Bad_data=NULL
  for (x in Name_of_files){
    Data_Subset=subset(Merge_File,Merge_File$Name==x)
    Data=Data_Subset$value
    print(Data)
    Shapiro=ad.test(Data)
    Shapiro=Shapiro$p.value
    print(Shapiro)
    if (Shapiro<.05){
      Bad_data=rbind(Bad_data,Data_Subset)
    }
  }
  #Bad_data$value <- as.numeric(Bad_data$value)
  return(Bad_data)
}



# function to extract read pair id from fastq file name or tar or fastqc result

miserableTestData=c("140618_I1069_FCC4MWGACXX_L4_LSCRHmM121ACAAAPEI-3_1.clean_fastqc", "AH07020812_ACACGA_L002_R1_002.fastq", "E021_0001_20140916_tumor_RNASeq_R1.clean.tar","C021_0001_20140916_tumor_RNASeq.tar", "LS3_GGCTAC_L002_R2_002.fq.gz", "AM163062014_2.fastq.gz", "476_R1.fastq.gz", "486_R1.trimmed.fastq.gz", "home/hbeale/BS35112812_CACACA_L003_001.tar", "C021_0003_001409_BR_Whole_T3_TSMRU_A07217_R2_fastqc")
readEndExtractionExpected=c(T, T, T, F, T, T, T, T, F, T)

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


options(stringsAsFactors=FALSE)
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
        data[[folder]][['Overrepresented_sequences.txt']]=data[[folder]][[Info]]=read.table(paste(File,"/",folder,"/",Name_of_folder,"/",Info,sep = ''),header=TRUE,comment.char = "",sep="\t",check.names = FALSE)
    }
    else
      data[[folder]][[Info]]=read.table(paste(File,"/",folder,"/",Name_of_folder,"/",Info,sep = ''),header=TRUE,comment.char = "",sep="\t",check.names = FALSE)
  }
}
name_of_txt=append(name_of_txt, "Total Duplicate Percentage")
#Name of of pdf file where all graphs will be saved
#pdf(file = "/Users/Ilian/Desktop/fastQC_treehouse_results/Merged.pdf", title="Graphs of Merged Data")
#ggplot Per_sequence_quality_scores
PSQS_Merged=NULL
i=1
for (Name in Name_of_files){
  PSQS=data [[paste('./',Name,sep='')]][['Per_sequence_quality_scores.txt']]
  PSQS$Name <- Name
  PSQS_Merged <- rbind(PSQS_Merged,PSQS)
}
PSQS_Merged=rename(PSQS_Merged, c("#Quality"="Quality"))
Quality_values_PSQS=c(2:40)
Conf_Int_Table=Calculate_Confidence_Intervals_of_Table(PSQS_Merged,Quality_values_PSQS)
Lower=Conf_Int_Table[[2]]
Upper=Conf_Int_Table[[3]]
Bad_samples=Conf_Int_Table[[4]]
Bad_samples$Lower <- Lower$value
Bad_samples$Upper <- Upper$value
Bad_samples=data.frame(Bad_samples)
colnames(Bad_samples)<-c("Quality","Count","Name","Lower","Upper")

Bad_Data=Verifying_Intervals(Bad_samples)


PSQS_PLOT_Conf_Int <- ggplot(Bad_Data, aes(x = Quality, y = Count,group=Name))
PSQS_PLOT_Conf_Int + geom_line(aes(color =Name))+ggtitle("Per Sequence Quality Scores 95% Confidence Intervals Excluding Upper Outliers Above Quality 35")+
  geom_ribbon(aes(ymin=Lower,ymax=Upper), alpha=0.01, linetype=1,colour="grey60", size=.05,fill="grey40")

PSQS_PLOT_Conf_Int_IND <- ggplot(Bad_Data, aes(x = Quality, y = Count)) 
PSQS_PLOT_Conf_Int_IND + geom_line(aes(color =Name)) +facet_wrap(~Name) + geom_ribbon(aes(ymin=Lower,ymax=Upper), alpha=0.1, linetype=1,colour="grey20", size=.05,fill="grey20")

#ggplot per base sequence content plot
PBSC_Merged=NULL
PBSC_Sub=NULL
PBSC_Subset=NULL
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
PBSC_Sub$GC_bot = ifelse(PBSC_Merged$C<PBSC_Merged$T | PBSC_Merged$C<PBSC_Merged$A
                         | PBSC_Merged$G<PBSC_Merged$T | PBSC_Merged$G<PBSC_Merged$A,1,0)
PBSC_Sub=data.frame(PBSC_Sub)

PBSC_Subset=PBSC_Sub
PBSC_Subset$Base=ifelse(PBSC_Sub$G_C<1 & PBSC_Sub$A_T<1 , PBSC_Sub$Base, NA )
PBSC_Subset=data.frame(PBSC_Subset)

PBSC_DIV=Splitting_data(PBSC_Subset,Name_of_files)
Good_data=Finding_Continuous_Line(PBSC_Subset,Name_of_files)
Good_data_DIV=Splitting_data(Good_data,Name_of_files)

for (i in length(PBSC_DIV)){
  PBSC_Subset_Plot<-ggplot(PBSC_DIV[[i]],aes(x=Base,y=Name)) 
  PBSC_Subset_Plot + geom_point(aes(color = GC_bot),alpha = 0.5,size = 1.5,position = position_jitter(width = 0.25, height = 0))+
    ggtitle("Per Base Sequence Content Complete")
}
for (i in length(Good_data_DIV)){
  Good_data_Plot<-ggplot(Good_data_DIV[[i]],aes(x=Base,y=Name)) 
  Good_data_Plot + geom_point(aes(color = GC_bot),alpha = 0.5,size = 1.5,position = position_jitter(width = 0.25, height = 0))+
    ggtitle("Per Base Sequence Content Longest Lines")
}

PBSC_Subset_Plot<-ggplot(PBSC_Subset,aes(x=Base,y=Name)) 
PBSC_Subset_Plot + geom_point(aes(color = GC_bot),alpha = 0.5,size = 1.5,position = position_jitter(width = 0.25, height = 0))+
  ggtitle("Per Base Sequence Content")

# Basic Statistics
Basic_Stat_Merged=NULL
BasicS=NULL
i=1
for (Name in Name_of_files){
  Basic_Stat=data [[paste('./',Name,sep='')]][['Basic_Statistics.txt']]
  BasicS$Total_Sequences=Basic_Stat[[2]][4]
  BasicS$Name <- Name
  Basic_Stat_Merged <- rbind(Basic_Stat_Merged,BasicS)
}
len=nrow(Basic_Stat_Merged)
row.names(Basic_Stat_Merged)<-c(1:len)


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
