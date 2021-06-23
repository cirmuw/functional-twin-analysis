library(readxl)
library(xlsx)


readRenviron("../settings.Renviron")
sexes<-Sys.getenv("HCP_information_sheet_unrestricted")
subjects<-Sys.getenv("HCP_information_sheet")

sexes=read_excel(sexes,1)
sexes=sexes[-113,] # in our case row 113 had do be removed to match the sheet containing zygosity information
sexes=sexes[,3,drop=FALSE]

subjects=read_excel(subjects,1)
subjects=subjects[,c(1,4,5),drop=FALSE]

joined=cbind(subjects,sexes)
joined=joined[joined$Zygosity!='NotTwin',, drop=FALSE]
joined=joined[order(joined$Zygosity,joined$Mother_ID), ,drop=FALSE]

#remove subjects for which fMRI data was not available
gzf='../../Deliveries/Subjects.csv'
data=read.csv(gzf,header=T)

joined=joined[joined$Subject %in% data$Subject,, drop=FALSE]
joined=joined[,4,drop=FALSE]
write.xlsx(joined,"../../Deliveries/sexes.xlsx")
