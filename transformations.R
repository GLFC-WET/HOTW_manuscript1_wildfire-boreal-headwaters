############### Bio Transformations #############################
# needs:
# ftms3
# data/transformations.csv
## should be able to run in totality once set up

library(readr)

ftms3 <- read_csv("data/ftms3.csv")

#A function to indicate if R returns "integer(0)" as output - Erin, you don't have to touch that. 
is.integer0 <- function(x)
{is.integer(x) && length(x) == 0L}

#The function that will extract the transformation in the dataset - Erin, you don't need to change that, it is just a function
Putative_Formulas <- function(x) {
  Formulas_Present <- which(!is.na(x))#Get the formulas present in the samples.
  Production <- lapply(paf$mz[Formulas_Present],FUN=function(x)(x+Transformation$Mass))#Get all possible Production
  Degradation <- lapply(paf$mz[Formulas_Present],FUN=function(x)(x-Transformation$Mass))#Get all possible degradation
  Putative_Production <- lapply(Production,FUN = function(x)(which(paf$mz[Formulas_Present]%in%x)))#Check which Production occurs
  Putative_Degradation <- lapply(Degradation,FUN = function(x)(which(paf$mz[Formulas_Present]%in%x)))#Check which degradation occurs
  
  Number_Production <- length(unlist(Putative_Production))#Number of Production 
  Number_Degradation <- length(unlist(Putative_Degradation))#Number of Degradation 
  Formulas_with_Transformation <- which(!sapply(Putative_Production, function(x) is.integer0(x))) #Which formulas were transformed
  
  Transformation_Type <- lapply(Degradation,FUN = function(x)(which(!is.na(match(x,paf$mz[Formulas_Present])))))#What type of transformation per formulas
  DIM_Lists <- unlist(sapply(Transformation_Type,function(x) length(x)))#The number of transformation per formulas
  DIM_Lists <- DIM_Lists[Formulas_with_Transformation]#The number of transformation for each formulas (without formulas with 0 transformation)
  
  # Number_Transformation_per_Mass <- (cbind(paf$mz[Formulas_Present][Formulas_with_Transformation],
  #                                             DIM_Lists[Formulas_with_Transformation]))#Number of transformation for each mass
  
  Transformation_Type <- unlist(Transformation_Type)
  Matrix_Transformation_Names <- matrix(nrow=length(Formulas_with_Transformation),ncol=nrow(Transformation))
  colnames(Matrix_Transformation_Names)<-Transformation$Name
  rownames(Matrix_Transformation_Names)<-paf$mz[Formulas_with_Transformation]
  if(length(Formulas_with_Transformation)>1){
    for (i in 1:length(Formulas_with_Transformation)){
      if(i==1)(Transformation_names <- Transformation_Type[1:DIM_Lists[i]])
      if(i==1)(j <-  DIM_Lists[i])
      if(i!=1)(j <- j+DIM_Lists[i])
      if(i!=1)(Transformation_names <- Transformation_Type[j:(j+DIM_Lists[i])])
      Transformation_names <- unique(Transformation$Name[Transformation_names])
      Matrix_Transformation_Names[i,colnames(Matrix_Transformation_Names)%in%Transformation_names]<-1
    }
  }else{Matrix_Transformation_Names<-0}
  
  return(Matrix_Transformation_Names)
  on.exit(stopCluster(cl))
}

library(parallel)

paf <- ftms3[,1:83] 
Crosstab.2 <- ftms3[,84:169] # my samples
dim(Crosstab.2)
# Crosstab.2 <- t(Crosstab.2) # transposed the table
Transformation <-  read.csv("data/transformations.csv")

# detectCores() to figure out how many cores you have
cl <- makeCluster(getOption("cl.cores", 7)) # I have 16 cores
clusterExport(cl=cl, varlist=c("paf","Transformation","is.integer0"))
Crosstab_List <- as.list(data.frame(Crosstab.2)) # Convert the matrix into a list - mandatory for parallel computing 
Putative_Degradation_Matrix_List <- parLapply(cl=cl,Crosstab_List,
                                              fun=Putative_Formulas)
stopCluster(cl)

saveRDS(Putative_Degradation_Matrix_List,"C:/Users/erinmatula/Downloads/Putative_Trannsform_Matrix_List.rds")
Putative_list <- readRDS("C:/Users/erinmatula/Downloads/Putative_Trannsform_Matrix_List.rds")


# Putative_list2 <- Putative_list[-c(5,122,123)] #for subsetting, took out a few columns

for(i in 1:length(Putative_list)){
  x <- unlist(Putative_list[[i]])
  dim(x)
  if(is.null(dim(x)))(print(i))
}

# take out any blanks or samples you dont want in
Putative_list <- Putative_list[-c(72,78)]

length(Putative_list)
Masses <- unlist(lapply(Putative_list,FUN = function(x) (rownames(x))))
Process <- (lapply(Putative_list,FUN = function(x) (rowSums(x,na.rm=T))))
Masses <- unique(Masses)
Process2 <- unname(Process)

Table_Process <- matrix(nrow=length(Masses),ncol=length(names(Putative_list)))
rownames(Table_Process)<- Masses
colnames(Table_Process)<- names(Putative_list)



for (i in 1:length(names(Putative_list))){
  Names <- names(unlist(Process2[i]))
  X <- which(Names%in%Masses)
  Y <- which(Masses%in%Names)  
  Table_Process[Y,i]<- unlist(Process[i])[X]
}


colSums(Table_Process,na.rm = TRUE) # transformation per samples
# can run row names to find formula

Process[[1]] # one sample, which formula has been transformed in this one sample
