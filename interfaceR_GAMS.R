#preamble
library(gdxrrw)
library(xtable)
library(parallel)
igdx(gamspath,silent=TRUE)

#create subfolder
files_id <<- sample(1:10^7,1)
dir.create(paste(files_id))
files2move <- list.files(pattern = 'gms')
file.copy(from = files2move, to = paste(files_id))
files2move <- list.files(pattern = 'Rda')
file.copy(from = files2move, to = paste(files_id))
files2move <- list.files(pattern = '.r')
file.copy(from = files2move, to = paste(files_id))
setwd(paste(files_id))
alldata <<- array(list(NULL),c(0)) 
result <<- array(list(NULL),c(0)) 
fdata1 <- file(".filedata1.txt", "w")
fdata2 <- file(".filedata2.txt", "w")
cat(file = fdata2, "$gdxin %infile% \n")
cat(file = fdata2, "$load ")
fsol1 <- file(".filesol1.r", "w")
fsol2 <- file(".filesol2.txt", "w")
cat(file = fsol2, "execute_unload '%outfile%'")

# input set
input_set <- function(name_set,value_set){
alldata[[length(alldata)+1]] <<- list(name=name_set, type='set', uels=list(paste(name_set,1:value_set, sep=""))) 
cat(file = fdata1, "SET ",name_set, ";\n", sep = "")
cat(file = fdata2, name_set," ", sep = "")
}

#input scalar
input_sca <- function(name_sca,value_sca){
alldata[[length(alldata)+1]] <<- list(name=name_sca, type='parameter', val=value_sca, dim=0, form='full') 
cat(file = fdata1,"PARAMETER ",name_sca,";\n",sep = "")
cat(file = fdata2, name_sca," ", sep = "")
}

#input parameter
input_par <- function(name_par,value_par){
labels<-c()
par_left <- regexpr(pattern ='\\(',name_par)
par_right <- regexpr(pattern ='\\)',name_par)
name_par2 <- substr(name_par,1,par_left-1)
num.indexes <- (par_right-par_left)/2
for (k in 1:num.indexes){
set_index <- which(sapply(alldata, function(X) any(match(X$name,substr(name_par, par_left - 1 + 2*k , par_left - 1 + 2*k)))))
labels<-c(labels,alldata[[set_index]]$uels)
}
alldata[[length(alldata)+1]] <<- list(name=name_par2, type='parameter', val=value_par, dim=num.indexes, form='full', uels=labels) 
cat(file = fdata1,"PARAMETER ",name_par,";\n",sep = "")
cat(file = fdata2, name_par2," ", sep = "")
}

# solveGAMS                        
solveGAMS <- function(namefile){
cat(file = fdata2, "\n$gdxin \n")
close(fdata2)
close(fdata1)
cat(file = fsol2, ";")
close(fsol2)
close(fsol1)
wgdx (".inputdata.gdx",alldata)
rc <- system (paste(gamspath,"/gams ",namefile," --infile= .inputdata.gdx --outfile= .outputsol.gdx lo=2",sep=""))
source(".filesol1.r")
if (0 != rc) { stop(paste("Bad return from gams: wanted 0, got",rc))
} else { print (paste("GAMS call succeeded (modelstat =",modelstat, ",solvestat=",
solvestat,")")) }
result <<- array(list(NULL),c(0))
}

# output
output <- function(name_out){
if (!grepl('\\(',name_out)){
cat(file = fsol2, paste(", ",name_out,sep=""))
}
if (grepl('\\(',name_out)){
cat(file = fsol2, paste(", ",substr(name_out,1,regexpr(pattern ='\\(',name_out)-1),sep=""))
}
cat(file = fsol1, paste("read_value('",name_out,"')\n",sep=""))
}

# read value
read_value <- function(name_out){
if (!grepl('\\(',name_out)){
lst <- list(name=name_out,form='full',compress=FALSE)
result[[length(result)+1]] <<- rgdx(".outputsol.gdx", lst)
assign(name_out,result[[length(result)]]$val,envir = .GlobalEnv)
}
if (grepl('\\(',name_out)){
labels<-c() 
par_left <- regexpr(pattern ='\\(',name_out)
par_right <- regexpr(pattern ='\\)',name_out)
name_out2 <- substr(name_out,1,par_left-1)
num.indexes <- (par_right-par_left)/2
for (k in 1:num.indexes){
set_index <- which(sapply(alldata, function(X) any(match(X$name,substr(name_out, par_left - 1 + 2*k , par_left - 1 + 2*k)))))
labels[k] <- alldata[[set_index]]$uels
}
lst <- list(name=name_out2, form='full',uels=labels)
result[[length(result)+1]] <<- rgdx(".outputsol.gdx", lst)
assign(name_out2,as.data.frame(result[[length(result)]]$val),envir = .GlobalEnv)
}
}

output('modelstat')
output('solvestat')
output('time')
output('optca')
output('optcr')




