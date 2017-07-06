## File containing the sets of the optimization problem

rm(list=ls())

gamspath <<- '/opt/gams/24.8'

options("width"=200)
options(digits=4)

bilevel <- function(n,s,sc,case,met){
setwd('data')
load(paste("data_","n",n,"_","m",n,"_","p",n/2,"_","q",n/2,"_","r",n/2,"_","s",s,"_","sc",sc,"_","case",case,".Rda",sep=""))
input <<- data
setwd('..')
if (met==1) {source('BB.r')}
if (met!=1) {
source("interfaceR_GAMS.R")
#SETS
input_set('n',input$n)
input_set('m',input$m)
input_set('p',input$p)
input_set('q',input$q)
input_set('r',input$r)
#SCALARS
input_sca('met',met)
#PARAMETERS
input_par('c(n)',input$c)
input_par('d(m)',input$d)
input_par('A(p,n)',input$A)
input_par('b(p)',input$b)
input_par('e(m)',input$e)
input_par('F(q,n)',input$F)
input_par('G(q,m)',input$G)
input_par('h(q)',input$h)
input_par('I(r,m)',input$I)
input_par('j(r)',input$j)
output('OF_UL')
output('u1(m)')
output('u2(q)')
output('u3(r)')
solveGAMS("Bilevel")
file.remove(list.files(pattern = 'Rda'))
setwd('..')
unlink(paste(files_id), recursive = TRUE)
}
print(paste('Method:',met))
print(paste('Objective function:',OF_UL))
print(paste('Time:',time,'seconds'))
if ((modelstat==1)|(modelstat==2)|(modelstat==8)) {print('Status: 1')} else {print('Status: 0')}
}
