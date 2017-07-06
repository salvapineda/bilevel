# Solving branch and bound in R

library('lpSolve')

tol <<- 10^-10
maxtime <- 6*3600

n <- input$n
m <- input$m
p <- input$p
q <- input$q
r <- input$r
l <- q+r+m

A1 <- rbind(input$A,-diag(n))
B1 <- matrix(0,p+n,n)
b1 <- c(input$b,rep(0,n))
A2 <- rbind(input$F,matrix(0,r,n),matrix(0,m,n))
B2 <- rbind(input$G,input$I,-diag(m))
b2 <- c(input$h,input$j,rep(0,m))

#function to solve lp
solve_lp <- function(A,b,c){   
sol <- lp ("min", c, A, rep("<=",nrow(A)), b)
return(list(ms=sol$status,of=sol$objval,opt=sol$solution))
}

#solve lp with w and z
solve_lp_wz <- function(w,z){
A <- rbind(cbind(A1,B1,matrix(0,p+n,l)),cbind(matrix(0,m,n),matrix(0,m,m),t(B2)),-cbind(matrix(0,m,n),matrix(0,m,m),t(B2)),cbind(A2,B2,matrix(0,l,l)),-cbind(matrix(0,l,n),matrix(0,l,m),diag(l)),cbind(matrix(0,l,n),matrix(0,l,m),diag(w)),-cbind(diag(z)%*%A2,diag(z)%*%B2,matrix(0,l,l)))
b <- c(b1,-input$e,input$e,b2,-rep(0,l),rep(0,l),-c(diag(z)%*%b2))
c <- c(input$c,input$d,rep(0,l))
sol <- solve_lp(A,b,c)
com <- diag(sol$opt[(n+m+1):(n+m+l)])%*%(b2-A2%*%sol$opt[1:n]-B2%*%sol$opt[(n+1):(n+m)])
index <- which.max(com)
return(list(ms=sol$ms,of=sol$of,opt=sol$opt,com=com,index=index))
}

#branching
branch_bound <- function(matrix,l){
if (length(which(matrix[,2*l+2] < tol))==0) {best_int <- Inf} else {
int_sol <- which(matrix[,2*l+2] < tol)
best_int <- min(matrix[int_sol,2*l+1])
}
max_com <- max(matrix[,2*l+2])
if (max_com > tol){
index <- which.max(matrix[,2*l+2])
branch <- matrix[index,2*l+3]
w1 <- matrix[index,1:l]
w1[branch] <- 1
z1 <- matrix[index,(l+1):(2*l)]
w2 <- matrix[index,1:l]
z2 <- matrix[index,(l+1):(2*l)]
z2[branch] <- 1
matrix <- matrix[-index,]
output <- solve_lp_wz(w1,z1)
if ((output$ms==0) & (output$of < best_int)){matrix <- rbind(matrix,matrix(c(w1,z1,output$of,max(output$com),output$index),nrow=1))}
output <- solve_lp_wz(w2,z2)
if ((output$ms==0) & (output$of < best_int)){matrix <- rbind(matrix,matrix(c(w2,z2,output$of,max(output$com),output$index),nrow=1))}
}
return(matrix)
}

#time start
start.time <- Sys.time()
current.time <- Sys.time()

#add infeasible solution
nodes <- matrix(c(rep(0,l),rep(0,l),Inf,0,0),nrow=1)

#relaxed problem
w <- rep(0,l)
z <- rep(0,l)
output <- solve_lp_wz(w,z)
if (output$ms==0) {nodes <- rbind(nodes,matrix(c(w,z,output$of,max(output$com),output$index),nrow=1))}
#feasible solution
sol <- solve_lp(B2,b2-A2%*%output$opt[1:n],c(input$e))
z <- as.numeric(abs(B2%*%sol$opt-b2+A2%*%output$opt[1:n]) < tol)
w <- as.numeric(abs(B2%*%sol$opt-b2+A2%*%output$opt[1:n]) > tol)
output <- solve_lp_wz(w,z)
if (output$ms==0) {nodes <- rbind(nodes,matrix(c(w,z,output$of,max(output$com),output$index),nrow=1))}

#iteration
while ((max(nodes[,2*l+2]) > tol) & (as.numeric(Sys.time() - start.time, units="secs") < maxtime)) {nodes <- branch_bound(nodes,l) }

if (nrow(nodes) > 1) {modelstat <- 1;solvestat <- 1} else {modelstat <- 4;solvestat <- 1}

time <- as.numeric(Sys.time() - start.time, units="secs")

OF_UL <- min(nodes[which(nodes[,2*l+2]<tol),2*l+1])


















