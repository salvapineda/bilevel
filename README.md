# bilevel

The purpose of this code is to solve linear bilevel programming problems using different methodologies. All the details about the solution methods is explained in [1] and the references therein:

Instructions:

1. Download and install the programming language R version 3 (https://cran.r-project.org/)

2. Download and install GAMS (https://www.gams.com/)

3. Install the R-Gams interface (https://www.gams.com/latest/docs/tools/gdxrrw/index.html). Make sure that the variable gamspath is set to the directory where the executable GAMS file is located. You can modify that in the file Bilevel.r

4. Download the zip of this repository and extract in a directory of your choice.

5. Start R, move to the directory where you extracted the code, and run source('Bilevel.r')

6. Run bilevel(n,s,sc,case,met) where
   - n is the size of the problem. You can choose 50, 100 or 200.
   - s is the sparcity of the matrices. You can choose 100 or 50.
   - sc is the scalability of the problem. You can choose 0 or 1.
   - case refers to each randomly generated number. You can choose 1,2,3,...,99,100.
   - met refers to the solution approach. You can choose one of the following options:
          1: Branch and bound algorithm
          2: Special order set algorithm
          3: Fortuny-Amat approach with big-M = 5 
          4: Fortuny-Amat approach with big-M = 10 
          5: Fortuny-Amat approach with big-M = 20 
          6: Fortuny-Amat approach with big-M = 50 
          7: Fortuny-Amat approach with big-M = 100 
          8: Fortuny-Amat approach with big-M = 200 
          9: Fortuny-Amat approach with big-M = 500 
          10: Fortuny-Amat approach with big-M = 1000 
          11: Fortuny-Amat approach with big-M = 5000 
          12: Fortuny-Amat approach with big-M = 10000 
          13: Fortuny-Amat approach with big-M = 100000 
          14: Regularization method
          15: Penalty method
          16: Solution approach proposed in [1] with M = 2
          17: Solution approach proposed in [1] with M = 5
          18: Solution approach proposed in [1] with M = 10

Remarks:

  - This repository only includes the data corresponding to the 900 random problems used in [1]
  - To use GAMS you need a license. The results in [1] have been calculated using CPLEX as the linear and mixed-integer solver, and CONOPT as the non-linear solver. However, other solvers may be used as well.  
  
References:
  
[1] S. Pineda, H. Bylling, J.M. Morales, "Efficiently solving linear bilevel programming problems using off-the-shelf optimization software," submitted to optimization and engineering.
