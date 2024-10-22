# ERGO CODE: IP_PMM
This is an Interior Point-Proximal Method of Multipliers, suitable for solving linear and convex quadratic
programming problems. The method takes as input a problem of the following form:


![equation](https://latex.codecogs.com/gif.latex?%24%5C%20%24%5C%5C%20%5Cmin%5C%20c%5ETx%20&plus;%20%5Cfrac%7B1%7D%7B2%7Dx%5ET%20Q%20x%2C%20%5C%5C%20%5Ctext%7Bs.t.%7D%20Ax%20%3D%20b%2C%5C%5C%20x_C%20%5Cgeq%200%2C%20%5Cforall%5C%20i%20%5Cin%20C%20%5Csubset%20%5C%7B1%2C%5Ccdots%2Cn%5C%7D%2C%5C%5C%20x_I%5C%20%5Ctext%7Bfree%7D%2C%20%5Cforall%5C%20i%20%5Cin%20F%20%3D%20%5C%7B1%2C%5Ccdots%2Cn%5C%7D%5Csetminus%20C)  
 
                     
and solves it to optimality, returning the primal and dual optimal solutions (a message indicating that the
optimal solution was not found or an infeasibility indicator).

INPUT PARAMETERS:

IP_PMM(c, A, Q, b): 

                     find the optimal solution of the problem, with an error tolerance of 10^(-6)
                     Upon success, the method returns x (primal solution), y (Lagrange multipliers) and
                     z  (dual optimal slack variables). If the run was unsuccessful, the method  either returns
                     a certificate of infeasibility, or terminates after 100 iterations. 
                     
IP_PMM(c, A, Q, b, free_variables): 

                      The last parameter is a matrix of indices, pointing to the free variables of the
                      problem. If not given, it is assumed that there are no free variables.
                                     
IP_PMM(c, A, Q, b, free_variables, tol): 

                      This way, the user can specify the tolerance to which the problem is solved.

IP_PMM(c, A, Q, b, free_variables, tol, max_it):

                      This way, the user can also specify the maximum number of iterations.

IP_PMM(c, A, Q, b, free_variables, tol, maxit, pc):


                      predictor-corrector option.

                      false: no predictor-corrector.
                                                     
                      true: Mehrotra's predictor-corrector.
                                                     
                                                     
IP_PMM(c, A, Q, b, free_variables, tol, max_it,pc, printlevel): 

                      sets the printlevel.
                                                              
                      0: turn off iteration output
                                                              
                      1: print primal and dual residual and duality measure
                                                              
                      2: print centering parameter and step length
                                                              
OUTPUT: [x,y,z,opt,iter], where:

                     x: primal solution
         
                     y: Lagrange multiplier vector
         
                     z: dual slack variables
         
                     opt: true if problem was solved to optimality, false if problem not solved or found infeasible.
         
                     iter: numeber of iterations to termination.
      
Author: Spyridon Pougkakiotis.

