# WavesHelmholtz

Here are 4 examples with  solvers for Helmholtz equation  in 2D:

**********************************************************

1) First solver does't uses computation of residual:

For compilation use:

make -f makefile_v2 runmaxwellv2

To run:

>runmaxwellv2 argv[1]  argv[2]

argv[1] - number of preconditioner (1 - Jacobi, 2 - Gauss-Seidel, 3 - SOR )

argv[2]  - number of nodes in x and y directions 


% number of nodes for meshes in convergence tests:
%  l= 6; h= 0.015625 = 1/2^6  ---> nx = ny = 1/h +1 = 65
%  l= 5; h =    0.03125 = 1/2^5 ---> nx = ny = 1/h + 1 = 33
%  l= 4; h =  0.0625 = 1/2^4   ---> nx = ny = 1/h + 1 = 17
%  l= 3 ; h = 0.125 = 1/(2^3) ----> nx = ny = 1/h + 1 = 9
  
h=[0.015625, 0.03125, 0.0625, 0.125];
The C++/PETSc code  computes finite difference approximation for the test problem

Δu(x,y)+ omega^2*epsilon*u(x,y) = f(x,y)  in (0, 1) x (0, 1) =: D,
		
	       u = 0  on dD,

with

f(x, y) = -(8*pi^2)*sin(2*pi*x)*sin(2*pi*y) - 2i*x*(1-x) - 2i*y*(1-y)
	  +omega*omega*epsilon(x)*(sin(2*pi*x)*sin(2*pi*y) +
	  i*x*(1 - x)*y*(1 - y))


with the exact solution

u(x, y) = sin(2*pi*x)*sin(2*pi*y) + i*x*(1 - x)*y*(1 - y).

and epsilon(x,y) =  2*exp(-((x-x_0)*(x-x_0)/(2*c_x*c_x) +(y-y_0)*(y-y_0)/(2*c_y*c_y)))

To compile:
make  -f makefile_v2 runmaxwellv2

To run solver:

./runmaxwellv2  argv[1]  argv[2]


argv[1] - number of preconditioner (1 - Jacobi, 2 - Gauss-Seidel, 3 - SOR )

argv[2]  - number of nodes in x and y directions 


where N is the number of nodes on each side for the finite difference grid.

Output: nodes.m, values.m, containing a list of the nodes, and corresponding real and imaginary parts of the approximation, respectively.

To visualize:
matlab -nodesktop -r viewer
(starts Matlab in terminal and runs the visualization script).

**********************************************************************************************************
2) 

To run solver A u = b without preconditioner but with computing residual at every KSP iteration: compile the program cplmaxwell.cpp  as

>make runmaxwell

Then to run use:

>runmaxwell argv[1]

Here, argv[1] is number of nodes in x and y directions

This program writes residual for every ksp iteration in the file residual.m

To do it we use the PETSC-command

 KSPGetResidualHistory(ksp, &hist, &nHist);

Gets the array used to save the residual history and the number of residuals

***************************************************************************************************************

3) Testing  Helmholtz solver for different frequencies. The code can recognize if the inout frequency is resonant frequency.


For compilation use:

>make -f makefile_v3 runmaxwellv3

To run:

>runmaxwellv3 argv[1]  argv[2]


argv[1]  - number of nodes in x and y directions 

argv[2] - frequency omega

4) Testing convergence of computed solution:

   To compile:
  >make -f makefileconv_v2  runconv2

   To run:

>runconv2 argv[1]  argv[2]


argv[1] - number of preconditioner (1 - Jacobi, 2 - Gauss-Seidel, 3 - SOR )

argv[2]  - number of nodes in x and y directions 

Output files for analyzing of results:
values.m  ---> contains computed solution
errorvalues.m ----> contains discrete error between exact and computed solution

Both files can be used in Matlab in order to compute convergence rate.




   
