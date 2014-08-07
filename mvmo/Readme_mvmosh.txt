------------------------------------------------------------------------
Hybrid variant of Mean Variance Mapping Optimization algorithm (MVMO-SH)
------------------------------------------------------------------------

Matlab Demo Version July 2013 

Algorithmic implementation and maintenance by Dr. José L. Rueda 
<jose.rueda@uni-duisburg-essen.de>
Demo's conceptual design and implementation by Sebastian Wildenhues, M.Sc. 
<sebastian.wildenhues@uni-due.de>

© 2013 Institute of Electrical Power Systems, University Duisburg-Essen
See http://www.uni-due.de/ean/mvmo/ for more info.

------------
Introduction
------------

The Mean-variance mapping optimization (MVMO), conceived by Prof. I. Erlich 
is a new optimization algorithm, whose conceptual framework has certain 
similarities to other heuristic approaches. 
Its novel distinguishing feature, however, is its use of a special mapping 
function applied for mutating the offspring on the basis of mean and variance 
of the set comprising of the n-best solutions attained so far, and saved in 
the archive. Based on simple mathematical relationships, the shape of the 
mapping curve is adjusted according to the progress of the search process. 
Since the mapping function is defined in the interval [0,1] for all 
optimization variables, original variables have to be rescaled to within 
this range. However, the fitness evaluation is performed using the actual 
values in the problem space (i.e. conversion is embedded in this task). 
The previous implementation of MVMO represents a single particle approach 
based on random search and its own experiences in the search process.

The hybrid variant of MVMO (MVMO-SH) is a swarm intelligence based procedure 
which incorporates local search and multi-parent crossover strategies to 
increase the search diversity while striving for a balance between exploration 
and exploitation.

--------------
Key references
--------------

[1] I. Erlich, G. K. Venayagamoorthy, and W. Nakawiro, "A mean-variance 
    optimization algorithm," in Proc. 2010 IEEE World Congress on Computational 
    Intelligence, Barcelona, Spain.
[2] J.L. Rueda, and I. Erlich, “Evaluation of the Mean-Variance Mapping Optimization
    for Solving Multimodal Problems,” in Proc. 2013 IEEE Symposium Series on 
    Computational Intelligence, Singapore, April 2013.
[3] J.L. Rueda, and I. Erlich, “Mean-Variance Mapping Optimization for Solving 
    the IEEE-CEC 2013 Competition Problems,” in Proc. 2013 IEEE Congress on 
    Evolutionary Computation, June 2013.

------------
Terms of use
------------

- MVMO-SH is distributed free of charge in the hope that it will be useful, 
  but WITHOUT ANY WARRANTY. 
- The original version or modified versions of MVMO-SH may not be redistributed 
  or distributed without prior written permission.
- It is requested that publications derived from the use of MVMOS explicitly 
  acknowledge that fact by citing [1]-[3].

-------------------
System requirements
-------------------

MATLAB(R) 7.9.0.529 or later
MATLAB(R) Optimization Toolbox (R) (if local search strategy is requested)

----------------
Running the demo
----------------
To explore the algorithm in an intuitive and accessible manner, open the demo by 
pressing F5 (demo source is active window) or by typing MVMOSH_demo_new at the prompt 
(file is in Matlab search path). On the upper-left side of the starting figure, 
the popup menu can be used to choose between various recognized continuous 
and discontinuous unconstrained minimization problems.

If 'Start Evaluation' is clicked, the procedure will run the respective problem 
with a set of predefined parameters, such as the 'number of particles' equal to 1
(previous implementation of MVMO). Only two optimization variables (X1,X2) are 
specified then, which results intermediate candidate solution samples to be plotted 
spatially (fitness = f(X1,X2)) while the mapping functions' shape associated with 
(X1,X2) from the (best) particle is illustrated at the same time. 
Furthermore, the figure offers two sections to observe 
1) the actual best archive, containing the n best-so-far solutions of the 
   global best particle
2) the search history in terms of global best fitness and best (X1,X2)
in the bottom part.

In contrast to directly running an optimization problem, the starting figure
offers two sections to more specifically adjust the procedural and algorithmic 
behaviour as a preliminary step:

Procedure Settings
-  Maximum number of function evaluations	
   <= 1000 for two optimization variables (illustrative purposes) 
   <= 100000 for up to 30 optimization variables (testing, benchmarking)
-  Print step size
   Determines how often intermediate results are drawn spatially 
   (2 optimization variables) and/or printed to the MATLAB command window 
   (up to 30 optimization variables).
-  Probability of local search (percentage/number of optimization variables)
   If set equal to 0, memory of particles within the swarm is used as the only 
   searching engine. 
   If set > 0, the local search strategy is employed in a complementary manner 
   but requires MATLAB(R) Optimization Toolbox (R) to be installed in order 
   to perform embedded nonlinear optimization.   
-  Number of optimization variables (2 up to 30)

Algorithm Settings
-  Particle setting [Number, Initial percentage 'good', Final percentage 'good']
   If Number is set to 1, the previous implementation of MVMO is used. 
   Otherwise, multple particles perform concurrent search within the swarm. 
   Note that the swarm size must be chosen appropriately to enable multi-parent 
   crossover strategies over the search progress. Any irrelevant setting 
   (in conjunction with the Initial percentage 'good' and Initial percentage 'bad') 
   will be intervened by the demo.
-  Archive size
   Number n of best-so-far solutions representing each particle's memory
-  Selection mode
   1 - Random selection   
   2 - Neighbor group - block stepping selection
   3 - Neighbor group - single stepping selection
   4 - Sequential-random selection
-  No. of variables selected for mutation [Initial(<total), Final (<=initial)]
-  Shape-scaling factor [fs_factor_start(<=1), fs_factor_end(>=1)]
   Initial shape-scaling factor fs_factor_start, which is quadratically increased 
   with iterations towards the final shape-scaling factor fs_factor_end.
-  [dddd(1-5), delta_dddd_start(0.01-0.4), delta_dddd_end(<delta_dddd_start)]
   Control values for exploiting the asymmetric characteristic of the mapping 
   function as well as for zero variance handling.
-  Decrement strategy of variables selected for mutation [1,2,3,41,42,43]
   1 - Linear
   2 - Quadratic progressive
   3 - Quadratic degressive
   plus 40 means variable randomly changed within the range
   This option takes effect only if the initial no. of variables selected for 
   mutation does not equal the final no. of variables selected for mutation.

For each chosen number of optimization variables, after running the problem, 
the global best fitness is plotted in a separate figure. The demo can be restartet 
using the refreshed set of parameters by clicking 'New Problem' in the bottom-
left part of that figure.
