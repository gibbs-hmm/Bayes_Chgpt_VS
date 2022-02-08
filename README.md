# Bayes_Chgpt_VS

"The Bayesian Change Point and Variable Selection Algorithm: Application to the d18O Record of the Plio-Pleistocene"
By E. Ruggieri and C.E. Lawrence

https://www.tandfonline.com/doi/full/10.1080/10618600.2012.707852

%************************************************************************/
%* The Bayesian Change Point and Variable Selection algorithm 		      */
%* - A program to caluclate the posterior probability of a change point */
%* in a time series and the set of predictor variables that best fits   */
%* each regime of the data set.					                                */
%*                                                                      */
%* Please acknowledge the program author on any publication of          */
%* scientific results based in part on use of the program and           */
%* cite the following article in which the program was described.       */
%*                                                                      */
%* E. Ruggieri and C.E. Lawrence.  The Bayesian Change Point and 	      */
%* Variable Selection Algorithm: Application to the d18O Record of 	    */
%* the Plio-Pleistocene, Journal of Computational and Graphical 	      */
%* Statistics							                                              */
%*                                                                      */
%* Program Author: Eric Ruggieri                                        */
%* Duquesne University                                                  */
%* Pittsburgh, PA 15282                                                 */
%* Email:  ruggierie@duq.edu                                            */
%*                                                                      */
%* Copyright (C) 2012  Duquesne University                              */
%*                                                                      */
%* The Bayesian Change Point and Variable Selection algorithn is free 	*/
%* software: you can redistribute it and/or modify it under the terms 	*/
%* of the GNU General Public License as published by the Free Software  */
%* Foundation, either version 3 of the License, or (at your option) 	  */
%* any later version.							                                      */
%*                                                                      */
%* The Bayesian Change Point and Variable Selection algorithm is 	      */
%* distributed in the hope that it will be useful, but WITHOUT ANY 	    */
%* WARRANTY; without even the implied warranty of MERCHANTABILITY or 	  */
%* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public 	      */
%* License for more details.					                                  */
%*                                                                      */
%* You should have received a copy of the GNU General Public License    */
%* along with Bayesian Change Point and Variable Selection. If not, see */
%* <http://www.gnu.org/licenses/> or write to                           */
%* Free Software Foundation, Inc.                                       */
%* 51 Franklin Street Fifth Floor                                       */
%* Boston, MA 02110-1301, USA.                                          */
%************************************************************************/


Contents:

A) Matlab Code:
   i) Bayes_Chgpt_VS_Script.m
  ii) DP_Exact.m
 iii) DP_d18O.m
  iv) partition_fn.m
   v) pick_k1.m
  vi) logsumlog.m
 vii) sample_model_given_chgpt_sim.m
viii) sample_model_given_chgpt_d18O.m

B) Data Sets
  i) Bayes_Simulation_Short.mat
 ii) Bayes_Simulation.mat
iii) d18O.mat
 iv) data5000.txt

Note: Each of the Matlab (.m) files contain detailed comments / step-by-step guidance



Description of Files:
A) Matlab Code:

i) Bayes_Chgpt_VS_Script.m
This is the main script for the Bayesian Change Point and Variable Selection Algoritm.  All other functions (.m files) are called from here.  
Outline:
	a) Load the data set to be used.
	b) Define the parameters
	c) Step 1: Calculate the Probability Density for each Possible Substring of the Data f(Y_i:j) 
	d) Step 2: Forward Recursion
	e) Step 3: Stochastic Backtrace
		Step 3.1: Sample a Number of Change Points
		Step 3.2: Recursively Sample the Locations of the Change Points
		Step 3.3: Sample a Sub-Model for the Interval Between adjacent Change Points
		Step 3.4: Sample the Regression Parameters for the Interval Between Adjacent Change Points
	f) Graphs of Results

Note: The data sets described in this manuscript utilized parallel processing on a computer cluster to reduce the compute time.  When run on a single processor,
the simulation [which has 1000 data points and 10 predictor variables] takes ~12 hours and 
the analysis of the d18O proxy record [which has 2115 data points and 10 predictor variables] takes ~24 hours to finish.  

Given a fixed set of variables (i.e. no variable selection), the Bayesian Change Point algorithm takes 10-12 minutes to analyze the d18O proxy record.  
When variable selection is added to the problem, the compute time will depend on the number of possible predictor variables.  For example, a model with 10 
variables implies that there are 2^10 possible subsets of variables to consider in each substring.  In essence, this forces you to run the Bayesian 
Change Point algorithm 2^10 times.  For this reason, it is not recommended that you attempt problems of this size on your personal computer, simply 
out of convenience.  A smaller simulation, not included in the manuscript that can quickly be run on a personal computer has been included in this .zip file.  
It contains m=4 regressors for N=250 data points and runs in ~30 seconds.

Please see the end of this README file for a complete list of variables outputted by the Bayesian Change Point and Variable Selection algorithm.

ii) DP_Exact.m
This function performs variable selection via the EBIR algorithm (Ruggieri and Lawrence, 2012) for every possible substring of the data, Y_i:j and
returns the density of the data, marginalizing over all possible sub-models.
Input: Parameters defind by the user, matrix of Regressors (X), Data to be analyzed (Y)
Output: f(Y_i:j)

iii) DP_d18O.m
Similar to DP_Exact.m except written specifically for the analysis of the d18O proxy record in which there were restrictions on varaibles being included
together in a single model.

iv) partition_fn.m
Performs the Forward Recursion step of the Algorithm.
Input: The density of every substring of the data [f(Y_i:j)]
Output: P_k(Y_1:j) - Probability Density of the data given k change points
 
v) pick_k1.m
A function that takes a probability vector as input and draws a random sample from that vector.  Used to sample both the number of change points and
their location within the time series.
Input: A probability vector
Output: Random Sample from that vector

vi) logsumlog.m
A function that calculates the sum of a set of numbers in log form and stores the answer as a logarithm.  Since Matlab would exponentiate before adding 
and then taking a logarithm [log(exp(x) + exp(y))] this brief function provides more numerical precision.
Input: A vector of logarithms
Output: The sum of all numbers in the vector stored as a logarithm
 
vii) sample_model_given_chgpt_sim.m
In order to do Step 3.4: Sample the Parameters of the Regression Model between Successive Change Points, we first need to draw a sample of the sub-model
to use from the posterior distribution.  This function performs that task, utilizing the EBIR algorithm as in DP_Exact.
Input: Location of two successive change points, matrix of Regressors (X), Data to be analyzed (Y)
Output: A sub-model sampled from the posterior distribution

viii) sample_model_given_chgpt_d18O.m
Similar to sample_model_given_chgpt_sim.m except written specifically for the analysis of the d18O proxy record in which there were restrictions on 
varaibles being included together in a single model.



B) Data Sets

Each of the .mat files contains the variables:
Y = The data set to be analyzed
N = Total number of data points
X = The matrix of regressors
m = Number of regressors included in the model.
x = Time points.  For the simulations time points are equally spaced (Ex 1:1000), but for the d18O proxy record, the distance between successive
	data points changes trough time --> Data are spaced every 1, 2, 2.5, and 5 kyr

The two simulations also contain:

num_chgpts = The true number of change points in the data set
true_chgpt = The true locations of the change points in the data set

So that you can compare the output to the true model


i) Bayes_Simulation_Short.mat
This data set was not included in the manuscript but can easily be run on a personal computer and so is included here.
N=250; num_terms = 4; x = 1:250; num_chgpts=1; chgpt_loc=[174];
X = [sin(2*pi/25) sin(2*pi/50) sin(2*pi/75) sin(2*pi/125)];
Y(t) = 	 0.96*X2(t) 	1   <= t <  175	
 	-0.89*X4(t) 	175 <= t <= 250
Noise is added to this data set in Bayes_Chgpt_VS_Script.m

ii) Bayes_Simulation.mat
The simulation described in the manuscript.
N = 1000; num_terms=10; x = 1:1000; num_chgpts=4; chgpt_loc = [308 504 647 749]
X = [sin(2*pi/20) sin(2*pi/30) sin(2*pi/40) sin(2*pi/50) sin(2*pi/60) sin(2*pi/70) sin(2*pi/80) sin(2*pi/100) sin(2*pi/150) sin(2*pi/200)];

Y(t) = 	 3.3115*X2(t)  -0.2222*X7(t) +0.7684*X8(t)			1   <= t < 309
	-0.1844*X2(t)  -0.6663*X3(t) +0.5550*X4(t)			309 <= t < 505
	 0.5957*X4(t)  -1.9639*X6(t) -1.3236*X10(t)			505 <= t < 648
	 0.0468*X1(t)  +2.2412*X5(t) -1.3785*X7(t)  +1.5937*X10(t)	648 <= t < 751
	 2.3593*X3(t)  -1.1507*X5(t) -1.3162*X7(t)  -0.3304*X9(t)	751 <= t < 1000
Noise is added to this data set in Bayes_Chgpt_VS_Script.m

iii) d18O.mat
The d18O record of the Plio-Pleistocene as described in Lisiecki and Raymo (2005) after removing the long-term cooling trend via an exponential function
N=2115; num_terms=10; x varies from every 1 kyr --> 2 kyr --> 2.5 kyr --> 5 kyr as you go back in time.
Y = d18O record.
X = [cos(2*pi/23) cos(2*pi/41) cos(2*pi/53) cos(2*pi/82) cos(2*pi/92) cos(2*pi/95) cos(2*pi/100) cos(2*pi/115) cos(2*pi/123.5) cos(2*pi/404) ... 
     sin(2*pi/23) sin(2*pi/41) sin(2*pi/53) sin(2*pi/82) sin(2*pi/92) sin(2*pi/95) sin(2*pi/100) sin(2*pi/115) sin(2*pi/123.5) sin(2*pi/404) constant];
Note:
Both sine and cosine functions are necessary because there is no guarantee that the proxy record is 'in phase' with time=0.  The sine and cosine 
function allow one to take changes in 'phase' into account.


iv) data5000.txt
The original d18O Proxy Record of the Plio-Pleistocene as found in Lisiecki and Raymo (2005).
The first column contains the time points
The second column contains the proxy values
The third column contains uncertainty estimates.   



Output varibles for the Bayesian Change Point and Variable Selection algorithm: 
N = the total number of data points.
m = the total number of predictor variables.
X = Nxm matrix of the predictor variables.
Y = Nx1 vector of the output variable.
x = Nx1 vector that indexes the data points.  x generally equals 1:N unless the data points are not equally spaced. 

noise = noise level (variance) added to the simulated data.
true_chgpt = The true locations of the change points in the data set.
num_chgpts = true number of change points in the simulation.  

k_i, k_e = parameters for the multivariate normal prior distribtuion on the regression coefficients, beta.
p_i, p_e = probability of including or excluding a predictor variable, respectively
sig_0, v_0 = parameters for the Scaled Inverse-Chi^2 distribution on the error variance, sigma^2.
k_max = maximum number of allowed change points.
min_chgpt = minimun distance between adjacent change points.
parameters = [min_chgpt, v_0, sig_0, k_i, k_e, p_i, p_e].  Parameters used in the likelihood calculation

Py = NxN matrix containing the probability density of the data for every possible subset of the data, i.e. f(Y_i:j), Equation (3).
P = P_k(Y_i:j) = (k_max) x (N) matrix containing the probabiity density of the first j observations of the data containing k change points, Equations (4) and (5).
num_comb = N_k in the manuscript. (k_max)x(N) matrix containing the number of permissible change point solutions with k change points in the first j data points.
k = posterior probability of selecting k change points. (k_max)x1 vector, f(k|Y), Equation (7).  Compare to num_chgpts, which is the true number of change points in the simulation.

num_samp = number of independent solutions sampled directly from the posterior distribution.

chgpt_loc = 1xN vector indicating the marginal probability of a change point being selected at a given data point.  Compate to true_chgpt which are their true locations in the simulation data sets.
chgpt_soln = (num_samp)x(k_max) matrix containing the individual change point sampled solutions.
samp_holder = Nxm matrix containing the number of times that each variable was selected for inclusion at a given data point.
BETA = mxN matrix containing the average regression coefficients inferred at each data point.
model = Nx1 vector containing the average inferred model.
centroid = Nxm matrix containing the centroid solution to the variable selection problem.  A variable is included in the centroid solution at each data point if it was selected for inclusion in more than half of the sampled solutions.

R_2 = R^2 value for the inferred model.
