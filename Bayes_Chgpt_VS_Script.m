%************************************************************************/
%* The Bayesian Change Point and Variable Selection algorithm           */
%* - A program to caluclate the posterior probability of a change point */
%* in a time series and the set of predictor variables that best fits   */
%* each regime of the data set.                                         */
%*                                                                      */
%* Please acknowledge the program author on any publication of          */
%* scientific results based in part on use of the program and           */
%* cite the following article in which the program was described.       */
%*                                                                      */
%* E. Ruggieri and C.E. Lawrence.  The Bayesian Change Point and        */
%* Variable Selection Algorithm: Application to the d18O Record of      */
%* the Plio-Pleistocene, Journal of Computational and Graphical         */
%* Statistics                                                           */
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
%* Foundation, either version 3 of the License, or (at your option) 	*/
%* any later version.                                                   */
%*                                                                      */
%* The Bayesian Change Point and Variable Selection algorithm is        */
%* distributed in the hope that it will be useful, but WITHOUT ANY      */
%* WARRANTY; without even the implied warranty of MERCHANTABILITY or 	*/
%* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public        */
%* License for more details.                                            */
%*                                                                      */
%* You should have received a copy of the GNU General Public License    */
%* along with Bayesian Change Point and Variable Selection. If not, see */
%* <http://www.gnu.org/licenses/> or write to                           */
%* Free Software Foundation, Inc.                                       */
%* 51 Franklin Street Fifth Floor                                       */
%* Boston, MA 02110-1301, USA.                                          */
%************************************************************************/
% This script will perform the Bayesian Change Point and Variable Selection
% algorithm for simulated data and with modification (see below) to the d18O
% record of the Plio-Pleistocene

% Outline of the Bayesian Change Point and Variable Selection algorithm:
% 1) Load the Data  
% 2) Define the Parameter Values
% 3) Calculate the Probability Density of the Data for Each Sub-Interval
% 4) Forward Recursion [Dynamic Programming]
% 5) Stochastic Backtrace via Bayes Rule
%       a) Sample a Number of Change Points
%       b) Sample the Location of the Change Points
%       c) Sample a Sub-Model for the Interval Between adjacent Change Points
%       c) Sample the Regression Parameters Between Adjacent Change Points
% 6) Plot the Results
%

clear;
%************(1)Load the Data*****************************
load Bayes_Simulation_Short
%load Bayes_Simulation
%load d18O

%{
Y = Nx1 column vector of data points
X = Nxm matrix of predictor variables
x = Nx1 column vector of time points.  X gives the values of the regressors at each value of x
m = number of predictor variables
N = number of data points
num_chgpt = true number of change points in the simulated data (Simulations Only)
chgpt_loc = true location of the change points in the data set (Simulations Only)
%}

%*************(2)Define the Parameter Values**************
noise=1.0;      %The level (variance) of noise to be added to the simulated data
                %Note: Only add noise to the simulations!
Q=sqrt(noise)*randn(N,1);
Y=Y+Q;

k_max=8;                % Maximum number of allowed change points
min_chgpt=25;           % Minimum distance between adjacent change points
v_0=5; sig_0= var(Y);   % Prior parameters for the inverse chi-square distribution on the error variance, sigma^2
k_i=0.01;  k_e = 100;   % Scale parameter for the milti-variate normal prior on the regression coefficients, beta
p_i=0.5;                % Probability of including a predictor variable
p_e=0.5;                % Probability of including a predictor variable
parameters = [min_chgpt, v_0, sig_0, k_i, k_e, p_i, p_e]; 
                        % Parameters used in the likelihood calculation

%******(3)Calculate the Probability Density of Data for Each Sub-Interval*********
%STEP 1: CALCULATING THE PROBABILITY DENSITY OF THE DATA, f(Y_i:j)   

%Py is the matrix that will hold the log-probability of each substring of the data, Y_i:j
for i=1:50:N
    % Groups of 50 are sent in to allow for parallel processing or to
    % provide updates on the progress of the algorithm
	Py1=DP_exact(i, i+49, x,X,Y,N, parameters);     %To be used with SIMULATIONS
    %Py1=DP_d18O(i, i+49, x,X,Y,N, parameters);     %To be used with d18O record
    if (i==1)
        Py=Py1;
    else
        Py=[Py;Py1];
    end
end

clear Py1;
disp('Probabilities Calculated');

%****************(4) Forward Recursion *********************************
%STEP 2: FORWARD RECURSION

[P num_comb]=partition_fn(Py, k_max, N);    
% num_comb is combinatorial factor used to obtain normalization constant, Equation (6)
% P = P_k(Y_1:j), Equations (4) and (5)
disp('Partition Function Built');

%****************(5) Stochastic Backtrace via Bayes Rule****************
%STEP 3: STOCHASTIC BACKTRACE
%NOTE: To sample from the d18O record, use commented code at bottom of the loop

k=[Py(1,N); P(:,N)];        % P(k|y)
    k(1) = k(1) + log(0.5); % P(K = 0)
for i=2:k_max+1  
    k(i) = k(i) + log(0.5) - log(k_max) - log(num_comb(i-1,N)); % P(K = k, c1, ..., ck)
end

total_k=logsumlog(k);        % logsumlog adds logarithms and puts answer in log form while keeping precision - Equation (6)
k(:)=exp(k(:)-total_k);      % Normalize the vector, f(K=k | y), Equation (7)

save BAYES_SIMULATED_output

%The above calculations need only be done ONCE.  After they are completed,
%the sampling procedure can be done as many times as desired without having
%to run the above code.
%You can draw as many samples as you want and output the figures and
%quantities of your choice.  A few possibilities are shown below:

% Variables used in sampling procedure:
num_samp=500;                   % Number of independent samples drawn directly from the posterior distribution
model=zeros(N,1);               % Holds the final model-averaged result
samp_holder=zeros(N,m);         %The model selected at each data point.
    %Can be used to find marginal probability of including a variable at a given point
chgpt_soln=zeros(num_samp, k_max);  % Contains each of the num_samp change point solutions
chgpt_loc=zeros(1,N);               % The locations/marginal probabilities of a change point
BETA = zeros(m,N);              % Holds the regression coefficients

for i=1:num_samp
    %****** Sampling from the Simulated Data Sets ******
    %STEP 3.1 - SAMPLE A NUMBER OF CHANGE POINTS
    num_chgpts = pick_k1(k)-1;  % Since we allow for 0 changepoints, function returns the index of the 'k' vector,
                                % which is offset from the number of change points by 1
    if(num_chgpts>0)
        %STEP 3.2 - RECURSIVELY SAMPLE THE LOCATIONS OF THE CHANGE POINTS
        back=N;
        for kk=num_chgpts:-1:2      % Start at the end of the time series and work backwards
            temp=zeros(1,back-1);   
            for v=1:back-1          % Build the vector to sample from
                temp(v)= P(kk-1,v)+Py(v+1,back);
            end                         % Equation (8)
            total=logsumlog(temp);      % Caculate the normalization constant
            temp(:)=exp(temp(:)-total); % Normalize the vector
            changepoint=pick_k1(temp);  % Sample the location of the change point
            chgpt_soln(i,kk)=changepoint;                      % Keep track of individual change point solutions
            chgpt_loc(changepoint)= chgpt_loc(changepoint) +1; % Keep track of change point locations
            
            %STEP 3.3- SAMPLE A SUB-MODEL FOR THE INTERVAL BETWEEN ADJACENT CHANGE POINTS
            %Samples a model from the posterior distribution, f(Am | Y), Equation (A7)
            
            %To be used with SIMULATIONS
            Am = sample_model_given_chgpt_sim(changepoint+1,back,X,Y,parameters);
            for j=changepoint+1:back
                samp_holder(j,:) = samp_holder(j,:)+Am;    %Keep track of model selected
            end
            
            %STEP 3.4- SAMPLE THE REGRESSION PARAMETERS FOR THE INTERVAL BETWEEN ADJACENT CHANGE POINTS
            % Regression Coefficients (Beta)
            I_Am=eye(m);
            for ii=1:m
                if (Am(ii)==1)
                    I_Am(ii,ii)=0.01;
                else
                    I_Am(ii,ii)=100;
                end
            end
           
            XTy=X(changepoint+1:back,:)'*Y(changepoint+1:back); 
            J=I_Am+X(changepoint+1:back,:)'*X(changepoint+1:back,:); J=inv(J);
            
            beta_hat=J*(XTy);
            for j=1:m
               BETA(j,changepoint+1:back) = BETA(j,changepoint+1:back) +beta_hat(j)*ones(1,back-changepoint);  
            end
            
            % Update the model averaged result
            model(changepoint+1:back)=model(changepoint+1:back)+X(changepoint+1:back,:)*beta_hat;
            
            back=changepoint;   % Now work with the smaller segment
        end
        
        % The Final Change Point
        %STEP 3.2 - SAMPLE THE LOCATIONS OF THE CHANGE POINTS
        kk=1;
        temp=zeros(1,back-1);
        for v=1:back-1
            temp(v)= Py(1,v)+Py(v+1,back);  % Build the vector to sample from
        end                                 % Equation (8)
        total=logsumlog(temp);              % Normalization constant
        temp(:)=exp(temp(:)-total);         % Normalize the vector
        changepoint=pick_k1(temp);          % Sample a change point location
        chgpt_soln(i,kk)=changepoint;                      % Keep track of individual change point solutions
        chgpt_loc(changepoint)= chgpt_loc(changepoint) +1; % Keep track of change point locations
        
        %STEP 3.3- SAMPLE A SUB-MODEL FOR THE INTERVAL BETWEEN ADJACENT CHANGE POINTS
        %Samples a model from the posterior distribution, f(Am | Y), Equation (A7)
        
        %To be used with SIMULATIONS
        Am = sample_model_given_chgpt_sim(changepoint+1,back,X,Y, parameters);
        for j=changepoint+1:back
            samp_holder(j,:) = samp_holder(j,:)+Am;    %Keep track of model selected
        end
        
        %STEP 3.4- SAMPLE THE REGRESSION PARAMETERS FOR THE INTERVAL BETWEEN ADJACENT CHANGE POINTS
        % Regression Coefficients (Beta)
        I_Am=eye(m);
        for ii=1:m
            if (Am(ii)==1)
                I_Am(ii,ii)=0.01;
            else
                I_Am(ii,ii)=100;
            end
        end
        XTy=X(changepoint+1:back,:)'*Y(changepoint+1:back); 
        J=I_Am+X(changepoint+1:back,:)'*X(changepoint+1:back,:); J=inv(J);
        
        beta_hat=J*(XTy);
        for j=1:m
               BETA(j,changepoint+1:back) = BETA(j,changepoint+1:back) +beta_hat(j)*ones(1,back-changepoint);  
        end
        
        %Update the model averaged result
        model(changepoint+1:back)=model(changepoint+1:back)+X(changepoint+1:back,:)*beta_hat;
        
    else        % 0 change points, so a single homogeneous segment
        changepoint=N;
    end
    
    %The final sub-interval
    %STEP 3.3- SAMPLE A SUB-MODEL FOR THE INTERVAL BETWEEN ADJACENT CHANGE POINTS
    %Samples a model from the posterior distribution, f(Am | Y), Equation (A7)
    
    %To be used with SIMULATIONS
    Am = sample_model_given_chgpt_sim(1,changepoint,X,Y,parameters);
    for j=1:changepoint
        samp_holder(j,:) = samp_holder(j,:)+Am;    %Keep track of model selected
    end
    
    %STEP 3.4- SAMPLE THE REGRESSION PARAMETERS FOR THE INTERVAL BETWEEN ADJACENT CHANGE POINTS
    % Regression Coefficients (Beta)
    I_Am=eye(m);
    for ii=1:m
        if (Am(ii)==1)
            I_Am(ii,ii)=0.01;
        else
            I_Am(ii,ii)=100;
        end
    end
    XTy=X(1:changepoint,:)'*Y(1:changepoint);
    J=I_Am+X(1:changepoint,:)'*X(1:changepoint,:); J=inv(J);
    
    beta_hat=J*(XTy);
    for j=1:m
        BETA(j,1:changepoint) = BETA(j,1:changepoint) +beta_hat(j)*ones(1,changepoint);
    end
    
    %update the model averaged result
    model(1:changepoint)=model(1:changepoint)+X(1:changepoint,:)*beta_hat;
    
    %{
    %*** Sampling from the d18O Record*******
    %STEP 3.1 - SAMPLE A NUMBER OF CHANGE POINTS
    num_chgpts = pick_k1(k)-1;  % Since we allow for 0 changepoints, function returns the index of the 'k' vector,
                                % which is offset from the number of change points by 1
    if(num_chgpts>0)
        %STEP 3.2 - RECURSIVELY SAMPLE THE LOCATIONS OF THE CHANGE POINTS
        back=N;
        for kk=num_chgpts:-1:2      % Start at the end of the time series and work backwards
            temp=zeros(1,back-1);   
            for v=1:back-1          % Build the vector to sample from
                temp(v)= P(kk-1,v)+Py(v+1,back);
            end                         % Equation (8)
            total=logsumlog(temp);      % Caculate the normalization constant
            temp(:)=exp(temp(:)-total); % Normalize the vector
            changepoint=pick_k1(temp);  % Sample the location of the change point
            chgpt_soln(i,kk)=changepoint;                      % Keep track of individual change point solutions
            chgpt_loc(changepoint)= chgpt_loc(changepoint) +1; % Keep track of change point locations
            
            %STEP 3.3- SAMPLE A SUB-MODEL FOR THE INTERVAL BETWEEN ADJACENT CHANGE POINTS
            %Samples a model from the posterior distribution, f(Am | Y), Equation (A7)
            
            % To be used with d18O record
            Am = sample_model_given_chgpt_d18O(changepoint+1,back,X,Y,parameters);
            for j=changepoint+1:back
                samp_holder(j,:) = samp_holder(j,:)+Am;    %Keep track of model selected
            end
            Am = [Am Am 1];        % Accounts for both sine and cosine function for each predictor variable
            
            %STEP 3.4- SAMPLE THE REGRESSION PARAMETERS FOR THE INTERVAL BETWEEN ADJACENT CHANGE POINTS
            % Regression Coefficients (Beta)
            I_Am=eye(2*m+1);
            for ii=1:(2*m+1)
                if (Am(ii)==1)
                    I_Am(ii,ii)=0.01;
                else
                    I_Am(ii,ii)=100;
                end
            end
           
            XTy=X(changepoint+1:back,:)'*Y(changepoint+1:back); 
            J=I_Am+X(changepoint+1:back,:)'*X(changepoint+1:back,:); J=inv(J);
            
            beta_hat=J*(XTy);
            for j=1:m
               BETA(j,changepoint+1:back) = BETA(j,changepoint+1:back) +beta_hat(j)*ones(1,back-changepoint);  
            end
            
            % Update the model averaged result
            model(changepoint+1:back)=model(changepoint+1:back)+X(changepoint+1:back,:)*beta_hat;
            
            back=changepoint;   % Now work with the smaller segment
        end
        
        % The Final Change Point
        %STEP 3.2 - SAMPLE THE LOCATIONS OF THE CHANGE POINTS
        kk=1;
        temp=zeros(1,back-1);
        for v=1:back-1
            temp(v)= Py(1,v)+Py(v+1,back);  % Build the vector to sample from
        end                                 % Equation (8)
        total=logsumlog(temp);              % Normalization constant
        temp(:)=exp(temp(:)-total);         % Normalize the vector
        changepoint=pick_k1(temp);          % Sample a change point location
        chgpt_soln(i,kk)=changepoint;                      % Keep track of individual change point solutions
        chgpt_loc(changepoint)= chgpt_loc(changepoint) +1; % Keep track of change point locations
        
        %STEP 3.3- SAMPLE A SUB-MODEL FOR THE INTERVAL BETWEEN ADJACENT CHANGE POINTS
        %Samples a model from the posterior distribution, f(Am | Y), Equation (A7)
        
        % To be used with d18O record
        Am = sample_model_given_chgpt_d18O(changepoint+1,back,X,Y,parameters);
        for j=changepoint+1:back
            samp_holder(j,:) = samp_holder(j,:)+Am;    %Keep track of model selected
        end
        Am = [Am Am 1];        % Accounts for both sine and cosine function for each predictor variable
        
        %STEP 3.4- SAMPLE THE REGRESSION PARAMETERS FOR THE INTERVAL BETWEEN ADJACENT CHANGE POINTS
        % Regression Coefficients (Beta)
        I_Am=eye(2*m+1);
        for ii=1:(2*m+1)
            if (Am(ii)==1)
                I_Am(ii,ii)=0.01;
            else
                I_Am(ii,ii)=100;
            end
        end
        XTy=X(changepoint+1:back,:)'*Y(changepoint+1:back); 
        J=I_Am+X(changepoint+1:back,:)'*X(changepoint+1:back,:); J=inv(J);
        
        beta_hat=J*(XTy);
        for j=1:m
               BETA(j,changepoint+1:back) = BETA(j,changepoint+1:back) +beta_hat(j)*ones(1,back-changepoint);  
        end
        
        %Update the model averaged result
        model(changepoint+1:back)=model(changepoint+1:back)+X(changepoint+1:back,:)*beta_hat;
        
    else        % 0 change points, so a single homogeneous segment
        changepoint=N;
    end
    
    %The final sub-interval
    %STEP 3.3- SAMPLE A SUB-MODEL FOR THE INTERVAL BETWEEN ADJACENT CHANGE POINTS
    %Samples a model from the posterior distribution, f(Am | Y), Equation (A7)
    
    % To be used with d18O record
    Am = sample_model_given_chgpt_d18O(1,changepoint,X,Y,parameters);
    for j=1:changepoint
        samp_holder(j,:) = samp_holder(j,:)+Am;    %Keep track of model selected
    end
    Am = [Am Am 1];        % Accounts for both sine and cosine function for each predictor variable
    
    %STEP 3.4- SAMPLE THE REGRESSION PARAMETERS FOR THE INTERVAL BETWEEN ADJACENT CHANGE POINTS
    % Regression Coefficients (Beta)
    I_Am=eye(2*m+1);
    for ii=1:(2*m+1)
        if (Am(ii)==1)
            I_Am(ii,ii)=0.01;
        else
            I_Am(ii,ii)=100;
        end
    end
    XTy=X(1:changepoint,:)'*Y(1:changepoint);
    J=I_Am+X(1:changepoint,:)'*X(1:changepoint,:); J=inv(J);
    
    beta_hat=J*(XTy);
    for j=1:m
        BETA(j,1:changepoint) = BETA(j,1:changepoint) +beta_hat(j)*ones(1,changepoint);
    end
    
    %update the model averaged result
    model(1:changepoint)=model(1:changepoint)+X(1:changepoint,:)*beta_hat;
    %}
end

BETA=BETA/num_samp;             % Average regression coefficient at each data point
chgpt_loc=chgpt_loc/num_samp;   % Marginal posterior probability of a change point
model=model/num_samp;           % Average of inferred model at each data point

clear Am I_Am J Q XTy back beta_hat changepoint i ii j kk r rr temp total total_k v

%**********(6) Plot the Results ********************************
% Adapt as Necessary

%Graph of Simulated Data with Inferred Model
h1 = figure('visible', 'off');
subplot(3,1,1);
plot(x(1:N), Y);
title('Input Data'); 
subplot(3,1,2); 
plot(x(1:N), model)
title('Inferred Model');
subplot(3,1,3);
plot(x(1:N), chgpt_loc);
title('Probability of a Change Point')
ylabel('Probability')
xlabel('Data Point')
axis([0 N 0 1]);
print(h1, 'fig1.pdf', '-dpdf');
close(h1);

%Heat Map displaying how many times each variable was selected at a given data point
h2 = figure('visible', 'off'); 
imagesc(samp_holder');    
title('The Number of Times Each Variable was Selected at Each Data Point');
xlabel('Data Point');
ylabel('Variables');
set(gca,'YTick',1:m)
print(h2, 'fig2.pdf', '-dpdf');
close(h2);

%Regression Coefficients - use as many as needed
h3 = figure('visible', 'off');
plot(BETA(1,:), 'b'); hold;
plot(BETA(2,:), 'g');plot(BETA(3,:), 'r');plot(BETA(4,:), 'c');
%plot(BETA(5,:), 'm');plot(BETA(6,:), 'y');plot(BETA(7,:), 'k');
%plot(BETA(8,:), 'b');plot(BETA(9,:), 'g');plot(BETA(10,:), 'r');hold;
title('Amplitudes of Each Forcing Function');
xlabel('Data Point')
print(h3, 'fig3.pdf', '-dpdf');
close(h3);

%Centroid Solution for Variable Selection - select a variable at a given
%data point if it was included in over 50% of the models selected
centroid=zeros(x(N),m);
for i=1:x(N)
    for j=1:m
        if (samp_holder(i,j)<num_samp/2)
            centroid(i,j)=0;
        else
            centroid(i,j)=1;
        end
    end
end
h4 = figure('visible', 'off');
imagesc(centroid'); 
title('Centroid Solution');
xlabel('Data Point');
ylabel('Variables');
set(gca,'YTick',1:m)
print(h4, 'fig4.pdf', '-dpdf');
close(h4);

%Calculate R^2 value
r=0; rr=0;
for i=1:N
    r = r+ (Y(i)-mean(Y))^2;
    rr= rr+ (Y(i)-model(i))^2;
end
R_2=1-rr/r;


%**If data points are not equally spaced, as in the d18O record, this    **
%**needs to be taken into account to avoid distorted graphs              **
%{
real_samp_holder=zeros(x(N),m); x=floor(x);
real_BETA=zeros(m, x(N));

for i=1:N-1
    count=x(i);
    while(count<x(i+1))
        real_samp_holder(count,:)=samp_holder(i,:);
        real_BETA(:,count) = BETA(:,i);
        
        count=count+1;
    end
end

%Heat Map displaying how many times each variable was selected at a given data point
figure(2); imagesc(real_samp_holder')
title('The Number of Times Each Variable was Selected at Each Data Point');
xlabel('Data Point');
ylabel('Variables');
set(gca,'YTick',1:m)

%Regression Coefficients - use as many as needed
figure(3); plot(real_BETA(1,:), 'b'); hold;
plot(real_BETA(2,:), 'g');plot(real_BETA(3,:), 'r');plot(real_BETA(4,:), 'c');
plot(real_BETA(5,:), 'm');plot(real_BETA(6,:), 'y');plot(real_BETA(7,:), 'k');
plot(real_BETA(8,:), 'b');plot(real_BETA(9,:), 'g');plot(real_BETA(10,:), 'r');hold;
title('Amplitudes of Each Forcing Function');
xlabel('Data Point')

%CENTROID Solution for Variable Selection - select a variable at a given
%data point if it was included in over 50% of the models selected
centroid=zeros(x(N),m);
for i=1:x(N)
    for j=1:m
        if (real_samp_holder(i,j)<num_samp/2)
            centroid(i,j)=0;
        else
            centroid(i,j)=1;
        end
    end
end
figure(4); imagesc(centroid'); 
title('Centroid Solution');
xlabel('Data Point');
ylabel('Variables');
set(gca,'YTick',1:m)
%}

f = fopen('output.txt', 'w');
fprintf(f, '%s\n', 'parameters for the multivariate normal prior distribtuion on the regression coefficients, beta.');
fprintf(f, '%f %f\n', k_i, k_e );
fprintf(f, '%s\n', 'probability of including or excluding a predictor variable, respectively');
fprintf(f, '%f %f\n', p_i, p_e );
fprintf(f, '%s\n', 'parameters for the Scaled Inverse-Chi^2 distribution on the error variance, sigma^2.');
fprintf(f, '%f %f\n', sig_0, v_0 );
fprintf(f, '%s\n', 'maximum number of allowed change points.');
fprintf(f, '%f\n', k_max );
fprintf(f, '%s\n', 'minimun distance between adjacent change points.');
fprintf(f, '%f\n', min_chgpt );
fprintf(f, '%s\n', '[min_chgpt, v_0, sig_0, k_i, k_e, p_i, p_e].  Parameters used in the likelihood calculation');
fprintf(f, '%f\n', parameters );
fprintf(f, '%s\n', 'number of independent solutions sampled directly from the posterior distribution.');
fprintf(f, '%f\n', num_samp );
fprintf(f, '%s\n', 'R^2 value for the inferred model.');
fprintf(f,'%f\n', R_2 );
fclose(f);

save('Py.txt', 'Py', '-ascii');
save('P.txt', 'P', '-ascii');
save('num_comb.txt', 'num_comb', '-ascii');
save('k.txt', 'k', '-ascii');
save('chgpt_loc.txt', 'chgpt_loc', '-ascii');
save('chgpt_soln.txt', 'chgpt_soln', '-ascii');
save('samp_holder.txt', 'samp_holder', '-ascii');
save('BETA.txt', 'BETA', '-ascii');
save('model.txt', 'model', '-ascii');
save('centroid.txt', 'centroid', '-ascii');

disp('Script Finished');
