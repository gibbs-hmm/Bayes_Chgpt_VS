function PROB = DP_d18O(start, stop,x,XX,Y,N, parameters)

% Step 1: Calculate the Probability Density of the Data
%
% The EBIR [Exact Bayesian Inference in Regression] algorithm (Ruggieri and 
% Lawrence 2012) calculates the exact probability of all possible
% subsets of predictor variables.  The recursion can be thought of as a tree  
% with the null model at the root.  For the analysis of the d18O record of 
% the Plio-Pleistocene, there are restrictions on which predictor variables
% are permitted to be in the same model.  Thus, the main difference of  
% this function and DP_exact.m is that here there are less than 2^m subsets
% of predictor variables to consider.
%
% Output: f(Y_i:j)

min_chgpt = parameters(1);  % Minimum distance between adjacent change points
v_0 = parameters(2);        % Prior parameter for the error variance, sigma^2
sig_0 = parameters(3);      % Prior parameter for the error variance, sigma^2
k_inc = parameters(4);      % Prior parameter for the regression coefficients, beta (k_i in manuscript)
k_exc = parameters(5);      % Prior parameter for the regression coefficients, beta (k_e in manuscript)
% The tree starts with the null model at its root, i.e. X'*X +k_exc*I 
k_chg = k_inc-k_exc;        % When we add a variable to the model, this is the amount that the matrix changes
p_i = parameters(6);        % Probability of including a varible
p_e = parameters(7);        % Probability of excluding a variable
[rows, cols] = size(XX);    
m = (cols-1)/2;             % m is the total number of predictor variables
    % Note: XX is the full matrix of regressors and contains both a sine and
    % cosine term for each predictor variable, as well as a constant term.
    % Only a subset will be used in each calculation

function PROB = tree(invG, detG, Am, C)
        C=C+1;                             % Depth in the tree
        if (C < m)
            % Exclusion = no update
            P1 = tree(invG, detG, Am,C);
            PROB=P1;
            
            % Inclusion = update our matrix
            Am(C)=1;                       % Include variable C
            E=invG(:,C)*k_chg; tr=E(C);    % EBIR step: E = Ck*Ek; tr = trace(Ck*Ek)
            invG=invG-E*invG(C,:)/(1+tr);  % Update matrix inverse
            detG=log(1+tr) +detG;          % Update determinant
            CC=C+m; % X contains both sine and cosine functions, so this will update the corresponding sine function
            E=invG(:,CC)*k_chg; tr=E(CC);  % EBIR step: E = Ck*Ek; tr = trace(Ck*Ek)
            invG=invG-E*invG(CC,:)/(1+tr); % Update matrix inverse
            detG=log(1+tr) +detG;          % Update determinant
            
            M1=(Am(6) || Am(9));           % Milankovitch Model
            M2=(Am(4) || Am(9));           % Harmonics of Obliquity
            M3=(Am(5) || Am(8));           % Harmonics of Precession
            
            if(M1+M2+M3+Am(7)-Am(9)<2)     % No conflicting models
                P2=tree(invG, detG, Am,C);
            else
                P2=[];
            end
            PROB=[PROB P2];
            
        else
            % The Final Step in the Recursion - Calculate the Probability
            % Exclusion = no update
            m_i=sum(Am);                    % m_i = Number of included variables
            hat_beta=invG*XTy;              % Regression coefficients
            I_Am= [Am Am 1];                % 1 is to include the constant term
            I_Am=I_Am*k_chg+k_exc;          % Diagonal matrix
            
            s_n = (y-X*hat_beta)'*(y-X*hat_beta)+I_Am.*hat_beta'*hat_beta+ v_0*sig_0;
            P1 = (m_i)*log(k_inc*p_i)+ (m-m_i)*log(k_exc*p_e) -v_n/2*log(s_n)-0.5*(detG);
            
            % Inclusion = update our matrix
            M3=(Am(5) || Am(8));             % Harmonics of Precession
           
            if(M3+Am(4)+Am(7)==0)            % If no model conflicts with M2, M3, or M4
                Am(C)=1;                     % C = 10, so we are adding a 'Milankovitch' (M1) term
                E=invG(:,C)*k_chg; tr=E(C);  % EBIR step: E = Ck*Ek; tr = trace(Ck*Ek)
                invG=invG-E*invG(C,:)/(1+tr);% Update matrix inverse
                detG=log(1+tr)+detG;          % Update matrix determinant
                CC=C+m;  % X contains both sine and cosine functions, so this will update the corresponding sine function
                E=invG(:,CC)*k_chg; tr=E(CC);% EBIR step: E = Ck*Ek; tr = trace(Ck*Ek)
                invG=invG-E*invG(CC,:)/(1+tr);% Update matrix inverse
                detG=log(1+tr)+detG;          % Update matrix determinant
                
                m_i=sum(Am);                 % Total number of included variables
                hat_beta=invG*XTy;           % Regression coefficients
                I_Am = [Am Am 1];
                I_Am = I_Am*k_chg+k_exc;     % Diagonal matrix
                
                s_n = (y-X*hat_beta)'*(y-X*hat_beta)+I_Am.*hat_beta'*hat_beta+ v_0*sig_0;
                P2 = (m_i)*log(k_inc*p_i)+ (m-m_i)*log(k_exc*p_e) -v_n/2*log(s_n)-0.5*(detG);
            
            else P2=[];
            
            end
            
            PROB=[P1 P2];
        
        end
        
    end

PROB=zeros(stop-start+1,N)-Inf;
for i=start:stop
    for j=i+2*cols:N	% 2*cols to ensure twice as many data points as parameters being estimated
        if(x(j)-x(i)>=min_chgpt)    % min_chgpt is the minimum distance between adjacent change points
            n = j-i+1;              % Length of data segment being considered
            X=XX(i:j,:); y=Y(i:j);  % Since we use only a substring of the data
            v_n = v_0+n;            % Parameter for posterior disdtribution on the error variance, sigma^2
            G=X'*X+k_exc*eye(2*m+1);% Create initial matrices - null model plus constant 
            G(2*m+1,2*m+1)=G(2*m+1,2*m+1)+k_chg; % Includes constant term
            invG=inv(G); XTy=X'*y; detG=log(det(G));
            Am=zeros(1,m);          % Vector indicating set of included predictor variables
            P=tree(invG, detG, Am,0);
            
            % Model Averaging
            PROB(i-start+1,j)=logsumlog(P)+gammaln(v_n/2)-n*log(pi)+v_0/2*log(v_0*sig_0) -gammaln(v_0/2); 
        end
    end
end





end
