function PROB = DP_exact(start, stop,x,XX,Y,N, parameters)

% Step 1: Calculate the Probability Density of the Data
%
% The EBIR [Exact Bayesian Inference in Regression] algorithm (Ruggieri and 
% Lawrence 2012) calculates the exact probability of all possible
% subsets of predictor variables.  The recursion can be thought of as a tree  
% with the null model at the root.
%
% Output: f(Y_i:j)

min_chgpt = parameters(1); % Minimum distance between adjacent change points
v_0 = parameters(2);    % Prior parameter for the error variance, sigma^2
sig_0 = parameters(3);  % Prior parameter for the error variance, sigma^2
k_inc = parameters(4);  % Prior parameter for the regression coefficients, beta (k_i in manuscript)
k_exc = parameters(5);  % Prior parameter for the regression coefficients, beta (k_e in manuscript)
% The tree starts with the null model at its root, i.e. X'*X +k_exc*I 
k_chg = k_inc-k_exc;    % When we add a variable to the model, this is the amount that the matrix changes
p_i = parameters(6);    % Probability of including a varible
p_e = parameters(7);    % Probability of excluding a variable
[rows, m] = size(XX);   % m is the total number of predictor variables
    % Note: XX is the full matrix of regressors, only a subset will be used
    % in each calculation

    function PROB = tree(invG, detG, Am, C)
        C=C+1;                          % Depth in the tree
        if (C < m)
            % Exclusion = no update
            P1 = tree(invG, detG, Am,C);
            PROB=P1;
            
            % Inclusion = update our matrix
            Am(C)=1;                     % Include variable C
            E=invG(:,C)*k_chg; tr=E(C);  % EBIR step: E = Ck*Ek; tr = trace(Ck*Ek)
            invG=invG-E*invG(C,:)/(1+tr);% Update matrix inverse
            detG=log(1+tr) +detG;        % Update determinant
            
            P2=tree(invG, detG, Am,C);
            
            PROB=[PROB P2];
            
        else
            % The Final Step in the Recursion - Calculate the Probability
            % Exclusion = no update
            m_i=sum(Am);                % m_i = Number of included variables
            hat_beta=invG*XTy;          % Regression coefficients
            I_Am= Am;
            I_Am=I_Am*k_chg+k_exc;      % Diagonal matrix
            
            s_n = (y-X*hat_beta)'*(y-X*hat_beta)+I_Am.*hat_beta'*hat_beta+ v_0*sig_0;
            P1 = (m_i/2)*log(k_inc*p_i^2)+ (m-m_i)/2*log(k_exc*p_e^2) -v_n/2*log(s_n)-0.5*(detG);
            
            % Inclusion = update our matrix
            Am(C)=1;                     % Include the final variable
            E=invG(:,C)*k_chg; tr=E(C);  % EBIR step: E = Ck*Ek; tr = trace(Ck*Ek)
            invG=invG-E*invG(C,:)/(1+tr);% Update matrix inverse
            detG=log(1+tr)+detG;         % Update determinant
            
            m_i=sum(Am);                 % Total number of included variables
            hat_beta=invG*XTy;           % Regression coefficients
            I_Am = Am;
            I_Am = I_Am*k_chg+k_exc;     % Diagonal matrix
            
            s_n = (y-X*hat_beta)'*(y-X*hat_beta)+I_Am.*hat_beta'*hat_beta+ v_0*sig_0;
            P2 = (m_i/2)*log(k_inc*p_i^2)+ (m-m_i)/2*log(k_exc*p_e^2) -v_n/2*log(s_n)-0.5*(detG);
            
            PROB=[P1 P2];
            
        end
        
    end

PROB=zeros(stop-start+1,N)-Inf;
for i=start:stop
    for j=i+min_chgpt:N         % min_chgpt is the minimum distance between adjacent change points
        if(x(j)-x(i)>=min_chgpt)
            n = j-i+1;          % Length of data segment being considered
            X=XX(i:j,:); y=Y(i:j);  % Since we use only a substring of the data
            v_n=v_0+n;          % Parameter for posterior disdtribution on the error variance, sigma^2
            G=X'*X+k_exc*eye(m);% Create initial matrices - null model
            invG=inv(G); XTy=X'*y; detG=log(det(G));
            Am=zeros(1,m);      % Vector indicating set of included predictor variables 
            P=tree(invG, detG, Am, 0);
            
            % Model Averaging
            PROB(i-start+1,j)=logsumlog(P)+gammaln(v_n/2)-n*log(pi)+v_0/2*log(v_0*sig_0) -gammaln(v_0/2);
        end
    end
end




end
