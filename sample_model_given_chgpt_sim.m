function V = sample_model_given_chgpt_sim(start, stop,XX,Y,parameters)

% Step 3.3: Sample a Sub-Model for the Interval Between adjacent Change Points
%
% The EBIR [Exact Bayesian Inference in Regression] algorithm (Ruggieri and 
% Lawrence 2012) calculates the exact probability of all possible
% subsets of predictor variables.  The recursion can be thought of as a tree  
% with the null model at the root.
%
% Output: V = A subset of the predictor variables sampled from the
% posterior distribution on the set of possible sub-models.

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

n = stop-start+1;               % Length of data segment being considered
v_n = v_0+n;                    % Parameter for posterior disdtribution on the error variance, sigma^2
X= XX(start:stop,:); y=Y(start:stop);   % Since we use only a substring of the data
G=X'*X+k_exc*eye(m);            % Create initial matrices - null model
invG=inv(G); XTy=X'*y; detG=log(det(G));
Am=zeros(1,m);                  % Vector indicating set of included predictor variables
P=tree(invG, detG, Am, 0);
total_P=logsumlog(P);           % Normalization constant
PROB=zeros(1,2^m);
PROB(:)=exp(P(:)-total_P);      % Normalize the vector - create the pdf P(Am | Y)

S=pick_k1(PROB);                % Draw a random sample from P(Am|Y)
%V=dec2binvec(S-1,m);           % Using Statistics Toolbox
V = zeros(1,m);
S=S-1;                          %Because the null model (all zeros) is in position 1 of vector
for i = m:-1:1
    V(i) = mod(S,2);
    S = floor(S/2);
end

end % of function
 