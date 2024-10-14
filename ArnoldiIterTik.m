function [Xapprox,RelResVec,RREvec,trackAlpha,trackExp,trackTikStep] = ArnoldiIterTik(A,bn,b,x,q_Orig,x0,NoiseLevel,maxIter)
% Description:
%   Modular function for Arnoldi iterative Tikhonov based on the 
%   Arnoldi-Tikhonov method published in 23' in Numerical Algorithms.
%   Method terminates according to the discrepancy principle.
%
%   Written by: Lucas Onisk, lucas.onisk@gmail.com
%   Last Edits: April 9th, 2023
%
% Inputs:
%           A - square matrix operator
%           bn - right-hand side vector of linear system with noise
%           b - right-hand side vector of linear system WITHOUT noise (used
%               for computing 'delta' in IAT algorithm based on 
%               Donatelli-Hanke scheme)
%           x - true solution of linear system (used in computing RRE - see
%               below)
%           q_Orig - a constant typically set at 0.6 - 0.8 which initially
%                    asserts how much one trusts the approximating model 
%                    and so is used to undersolve the reg. param. at each 
%                    iteration. As iterations proceed towards termination, 
%                    this number --> 1.
%           x0 - An initial guess for the IAT method. Numerical studies
%                suggest that the zero vector is often a good choice.
%           NoiseLevel - A priori knowldge about the noise/error
%                        contaminating the linear system (witten 
%                        in form 0.01 - for 1% noise)
%           maxIter - the maximum number of iterations to attempt
%
% Outputs:
%           Xapprox - final approximate solution vector of the method
%           RelResVec - a vector containing the iterative progression of
%                       the approximate relative residual errors 
%                       (used in method progression tracking)
%           RREvec - a vector containing the iterative progression of
%                       the relative reconstructive errors
%                       (used in method progression tracking)
%           trackAlpha - a vector containing the regularization parameters
%                        of the IAT method for each Tikhonov iterative step
%           trackExp - the number of Arnoldi expansions by the method
%                      (expect this number to be greater than or equal 
%                       to the number of Tikhonov type steps - see next 
%                       line down)
%           trackTikStep - the number of Tikhonov type steps the method
%                          takes (i.e. how many improvement vectors are 
%                          added to previous approximate solutions)
%
%% General alg. details
trackAlpha = [];
rho = 0.001; %used in AIT alg.
tau = (1+2*rho)/(1-2*rho);
delta = norm(bn-b,2);
breakout = tau*NoiseLevel;
Xapprox = x0;
r = bn - A*Xapprox;% make a starting residual vector
RelResVec = norm(r,2)/norm(bn,2);
normX = norm(x,2);
RREvec = norm(Xapprox - x,2)/normX; % RRE - reconstructive restoration error

%% Arnoldi Iterative Tikhonov

Q = bn/norm(bn,2); %start Q matrix
H = zeros(2,1); %start Hessenberg matrix

trackExp = 1; %each tik. iter. expand this vector equal to size of residual tracker; then add 1 for each expansion necessary for residual number
trackTikStep = 0; %tracking number of times Tikhonov iteration occurs

for i = 1:maxIter
    v = A*Q(:,i);
    for j = 1:i
        H(j,i) = Q(:,j)'*v;
        v = v - H(j,i)*Q(:,j);
    end
    
    H(i+1,i) = norm(v,2); %after this line we have H_{i+1,i}
    Q = [Q, v/H(i+1,i)]; %after this line we have Q_{i+1} where i = # columns
    
    tauK = norm(r(:,size(r,2)),2)/delta; %we want the last residual we have
    q = max(q_Orig,(2*rho + (1 + rho)/tauK));
        
    % Newton Method Section
    [W,S,~]=svd(H);
    singVals = diag(S*S'); %will now be (i+1)x1 instead of ix1
    rHat = W'*Q'*r(:,size(r,2)); %we want the last residual we have
    rHatsq = rHat.^2;
    C = q^2*sum(rHatsq);
    
    if rHatsq(length(rHatsq)) < q^2*sum(rHatsq)
        %if this passes then we should expect a unique reg. parameter
    else
        trackExp(length(trackExp)) = trackExp(length(trackExp)) + 1; % saying here that for said residual we needed to expand once more outward
        %expand H to continue Arnoldi
        H = [H ; zeros(1,i)]; %append single row of zeros under columns
        H = [H zeros(i+2,1)]; %append new column of zeros
        continue
    end

    beta = 0;
    for y = 1:8 %an arbitrary no. of newton steps
        f = rHat(length(rHat))^2 - C; %last entry of rHat needs to be accounted for
        fprime = 0;

        for k = 1:i %should only be "i" non-zero sing. vals for S*S'
            f = f + (rHat(k)^2)*((beta(length(beta))*singVals(k) + 1)^(-2));
            fprime = fprime + ((rHat(k)^2)*singVals(k))*((beta(length(beta))*singVals(k) + 1)^(-3));
        end

        fprime = (-2)*fprime;
        betaNew = beta(length(beta)) - f/fprime;
    end

    alpha = 1/betaNew;

    trackAlpha = [trackAlpha,alpha]; %track regularization parameter
    
    %%-------------------------------------%%
    
    % iterative tikhonov step

    h = Q(:,1:size(Q,2)-1)*((H'*H + alpha(length(alpha))*eye(i))\(H'*Q'*r(:,size(r,2))));
    trackTikStep = trackTikStep + 1; %tracking number of times Tikhonov iteration occurs
    Xapprox = Xapprox + h;
    
    newRes = bn - A*Xapprox;
    r = [r,newRes]; %build residual matrix;     %expand on last resiudal matrix
    RelResVec = [RelResVec,norm(r(:,size(r,2)),2)/norm(bn,2)];
    RREvec = [RREvec,norm(Xapprox - x,2)/normX];

    % Check if need to breakout according to discrepancy principle
    if norm(r(:,size(r,2)),2)/norm(bn,2) < breakout
        fprintf('\n \n ***Breakout condition was satisfied***')
        break
    end
    
    % If residual not satisfied, continue...
    H = [H ; zeros(1,i)]; %append single row of zeros under columns
    H = [H zeros(i+2,1)]; %append new column of zeros
    
    trackExp = [trackExp 1]; %only adds another if we're still going
   
end
end