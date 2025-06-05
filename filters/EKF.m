function [X_plus,Qxx_plus] = EKF(TS1,TS2,l,X,Qxx,Qww,Phi,trueStates,constraint,ep_first,ep_last)
%   function of the Extended Kalman Filter (EKF) algorithm
% 
%   INPUT:
%       TS1: position of the first total station (TS1)
%       TS2: positiion of the second total station (TS2)
%       l: matrix of the measured distance and angle values by means of the
%       total stations
%       X: state vector
%       Qxx: VCM of the states
%       Qww: VCM of the process noise 
%       Phi: transition matrix 
%       trueStates: matrix of the true trajectory values -> [x y constraint]
%       constraint: true -> constraint applies
%                   false -> constraint does not apply -> for the first
%                   exercise, set it to false
%       ep_first: the first epoch number
%       ep_last: the last epoch number
%
%   OUTPUT:
%       X_plus: filtered state vector 
%       Qxx_plus: VCM of the filtered states 

X_plus = zeros(size(l,1),size(X,1));
Qxx_plus = cell(size(l,1),1);
for ep = ep_first:ep_last
    %% measurement accuracies
    % store the measurement accuracies in a structure called sigma ->
    % sigma.d1, sigma.a1, sigma.d2, sigma.a2
    l_ep=l(ep,:);
    sigma.d1 = 40/1000 + 20*10^(-6)*l_ep(1); % [m]  40.0 mm + 20 ppm
    sigma.a1 = 20/1000*0.005*pi; % 1 gon = 0.005*pi [rad] 20.0 mgon
    sigma.d2 = 25/1000 + 15*10^(-6)*l_ep(3); % [m] 25.0 mm + 15 ppm
    sigma.a2 = 12/1000*0.005*pi; % [rad] 12.0 mgon
    
    %% form the VCM of the measurement accuracies and call it Qll %%
    v = [sigma.d1,sigma.a1,sigma.d2,sigma.a2];
    Qll = diag(v.^2);

    %% prediction step of the filter %%
    % fill out the "ffun" function to derive the predicted state vector
    % (predX) and VCM of the states (predQxx)

    if ep ~= 1                % update state and cov_matrix except epoch 1
        X = X_plus_ep;
        Qxx = Qxx_plus_ep;
    end
    [predX,predQxx] = ffun(X,Qxx,Phi,Qww);   % kinematic model, use it also at epoch 1

    %% update step of the filter %%
    % fill out the "design_matrix" function to derive the design matrix
    
    A = design_matrix(predX,TS1,TS2);

    % calculate the Kalman gain (K)
    K = predQxx * A' * inv(A * predQxx * A' + Qll);

    % fill out the "hfun" function to estimate the observations based on
    % the most recent state vector
    h = hfun(TS1,TS2,predX);
    
    % estimate the innovation vector (dl)
    dl = l(ep,:) - h;

    % update the predicted state vector and its corresponding VCM
    size_ide=size(K * A); % size of identity matrix
    X_plus_ep = predX + K * dl';
    Qxx_plus_ep = (eye(size_ide) - K * A) * predQxx;

    X_plus(ep,:) = X_plus_ep';
    Qxx_plus{ep,1} = Qxx_plus_ep;
    
    % apply the constraints if necessary -> NOT for the first exercise
    if (constraint == true) && trueStates(ep,3) == 1
        R = 9.00;   % radious of a sector of the trajectory
        D = D_Matrix(predX,TS1);  % xk_minus as input because of Taylor expansion 
        d = d_Vector(R,D,predX,TS1);
        W = eye(size(D,2));
        X_plus_ep = X_plus_ep - inv(W) * D' * inv(D * inv(W) * D' ) * (D * X_plus_ep - d);
        Qxx_plus_ep = Qxx_plus_ep - Qxx_plus_ep * D' * inv(D * Qxx_plus_ep * D' ) * D *Qxx_plus_ep;

        X_plus(ep,:) = X_plus_ep';
        Qxx_plus{ep,1} = Qxx_plus_ep;

    end
    
end
end