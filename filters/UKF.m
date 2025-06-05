function [X_plus,Qxx_plus] = UKF(TS1,TS2,l,X,Qxx,Qww,Phi,gamma,ep_first,ep_last)
%   function of the Unscented Kalman Filter (UKF) algorithm
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
%       gamma: scale factor in UKF
%       ep_first: the first epoch number
%       ep_last: the last epoch number
%
%   OUTPUT:
%       X_plus: filtered state vector
%       Qxx_plus: VCM of the filtered states
X_plus = zeros(size(l,1),size(X,1));
Qxx_plus = cell(size(l,1),1);

%% calculte the number of the necessary Sigma points and store them in a variable called NoSig %%
NoSig = 2 * size(X,1) + 1;

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
    
    %% Generating the Sigma points of the previous time step %%
    % fill out the "UT_transform" function to derive the Sigma points
    % (SigmaPoints) and their corresponding weights (WSigma)
    n = size(X,1);
    [SigmaPoints,WSigma] = UT_transform(NoSig,X,Qxx,gamma,n);

    % plot(SigmaPoints(1,:),SigmaPoints(2,:),'*')
    %% prediction step of the filter %%
    % fill out the "ffunStates" function to derive the predicted Sigma points 
    % (predSigmaP)
    predSigmaP = ffunStates(SigmaPoints,Phi); %[4*4] * [4*9]

    % fill out the "WMean" function to derive the predicted state vector
    % (predX) that is a weighted mean of the predicted Sigma points
    % (predSigmaP)
    predX = WMean(predSigmaP,WSigma); % [4*9] * [9*1] X_minus
    
    % fill out the "QMean" function to derive the predicted VCM of the
    % states (predQxx)
    delta = predSigmaP - predX; %[4*9]
    predQxx = QMean(delta,WSigma) + Qww; % Qxx_minus

    % derive the new Sigma points based on the predicted states by using
    % the "UT_transform" function and call them newPredSigmaP with new
    % corresponding weights (newWSigmaP)
    [newPredSigmaP,newWSigmaP] = UT_transform(NoSig,predX,predQxx,gamma,n);
    %% update step of the filter %%
    % fill out the "hfun" function to estimate the observations (updSigmaP)
    % based on the most recent Sigma point
    for i = 1:size(newPredSigmaP,2)
        updSigmaP(:,i) = hfun(TS1,TS2,newPredSigmaP(:,i))';
    end
    % estimate the predicted observation vector based on the estimated
    % observations (updSigmaP) and by using the "WMean" function
    predL = WMean(updSigmaP,newWSigmaP);
    
    % estimate the predicted VCM of the observations (predQll) by using the
    % "QMean" function
    deltaL = updSigmaP - predL;
    predQll = QMean(deltaL,newWSigmaP) + Qll;
    
    % fill out the "crossQMean" to derive the cross VCM between the state
    % vector and observation vector and call it Qxl
    Qxl = crossQMean(delta,deltaL,newWSigmaP);

    % derive the Kalman gain (K)
    K = Qxl * inv(predQll);
    %% derive the updated state vector along with its corresponding VCM %% 
    X_plus_ep = predX + K * (l(ep,:)' - predL);
    Qxx_plus_ep = predQxx - K * predQll * K';
    
    X = X_plus_ep; % update state
    Qxx = Qxx_plus_ep;
    
    X_plus(ep,:) = X_plus_ep';
    Qxx_plus{ep,1} = Qxx_plus_ep;

end
end