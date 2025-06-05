function [X_plus,Qxx_plus] = EnKF(TS1,TS2,l,X,Qxx,Qww,Phi,N,ep_first,ep_last)
%   function of the Ensemble Kalman Filter (EnKF) algorithm
%
%   INPUT:
%       TS1: position of the first total station (TS1)
%       TS2: position of the second total station (TS2)
%       l: matrix of the measured distance and angle values by means of the
%       total stations
%       X: state vector
%       Qxx: VCM of the states
%       Qww: VCM of the process noise
%       Phi: transition matrix
%       N: number of the generated ensembles
%       ep_first: the first epoch number
%       ep_last: the last epoch number
%
%   OUTPUT:
%       X_plus: filtered state vector
%       Qxx_plus: VCM of the filtered states
X_plus = zeros(size(l,1),size(X,1));
Qxx_plus = cell(size(l,1),1);

%% generate the ensembles and store them in a variable called "S" %%

 S = mvnrnd(X,Qxx,N); % generate points (N*4)

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
    % fill out the "ffunStates" function to derive the predicted ensemble
    % points (continue with the name of the variable being "S")
    
    if ep ~= 1
        W = mvnrnd(zeros(size(X)),Qww,N); % generate noise (N*4)
        S = ffunStates(S',Phi)'+ W; 
    end

    % derive the mean value of each state from the ensembles and call it XPred 
    XPred = mean(S);  % (1*4)

    % derive the error matrix of the ensembles and call it "Ex" 
    Ex = S - XPred;
    %% update step of the filter %%
    % fill out the "hfun" function to estimate the observations based on
    % the most recent derived ensembles. call the predicted observations
    % "SigmaP_PredObv".
    for i = 1:size(S,1)
        SigmaP_PredObv(i,:) = hfun(TS1,TS2,S(i,:)');
    end
    hfun_S = SigmaP_PredObv;
    % add noise to the derived deterministic predcited observations by
    % taking into account the measurement accuracies (Qll) and using the
    % multivariate normal distribution. continue with the name of the
    % variable being "SigmaP_PredObv".
    Q_ll = mvnrnd(zeros(size(X)),Qll,N);
    SigmaP_PredObv = SigmaP_PredObv + Q_ll;

    % derive the mean of the predicted observations
    lPred = mean(SigmaP_PredObv);

    % derive the error matrix of the observations and call it "El" 
    El = SigmaP_PredObv - lPred;

    % derive the VCM of the ensembles and call it Qxx 
    Qxx = 1/(N-1) * Ex' * Ex ;

    % derive the VCM of the observations and call it Qll 
    Qll = 1/(N-1) * El' * El ;

    % derive the cross VCM between the ensembles and observations and call it "Qxl" 
    Qxl = 1/(N-1) * Ex' * El ;

    % derive the Kalman gain and call it "K" 
    K = Qxl * inv(Qll);  %  4 * 4

    % update the ensembles based on the derived Kalman gain 
    l_matrix = repmat(l(ep,:),N,1); % N*4 matrix for 1*4 observation l
    update_S = S + (K * (l_matrix  - hfun_S)')';

    %% derive the state vector "X" along with its VCM "Qxx" %%
    X_plus_ep = mean(update_S);
    update_Ex = update_S - X_plus_ep;
    Qxx_plus_ep = 1/(N-1) * update_Ex' * update_Ex ;

    X = X_plus_ep;      % previous state as input (Markov's order)
    Qxx = Qxx_plus_ep;
    S = mvnrnd(X,Qxx,N); % generate points (N*4)

    X_plus(ep,:) = X_plus_ep';
    Qxx_plus{ep,1} = Qxx_plus_ep;
    
end
end

