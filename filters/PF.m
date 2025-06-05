function [X_plus,Qxx_plus] = PF(TS1,TS2,l,X,Qxx,Qww,Phi,NP,ep_first,ep_last)
%   function of the Particle Filter (PF) algorithm
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
%       NP: number of the particles
%       ep_first: the first epoch number
%       ep_last: the last epoch number
%
%   OUTPUT:
%       X_plus: filtered state vector
%       Qxx_plus: VCM of the filtered states
X_plus = zeros(size(l,1),size(X,1));
Qxx_plus = cell(size(l,1),1);

%% generate the particles and store them in a variable called "S" %%
S = mvnrnd(X,Qxx,NP); % generate points (N*4)

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
    % fill out the "ffunStates" function to derive the predicted sample
    % points (continue with the name of the variable being "S")
    if ep ~= 1
        W = mvnrnd(zeros(size(X)),Qww,NP); % generate noise (N*4)
        S = ffunStates(S',Phi)'+ W;
    end
    
    %% update step of the filter %%
    % fill out the "hfun" function to estimate the observations based on
    % the most recent derived particles. then, derive the difference
    % between these estimated observations and the real sensor data
    % (measurements). store the differences in a variable called "res".
    for i = 1:size(S,1)
        h_PredObv(i,:) = hfun(TS1,TS2,S(i,:)');
    end
    res = l(ep,:) - h_PredObv;
    % estimate the likelihood of the particles based on the derived
    % differences (res) and by using "mvnpdf" MATLAB library and store
    % these likelihoods in a vector called "q"

    q = mvnpdf(res,zeros(size(X))',Qll) + 10^(-99); 
   
    % derive the importance weights by using the derived likelihoods and
    % store the normalized importance weights in a vector called "w_tilde"
    w = (1 / NP) * q;
    w_tilde = w / sum(w);

    % resample the particles by using the "ResidualResampling" function 
    resample_Index = ResidualResample(1:NP, w_tilde');

    %% derive the state vector "X" along with its VCM "Qxx" %%
    resample_S = S(resample_Index',:); % X_new  
    X_plus_ep = mean(resample_S);
    Qxx_plus_ep = cov(resample_S);

    S = resample_S;      % previous state as input (Markov's order)

    X_plus(ep,:) = X_plus_ep';
    Qxx_plus{ep,1} = Qxx_plus_ep;
end
end

