function outIndex = ResidualResample(inIndex, weights);

% RESIDUALRESAMPLE  Residual Resampling Technik f�r den SIR-Filter.
% F�r Details siehe bitte Doucet. et.al. Da werden verschieden
% Resamplingstechnicken besprochen. Hier wurde die Technick auf S. 467 (Optimal State Estimation)
% Diese Technik ist bekannt als Residual Resampling. Es gibt noch vielzahl von
% varianten, die noch nicht programmiert worden
%   outIndex = ResidualResample(inIndex, weights)
%
%   INPUT
%          inIndex        (r-vector) Input Partikel-Indexies.
%          weights        (r-vector) Normalisieret Importance Gewichte
%   OUTPUT
%          outIndex          Resamplte Indicies der Partikel.
% Diese Routine wurde von Doucet. Et. all �bernommen
% Hamza Alkhatib 18.03.2008
%=============================================================================================

if (nargin ~= 2),
    error(' [ ResidualResample ] Nicht genuegende Input-Variablen.');
end

% Anzahl der Partikel
S = length(weights);       

%Intialisere Output-Index-Vec
outIndex = zeros(1,S);  

%=== RESIDUAL RESAMPLING  ==========================================================

N_kind= zeros(1,S);

% Zuerst Ineger Teil
weights_res = S*weights;
N_kind = fix(weights_res);

% Residualanzahl der Partikel zum Generieren
N_res = S-sum(N_kind);

if N_res

    weights_res = (weights_res-N_kind)/N_res;
    cumDist = cumsum(weights_res);   % Cumulative sum

    % Generate N_res sortierte  Zufallsvariablen aus der Uniform-Verteilung [0,1]
    u = fliplr(cumprod(rand(1,N_res).^(1./(N_res:-1:1))));
    j=1;
    for i=1:N_res
        while (u(1,i)>cumDist(1,j))
            j=j+1;
        end
        N_kind(1,j)=N_kind(1,j)+1;
    end;

end;


index=1;
for i=1:S
    if (N_kind(1,i)>0)
        for j=index:index+N_kind(1,i)-1
            outIndex(j) = inIndex(i);
        end;
    end;
    index = index+N_kind(1,i);
end










