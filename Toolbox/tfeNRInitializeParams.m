function NRParams = tfeNRInitializeParams(nIndDirections)
% Initialize the NR parameter struct array
%
% See also: tfeNRParamsToVec, tfeNRVecToParams
%

% Reasonable defaults
nParams = 4;
Rmax = 1;
sigma = 0.06;
n = 2;
offset = 0;
expFalloff = 0.3;
noiseSd = 0.2;

% Set up the struct array
for ii = 1:nIndDirections
    NRParams(ii).crfAmp = Rmax;
    NRParams(ii).crfSemi = sigma;
    NRParams(ii).crfExponent = n;
    NRParams(ii).crfOffset = offset;
    NRParams(ii).expFalloff = expFalloff;
    NRParams(ii).noiseSd = noiseSd;
end

end