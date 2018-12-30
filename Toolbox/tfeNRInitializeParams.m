function NRParams = tfeNRInitializeParams(nDirections)
% Initialize the NR parameter struct array
%
% Syntax:
%    NRParams = tfeNRInitializeParams(nDirections)     
%
% Description
%    Set up struct array of Naka-Rushton parameters for the given number of
%    directions.
%
% Inputs:
%    nDirections     - Number of color directions.
%
% Outputs:
%    NRParams        - Struct array of parameters.  Fields are:
%                         crfAmp - Max amplitude (default 1)
%                         crfSemi - Semi saturation constant (default 0.06)
%                         crfExponent - Exponent in NR function (default 2)
%                         crfOffset - Additive offest (default 0)
%                         expFalloff - Exponential falloff constatn for tfe neural
%                                      response modeling (default 0.3)
%                         noiseSd 
%                    
%
% Optional key/value pairs:
%    None.
%
% See also: tfeNRParamsToVec, tfeNRVecToParams
%

% History:
%   12/10/18  dhb  Header comments.

% Reasonable defaults
Rmax = 1;
sigma = 0.05;
n = 2;
offset = 0;
expFalloff = 0.3;
noiseSd = 0.2;

% Set up the struct array
for ii = 1:nDirections
    NRParams(ii).crfAmp = Rmax;
    NRParams(ii).crfSemi = sigma;
    NRParams(ii).crfExponent = n;
    NRParams(ii).crfOffset = offset;
    NRParams(ii).expFalloff = expFalloff;
    NRParams(ii).noiseSd = noiseSd;
end

end