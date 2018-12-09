function [params,paramsLb,paramsUb] = defaultParams(obj,varargin)
% [params,paramsLb,paramsUb] = defaultParams(obj,varargin)
%
% Set objects params to default values as well as provide reasonable lower
% and upper bournds.
%
% All three returns are in struct array form, use paramsToVec on the structs to
% get vector form.
%
% Optional key/value pairs
%    'defaultParams'     - (struct, default empty). If not empty, take the
%                          passed structure array as defining the default
%                          parameters.

% History:
%   12/09/18  dhb          Wrote it.

% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;
p.addParameter('defaultParams',[],@(x)(isempty(x) | isstruct(x)));
p.parse(varargin{:});

if (isempty(p.Results.defaultParams))
    %% Default parameters for Q
    %
    % This is number of directions specific
    params = tfeNRInitializeParams(nDirections);

    
    %% The Naka-Rushton
    params.crfAmp = 1;
    params.crfSemi = 1;
    params.crfExponent = 2;
    params.crfOffset = 0;
    
    % Exponential falloff (not used?)
    params.expFalloff = 0.3;
    
    % Noise level (used for simulations)
    params.noiseSd = 0.2;
    
else
    params = p.Results.defaultParams;
end

%% Set up search bounds
ampLowBound = 0; ampHighBound = 5;
semiLowBound = 0.01; semiHighBound = 10;
expLowBound = 0.01; expHighBound = 10;
if (obj.lockOffsetToZero)
    offsetLowBound = 0;
    offsetHighBound = 0;
else
    offsetLowBound = -ampHighBound;
    offsetHighBound = ampHighBound;
end

% These don't currently get searched over,
% so we just set them at their initial values.
expFalloffLowBound = params.expFalloff;
expFalloffHighBound = params.expFalloff;
noiseSdLowBound = params.noiseSd;
noiseSdHighBound = params.noiseSd;

% Pack bounds into vector form of parameters.
for ii = 1:nIndDirections
    NRParamsLow(ii).crfAmp = ampLowBound;
    NRParamsLow(ii).crfSemi = semiLowBound;
    NRParamsLow(ii).crfExponent = expLowBound;
    NRParamsLow(ii).crfOffset = offsetLowBound;
    NRParamsLow(ii).expFalloff = expFalloffLowBound;
    NRParamsLow(ii).noiseSd = noiseSdLowBound;
end
for ii = 1:nIndDirections
    NRParamsHigh(ii).crfAmp = ampHighBound;
    NRParamsHigh(ii).crfSemi = semiHighBound;
    NRParamsHigh(ii).crfExponent = expHighBound;
    NRParamsHigh(ii).crfOffset = offsetHighBound;
    NRParamsHigh(ii).expFalloff = expFalloffHighBound;
    NRParamsHigh(ii).noiseSd = noiseSdHighBound;
end

end