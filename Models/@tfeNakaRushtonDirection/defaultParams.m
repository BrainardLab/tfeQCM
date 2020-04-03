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

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;
p.addParameter('defaultParams',[],@(x)(isempty(x) | isstruct(x)));
p.parse(varargin{:});

if (isempty(p.Results.defaultParams))
    % Set up structure
    params = tfeNRInitializeParams(obj.nDirections);
    
else
    params = p.Results.defaultParams;
end

%% Check on parameters
for ii = 1:obj.nDirections
    if (params(ii).noiseSd ~= params(1).noiseSd)
        error('Noise sd parameter not matched across directions in parameters struct array');
    end
    if (params(ii).expFalloff ~= params(1).expFalloff)
        error('Exp falloff parameter not matched across directions in parameters struct array');
    end
end

%% Set up search bounds
ampLowBound = 0; ampHighBound = 5;
semiLowBound = 0.01; semiHighBound = 10;
expLowBound = 0.01; expHighBound = 3;
if (obj.lockOffsetToZero)
    offsetLowBound = 0;
    offsetHighBound = 0;
else
    offsetLowBound = -ampHighBound;
    offsetHighBound = ampHighBound;
end

% These don't currently get searched over,
% so we just set them at their initial values.
expFalloffLowBound = params(1).expFalloff;
expFalloffHighBound = params(1).expFalloff;
noiseSdLowBound = params(1).noiseSd;
noiseSdHighBound = params(1).noiseSd;

% Pack bounds into vector form of parameters.
for ii = 1:obj.nDirections
    paramsLb(ii).crfAmp = ampLowBound;
    paramsLb(ii).crfSemi = semiLowBound;
    paramsLb(ii).crfExponent = expLowBound;
    paramsLb(ii).crfOffset = offsetLowBound;
    paramsLb(ii).expFalloff = expFalloffLowBound;
    paramsLb(ii).noiseSd = noiseSdLowBound;
end
for ii = 1:obj.nDirections
    paramsUb(ii).crfAmp = ampHighBound;
    paramsUb(ii).crfSemi = semiHighBound;
    paramsUb(ii).crfExponent = expHighBound;
    paramsUb(ii).crfOffset = offsetHighBound;
    paramsUb(ii).expFalloff = expFalloffHighBound;
    paramsUb(ii).noiseSd = noiseSdHighBound;
end

end