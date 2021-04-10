function [params,paramsLb,paramsUb] = defaultParams(obj,varargin)
% [params,paramsLb,paramsUb] = defaultParams(obj,varargin)
%
% Set objects params to default values as well as provide reasonable lower
% and upper bournds.
%
% All three returns are in struct form, use paramsToVec on the structs to
% get vector form.
%
% Optional key/value pairs
%    'defaultParams'     - (struct, default empty). If not empty, take the
%                          passed structure as defining the default
%                          parameters.
%
% History:
%   03/01/21  dhb          Keep LCM version from QCM version.

% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;
p.addParameter('defaultParams',[],@(x)(isempty(x) | isstruct(x)));
p.parse(varargin{:});

if (isempty(p.Results.defaultParams))
  
    % Initial weights
    params.channelWeightsPos = ones(obj.nChannels/2,1)';
    params.channelWeightsPos = params.channelWeightsPos/norm(params.channelWeightsPos);
    
    % The Naka-Rushton
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

%% Lower bounds
%
% By default we set lower and upper bounds of channel weights to one, to
% lock scale of isoresponse countour.
paramsLb.channelWeightsPos = [-1e6*ones(obj.nChannels/2,1)'];
%paramsLb.channelWeightsPos = [1e-6*ones(obj.nChannels/2,1)'];
paramsLb.crfAmp = 1e-1;
paramsLb.crfSemi = 1e-2;
paramsLb.crfExponent = 1e-2;
paramsLb.expFalloff = 1e-1;
paramsLb.noiseLevel = 0;
paramsLb.crfOffset = -2;

%% Upper bounds
paramsUb.channelWeightsPos = [1e6*ones(obj.nChannels/2,1)'];
paramsUb.crfAmp = 3;
paramsUb.crfSemi = 1e2;
paramsUb.crfExponent = 10;
paramsUb.expFalloff = 1e1;
paramsUb.noiseLevel = 100;
paramsUb.crfOffset = 2;

end