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

% Dimension check
if (obj.dimension ~= 2)
    error('LCM only implemented in 2 dimensions');
end

if (isempty(p.Results.defaultParams))
    % LCM parameters
    params.nChannels = 6;
    params.startCenter = 0;
    
    % Check
    if (rem(nChannels,2) ~= 0)
        error('nChannels must be even');
    end
    
    % Criterion response
    params.criterionResp = 2;
    
    % Angle support
    params.angleSupport = angleSupport = 0:1:360;
    
    % Create the channels
    %
    % These have tuning described as half wave rectified sinusoids
    % squared. Compute response to unit contrast in each color
    % direction by regarding these as linear channels tuned to hue
    % angle.
    %
    % Set channel center points
    centerSpacing = 360/nChannels;
    centerLocations = startCenter:centerSpacing:360-centerSpacing+params.startCenter;
    for ii = 1:params.nChannels
        params.underlyingChannels(ii,:) = cosd(params.angleSupport-centerLocations(ii));
        params.underlyingChannels(ii,sign(params.underlyingChannels(ii,:)) == -1) = 0;
        params.underlyingChannels(ii,:) = params.underlyingChannels(ii,:).^2;
    end
    
    % Initial weights
    params.channelWeightsPos = ones(nChannels/2,1)';
    
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
paramsLb = params;
paramsLb.channelWeightsPos = -1e6*ones(nChannels/2,1)';
paramsLb.crfAmp = 1e-1;
paramsLb.crfSemi = 1e-2;
paramsLb.crfExponent = 1e-2;
paramsLb.expFalloff = 1e-1;
paramsLb.noiseLevel = 0;
paramsLb.crfOffset = -2;

%% Upper bounds
paramsUb = params;
paramsUb.channelWeightsPos = 1e6*ones(nChannels/2,1)';
paramsUb.crfAmp = 3;
paramsUb.crfSemi = 1e2;
paramsUb.crfExponent = 10;
paramsUb.expFalloff = 1e1;
paramsUb.noiseLevel = 100;
paramsUb.crfOffset = 2;

end