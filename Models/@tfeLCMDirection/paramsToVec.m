function x = paramsToVec(obj,params,varargin)
% x = paramsToVec(obj,params,varargin)
%
% Convert vector form of parameters to struct
%
% Key/value pairs
%   'UseNoiseParam' - true/false (default false).  Use the noise parameter?

% Parse input.
p = inputParser;
p.addRequired('params',@isstruct);
p.addParameter('UseNoiseParam',false,@islogical);
p.parse(params,varargin{:});
params = p.Results.params;

% Dimension check
if (obj.dimension ~= 2)
    error('LCM only implemented in 2 dimensions');
end

% Expand-eroo!
for i = 1:params.nChannels/2
    x(i) = params.channelWeightsPos(i);
end
x(params.nChannels/2+1) = params.crfAmp;
x(params.nChannels/2+2) = params.crfExponent;
x(params.nChannels/2+3) = params.crfSemi;
x(params.nChannels/2+4) = params.expFalloff;
x(params.nChannels/2+5) = params.crfOffset;

% Optional inclusion of noise
if (p.Results.UseNoiseParam)
    x(params.nChannels/2+6) = params.noiseSd;
end

% Transpose to match convention
x = x';

end