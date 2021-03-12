function x = paramsToVec(obj,params,varargin)
% x = paramsToVec(obj,params,varargin)
%
% Convert vector form of parameters to struct
%
% Key/value pairs
%   'UseNoiseParam' - true/false (default false).  Use the noise parameter?
%
% 03/01/21  dhb  Wrote it from QCM version.

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
for i = 1:obj.nChannels/2
    x(i) = params.channelWeightsPos(i);
end
x(obj.nChannels/2+1) = params.crfAmp;
x(obj.nChannels/2+2) = params.crfExponent;
x(obj.nChannels/2+4) = params.crfSemi;
x(obj.nChannels/2+5) = params.expFalloff;
x(obj.nChannels/2+5) = params.crfOffset;

% Optional inclusion of noise
if (p.Results.UseNoiseParam)
    x(obj.nChannels/2+6) = params.noiseSd;
end

% Transpose to match convention
x = x';

end