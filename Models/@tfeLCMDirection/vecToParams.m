function params = vecToParams(obj,x,varargin)
% params = vecToParams(obj,x,varargin)
%
% Convert vector form of parameters to struct
%
% Key/value pairs
%   'UseNoiseParam' - true/false (default false).  Use the noise parameter?
%
% 03/01/21  dhb  Wrote it from QCM version.

% Parse input. At the moment this does type checking on the params input
% and has an optional key value pair that does nothing, but is here for us
% as a template.
p = inputParser;
p.addRequired('x',@isnumeric);
p.addParameter('UseNoiseParam',false,@islogical);
p.parse(x,varargin{:});
x = p.Results.x;

% Dimension check
if (obj.dimension ~= 2)
    error('LCM only implemented in 2 dimensions');
end

params.channelWeightsPos = x(1:obj.nChannels/2)';
params.crfAmp = x(obj.nChannels/2+1);
params.crfExponent = x(obj.nChannels/2+2);
params.crfSemi = x(obj.nChannels/2+3);
params.expFalloff = x(obj.nChannels/2+4);
params.crfOffset = x(obj.nChannels/2+5);

% Optional inclusion of noise
if (p.Results.UseNoiseParam)
    params.noiseSd = x(obj.nChannels/2+6);
end
        
end