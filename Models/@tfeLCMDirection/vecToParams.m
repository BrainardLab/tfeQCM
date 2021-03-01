function params = vecToParams(obj,x,varargin)
% params = vecToParams(obj,x,varargin)
%
% Convert vector form of parameters to struct
%
% Key/value pairs
%   'UseNoiseParam' - true/false (default false).  Use the noise parameter?

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

params.Qvec(1:2) = x(1:2)';
params.crfAmp = x(3);
params.crfExponent = x(4);
params.crfSemi = x(5);
params.expFalloff = x(6);
params.crfOffset = x(7);

% Optional inclusion of noise
if (p.Results.UseNoiseParam)
    params.noiseSd = x(8);
end
        
end