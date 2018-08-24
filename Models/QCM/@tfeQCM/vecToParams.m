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

switch obj.dimension
    case 3
        params.Qvec(1:5) = x(1:5)';
        params.crfAmp = x(6);
        params.crfExponent = x(7);
        params.crfSemi = x(8);
        params.expFalloff = x(9);
        params.offset = x(10);
        % Optional inclusion of noise
        if (p.Results.UseNoiseParam)
            params.noiseSd = x(11);
        end
    case 2
        params.Qvec(1:3) = x(1:3)';
        params.crfAmp = x(4);
        params.crfExponent = x(5);
        params.crfSemi = x(6);
        params.expFalloff = x(7);
        params.offset = x(8);
        % Optional inclusion of noise
        if (p.Results.UseNoiseParam)
            params.noiseSd = x(9);
        end
end
        
end