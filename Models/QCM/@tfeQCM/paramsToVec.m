function x = paramsToVec(obj,params,varargin)
% x = paramsToVec(obj,params,varargin)
%
% Convert vector form of parameters to struct
%
% Key/value pairs
%   'UseNoiseParam' - true/false (default false).  Use the noise parameter?

% Parse input. At the moment this does type checking on the params input
% and has an optional key value pair that does nothing, but is here for us
% as a template.
p = inputParser;
p.addRequired('params',@isstruct);
p.addParameter('UseNoiseParam',false,@islogical);
p.parse(params,varargin{:});
params = p.Results.params;

switch obj.dimension
    case 3
        % Take the parameter structure into a vector
        for i = 1:5
            x(i) = params.Qvec(i);
        end
        x(6) = params.crfAmp;
        x(7) = params.crfExponent;
        x(8) = params.crfSemi;
        x(9) = params.expFalloff;
        x(10) = params.offset;
        % Optional inclusion of noise
        if (p.Results.UseNoiseParam)
            x(11) = params.noiseSd;
        end
    case 2
        for i = 1:3
            x(i) = params.Qvec(i);
        end
        x(4) = params.crfAmp;
        x(5) = params.crfExponent;
        x(6) = params.crfSemi;
        x(7) = params.expFalloff;
        x(8) = params.offset;
        % Optional inclusion of noise
        if (p.Results.UseNoiseParam)
            x(9) = params.noiseSd;
        end
end


end