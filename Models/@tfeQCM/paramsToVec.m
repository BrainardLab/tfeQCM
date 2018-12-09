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
        x(10) = params.crfOffset;
        
        % Optional inclusion of noise
        if (p.Results.UseNoiseParam)
            x(11) = params.noiseSd;
        end
    case 2
        for i = 1:2
            x(i) = params.Qvec(i);
        end
        x(3) = params.crfAmp;
        x(4) = params.crfExponent;
        x(5) = params.crfSemi;
        x(6) = params.expFalloff;
        x(7) = params.crfOffset;
        
        % Optional inclusion of noise
        if (p.Results.UseNoiseParam)
            x(8) = params.noiseSd;
        end
end


end