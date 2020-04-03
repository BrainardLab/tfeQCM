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
%   11/20/18  dhb, mab     Keep minor axis smaller than major in limits,
%                          for first two. 
%   11/24/18  dhb          Don't use degenerate default ellipse parameters.
%   12/09/18  dhb          For two-dimensional model, angle in range [-90,90]

% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;
p.addParameter('defaultParams',[],@(x)(isempty(x) | isstruct(x)));
p.parse(varargin{:});

if (isempty(p.Results.defaultParams))
    %% Default parameters for Q
    %
    % This is dimension specific
    switch obj.dimension
        case 3
            % Quadratic parameters
            %
            % We only store 5 because we handle the amplitude of the response
            % in the amplitude parameter below.  The first axis of the ellipse
            % has an implicit value of 1.
            params.Qvec = [0.5 0.25 0 0 0];
            
        case 2
            params.Qvec = [0.5 0];
    end
    
    %% The Naka-Rushton
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

switch obj.dimension
    case 3
        % Quadratic parameters
        %
        % We only store 5 because we handle the amplitude of the response
        % in the amplitude parameter below.  The first axis of the ellipse
        % has an implicit value of 1.  We want that to be the major axis,
        % so we bound the other axes at 1.
        paramsLb.Qvec = [1e-2 1e-2 -360 -360 -360];
        paramsUb.Qvec = [1 1 360 360 360];
           
    case 2
        paramsLb.Qvec = [1e-2 -180];
        paramsUb.Qvec = [1 180];
end

%% Lower bounds
paramsLb.crfAmp = 1e-1;
paramsLb.crfSemi = 1e-2;
paramsLb.crfExponent = 1e-2;
paramsLb.expFalloff = 1e-1;
paramsLb.noiseLevel = 0;
paramsLb.crfOffset = -2;

%% Upper bounds
paramsUb.crfAmp = 3;
paramsUb.crfSemi = 1e2 ;
paramsUb.crfExponent = 10;
paramsUb.expFalloff = 1e1;
paramsUb.noiseLevel = 100;
paramsUb.crfOffset = 2;

end