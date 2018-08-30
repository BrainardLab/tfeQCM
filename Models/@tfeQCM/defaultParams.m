function [params,paramsLb,paramsUb] = defaultParams(obj,varargin)
% [params,paramsLb,paramsUb] = defaultParams(obj,varargin)
%
% Set objects params to default values as well as provide reasonable lower
% and upper bournds.
%
% All three returns are in struct form, use paramsToVec on the structs to
% get vector form.

%% COMMENT DEFAULT PARAMS KEY/VALUE PAIR AS WELL NEW CODE.

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
            params.Qvec = [1 1 0 0 0];
            
        case 2
            params.Qvec = [1 0 ];
    end
    
    %% The Naka-Rushton
    params.crfAmp = 1;
    params.crfSemi = 1;
    params.crfExponent = 2;
    params.offset = 0;
    
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
        % has an implicit value of 1.
        paramsLb.Qvec = [1e-2 1e-2 -360 -360 -360];
        paramsUb.Qvec = [1e2 1e2 360 360 360];
        
        
    case 2
        paramsLb.Qvec = [1e-2 -360];
        paramsUb.Qvec = [1e2 360];
end

%% Lower bounds
paramsLb.crfAmp = 1e-1;
paramsLb.crfSemi = 1e-2;
paramsLb.crfExponent = 1e-2;
paramsLb.expFalloff = 1e-1;
paramsLb.noiseLevel = 0;
paramsLb.offset = -2;

%% Upper bounds
paramsUb.crfAmp = 1e1;
paramsUb.crfSemi = 1e2 ;
paramsUb.crfExponent = 10;
paramsUb.expFalloff = 1e1;
paramsUb.noiseLevel = 100;
paramsUb.offset = 2;
end