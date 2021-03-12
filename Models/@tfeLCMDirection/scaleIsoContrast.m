function paramsOut = scaleIsoContrast(obj,params,scalar,varargin)
% paramsOut = scaleIsoContrast(obj,params,varargin)
%
% Description:
%    Scale the parameters so as to change the size of the isoresponse
%    contour by the passed scalar,without changing the model's predictions.
%
% Key/value pairs
%
% 03/12/21  dhb  Wrote it.

% Parse input.
p = inputParser;
p.addRequired('params',@isstruct);
p.addRequired('scalar',@isnumeric);
p.parse(params,scalar,varargin{:});

% Dimension check
if (obj.dimension ~= 2)
    error('LCM only implemented in 2 dimensions');
end

% Check
if (scalar < 1e-6)
    error('That''s an awfully small scalar');
end

% Adjust weights and compensate with semi-satauration constant
paramsOut = params;
paramsOut.channelWeightsPos = params.channelWeightsPos/scalar;
paramsOut.crfSemi = params.crfSemi/scalar;

end