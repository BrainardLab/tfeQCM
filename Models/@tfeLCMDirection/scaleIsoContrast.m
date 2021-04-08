function paramsOut = scaleIsoContrast(obj,params,scaleFactor,varargin)
% Scale isocontrast contour
%
% Syntax:
%    paramsOut = scaleIsoContrast(obj,params,scaleFactor,varargin)
%
% Description:
%    Scale the parameters so as to change the size of the isoresponse
%    contour by the passed scalar,without changing the model's predictions.
%
% Optional key/value pairs
%   None.

% History
%   03/12/21  dhb  Wrote it.

% Parse input.
p = inputParser;
p.addRequired('params',@isstruct);
p.addRequired('scalar',@isnumeric);
p.parse(params,scaleFactor,varargin{:});

% Dimension check
if (obj.dimension ~= 2)
    error('LCM only implemented in 2 dimensions');
end

% Check
if (scaleFactor < 1e-7 | scaleFactor > 1e7)
    error('That''s an awfully small or large scalar');
end

% Adjust weights and compensate with semi-satauration constant
paramsOut = params;
paramsOut.channelWeightsPos = params.channelWeightsPos/scaleFactor;
paramsOut.crfSemi = params.crfSemi/scaleFactor;

end