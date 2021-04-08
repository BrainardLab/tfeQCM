function scaleFactor = regressToFitIsoContrast(obj,params,isoContrast,varargin)
% Scale parameters so that LCM isocontrast contour best fits a given one.
%
% Syntax:
%    paramsOut = regessIsoContrast(obj,params,varargin)
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
p.addRequired('isoContrast',@isnumeric);
p.parse(params,isoContrast,varargin{:});

% Dimension check
if (obj.dimension ~= 2)
    error('LCM only implemented in 2 dimensions');
end

% Check
if (length(isoContrast) ~= length(obj.angleSupport))
    error('Wrong length for passed isocontrast');
end

% Get the LCM model isocontrast contour
LCMIsoContrast = obj.getIsoContrast(params);

% Find scalar
scaleFactor = LCMIsoContrast'\isoContrast';


end