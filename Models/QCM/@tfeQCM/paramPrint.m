function paramPrint(obj,params,varargin)
% paramPrint(obj,params,varargin)
%
% Print out useful things about parameters to the command window
%
% Key/value pairs
%   'PrintType' - string (default 'standard').  What to print.
%     'standard' - Standard print of parameters.

% Parse input. At the moment this does type checking on the params input
% and has an optional key value pair that does nothing, but is here for us
% as a template.
p = inputParser;
p.addRequired('params',@isstruct);
p.addParameter('PrintType','standard',@ischar);
p.parse(params,varargin{:});
params = p.Results.params;

% Quadratic parameters
switch (p.Results.PrintType)
    case 'standard'
        fprintf('Quadratic ellipse lengths: 1.00, %0.2f, %0.2f\n',params.Qvec(1),params.Qvec(2));
        fprintf('Quadratic ellipse angles (degs): %0.1f, %0.1f %0.1f\n',(180/pi)*params.Qvec(3),(180/pi)*params.Qvec(4),(180/pi)*params.Qvec(5));
        fprintf('CRF amplitude: %0.2f, CRF semi-saturation: %0.2f, CRF exponent: %0.2f\n',params.crfAmp,params.crfSemi,params.crfExponent);
        fprintf('Exponential filter time constant: %0.2f\n',params.expFalloff);
        fprintf('Offset constant: %0.2f\n',params.offset);
        fprintf('\n');
    otherwise
        error('Unknown parameter print type passed')
        
end