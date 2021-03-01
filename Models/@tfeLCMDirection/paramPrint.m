function paramPrint(obj,params,varargin)
% paramPrint(obj,params,varargin)
%
% Print out useful things about parameters to the command window
%
% Key/value pairs
%   'PrintType' - string (default 'standard').  What to print.
%     'standard' - Standard print of parameters.
%
% 03/01/21  dhb  Wrote it from QCM version.

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
        fprintf('Channel weights (positive arm):\n');
        for ii = 1:obj.nChannels/2
            fprintf('Weight %d: %0.2g\n',params.channelWeightsPos(ii));
        end

        fprintf('CRF amplitude: %0.2f, CRF semi-saturation: %0.2f, CRF exponent: %0.2f\n',params.crfAmp,params.crfSemi,params.crfExponent);
        fprintf('Exponential filter time constant: %0.2f\n',params.expFalloff);
        fprintf('Offset constant: %0.2f\n',params.crfOffset);
        fprintf('\n');
    otherwise
        error('Unknown parameter print type passed')
        
end


