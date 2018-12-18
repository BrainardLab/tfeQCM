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

% N-R parameters
switch (p.Results.PrintType)
    case 'standard'
        for ii = 1:obj.nDirections
            fprintf('Naka-Rushton parameters for direction %d\n',ii);
            fprintf('\tCRF amplitude: %0.2f, CRF semi-saturation: %0.2f, CRF exponent: %0.2f, CRF offset: %0.2f\n',params(ii).crfAmp,params(ii).crfSemi,params(ii).crfExponent,params(ii).crfOffset);

        end
        %fprintf('Exponential filter time constant: %0.2f\n',params.expFalloff);
        fprintf('\n');
    otherwise
        error('Unknown parameter print type passed')    
end


