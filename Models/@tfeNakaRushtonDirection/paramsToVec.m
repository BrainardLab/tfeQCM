function x = paramsToVec(obj,params,varargin)
% Convert vector form of parameters to struct array
%
% Syntax:
%    x = paramsToVec(obj,params,varargin)
%
% Optoinal key/value pairs
%    None.

% Parse input.
p = inputParser; p.KeepUnmatched = true;
p.addRequired('params',@isstruct);
p.parse(params,varargin{:});
params = p.Results.params;

% Do it
paramsvec = tfeNRParamsToVec(NRParams);

end