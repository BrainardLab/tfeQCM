function params = vecToParams(obj,x,varargin)
% params = vecToParams(obj,x,varargin)
%
% Convert vector form of parameters to struct
%
% Optional key/value pairs
%   

% Parse input. At the moment this does type checking on the params input
% and has an optional key value pair that does nothing, but is here for us
% as a template.
p = inputParser; p.KeepUnmatched = true;
p.addRequired('x',@isnumeric);
p.parse(x,varargin{:});
x = p.Results.x;

% Do the conversion
params = tfeNRVecToParams(x,obj.nDirections);

end
