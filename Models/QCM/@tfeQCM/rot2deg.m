function deg = rot2deg(deg)
% deg2rot
%
% Description:
%   Takes in a 2x2 rotation matrix and returns an angle in degrees.
%
% Inputs:
%  rot             = A 2x2 rotation matrix
%
% Outputs:
%  deg             = The desired angle of rotation (in degrees) 
%
% Optional key/value pairs:
%   none
%
% Examples are provided in the source code.
%
% See also:
%

% History
%  8/23/18  mab  Created.

% Examples:
%{
    deg = rot2deg(deg);
%}

deg = acosd(rot(1,1));

end