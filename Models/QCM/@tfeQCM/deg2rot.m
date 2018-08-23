function rot = deg2rot(deg)
% deg2rot
%
% Description:
%   Takes in an angle (degrees) and returns a 2x2 rotation matrix
%
% Inputs:
%  deg              = The desired angle of rotation (in degrees)
%
% Outputs:
%   rot             = A 2x2 rotation matrix 
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
    rot = deg2rot(45);
%}

rot = [cosd(deg),-1*sind(deg);sind(deg),cosd(deg)];

end