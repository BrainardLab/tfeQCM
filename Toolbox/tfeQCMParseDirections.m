function [uniqueDirections,directionIndices] = tfeQCMParseDirections(directions)
% Parse a set of directions into unique directions
%
% Syntax:
%    [uniqueDirections,directionIndices] = tfeQCMParseDirections(directions)
%
% Desription:
%    Take stimulus directions and find unique directions as well as indices
%    into the columns of the passed directions matrix that correpond to
%    each unique direction.
%
%    Requires exact match as determined by unique. Could add some tolerance
%    if we discover that we need it.
%
% Inputs:
%    directions               - Array of directions.  Each should be a unit
%                               vector in a column of the array.
%
% Outputs:
%    uniqueDirections         - Array of unique directions, in columns.
%    directionIndices         - Cell array of indices, one vector for each
%                               entry of the cell array.  These index into the
%                               columns of the input that match each unique direction.
%
% Optional key/value pairs:
%    None.
%
% See also:
%

% History:
%   12/10/18  dhb  Pulled out into function of its own


%% Parse the stimuli into indDirections, contrasts and responses for each individual direction
% 
% The commented out call to uniquetol is almost good, but scrambles the
% order of the unique outputs.  One could fix that with a little work, if
% tolerance turns out to be an issue.
[indDirectionsTemp,~,whichColumnsOut] = unique(directions','rows','stable');  % uniquetol(directions',unique_tolerance,'ByRows',true);
indDirectionsTemp = indDirectionsTemp';
nIndDirections = size(indDirectionsTemp,2);
for ii = 1:nIndDirections
    whichColumns = find(whichColumnsOut == ii);
    directionIndices{ii} = whichColumns;
    uniqueDirections(:,ii) = indDirectionsTemp(:,ii);
end