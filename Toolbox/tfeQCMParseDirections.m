function [uniqueDirections,directionIndices] = tfeQCMParseDirections(directions,varargin)
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
%    'precision'             - Number of places to round directions to.
%                              Default 4.
%
% See also:
%

% History:
%   12/10/18  dhb  Pulled out into function of its own

%% Parse key/value pairs
p = inputParser; 
p.addRequired('directions',@isnumeric);
p.addParameter('precision',4,@isnumeric);
p.parse(directions,varargin{:});

%% Round directions to some reasonable precision
directions = round(directions,p.Results.precision);

%% Find zero length stimuli and replace with first non-zero stimulus
%
% This has the effect of making sure we don't return a zero vector as
% a unique direction, and also that we assign an index to each zero vector
% that points at some direction. This is OK because eventually the contrast
% of the corresponding stimulus will be at 0 and the correct stimulus will 
% be synthesized no matter what the direction is.
vectorNorms = vecnorm(directions);
zeroIndex = find(vectorNorms == 0);
nonZeroIndex = find(vectorNorms >0);

if (isempty(nonZeroIndex))
    error('No non-zero directions passed');
end
if (~isempty(zeroIndex))
    for ii = 1:length(zeroIndex)
        directions(:,zeroIndex(ii)) = directions(:,nonZeroIndex(1));
    end
end

%% Parse the stimuli into indDirections
[indDirectionsTemp,~,whichColumnsOut] = unique(directions','rows','stable');
indDirectionsTemp = indDirectionsTemp';
nIndDirections = size(indDirectionsTemp,2);

%% Set up return.
% Zero length directions in input are arbitrarily assigned to
% the first direction.
for ii = 1:nIndDirections
    whichColumns = find(whichColumnsOut == ii);
    directionIndices{ii} = whichColumns;
    uniqueDirections(:,ii) = indDirectionsTemp(:,ii);
end