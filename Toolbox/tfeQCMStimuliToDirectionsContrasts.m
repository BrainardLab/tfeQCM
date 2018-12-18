function [directions,contrasts] = tfeQCMStimuliToDirectionsContrasts(stimuli,varargin)
% Compute directions and contrasts from stimuli
%
% Synopsis
%    [directions,contrasts] = tfeQCMStimuliToDirectionsContrasts(directions,stimuli)
%
% Description:
%    Convert stimuli to directions and contrasts.
%
%    There is some ambiguity for stimuli of length 0.  See key/value pairs
%    below for how this is handled.
%
%    Directions are rounded to precision places to avoid problems when
%    determining direction matches elsewhere.
%
% Inputs:
%    stimuli       - Stimuli in columns of the passed matrix.
% 
% Outputs:
%    directions    - Directions in columns and each
%                    column vector having a vector lenght of 1.
%    contrasts     - Row vector of contrasts. Multiply directions by
%                    corresponding contrasts to get stimuli.
%
% Optional key/value pairs
%    'zeroContrastDirection' - Vector with same row dimension as stimuli.
%                              This is returned as the direction for zero
%                              contrast entries.  Default is empty. If
%                              empty, the zero contrast direction
%                              corresponds to the first direciton found in
%                              stimuli with non-zero contrast.
%    'precision'             - Number of places to round directions to.
%                              Default 4.
%
% See also: tfeQCMDirectionsContrastsToStimuli
%

% History:
%   11/24/19  dhb    Wrote it for modularity.

%% Figure out dimension from stimuli
dimension = size(stimuli,1);

%% Parse key/value pairs
p = inputParser; 
p.addRequired('stimuli',@isnumeric);
p.addParameter('zeroContrastDirection',[],@isnumeric);
p.addParameter('precision',4,@isnumeric);
p.parse(stimuli,varargin{:});

%% Check
if (~isempty(p.Results.zeroContrastDirection) && size(p.Results.zeroContrastDirection,1) ~= dimension)
    error('Hey buddy. You passed a default zero contrast direction vector with the wrong dimension');
end

%% Convert stimuli to directions/contrasts
nStimuli = size(stimuli,2);
directions = zeros(size(stimuli));
contrasts = zeros(1,nStimuli);
gotNonZero = false;
for ii = 1:nStimuli
    contrasts(ii) = norm(stimuli(:,ii));
    if (contrasts(ii) > 0)
        gotNonZero = true;
        if (isempty(p.Results.zeroContrastDirection))
            zeroContrastDirection = stimuli(:,ii)/contrasts(ii);
        end
    end
end
if (~gotNonZero)
    error('At least one stimulus needs non-zero contrast');
end
if (~isempty(p.Results.zeroContrastDirection))
    zeroContrastDirection = p.Results.zeroContrastDirection;
end

%% Fill in directions
for ii = 1:nStimuli
    if (contrasts(ii) == 0)
        directions(:,ii) = zeroContrastDirection;
    else
        directions(:,ii) = stimuli(:,ii)/contrasts(ii);
    end
end

%% Round directions to some reasonable precision
directions = round(directions,p.Results.precision);

