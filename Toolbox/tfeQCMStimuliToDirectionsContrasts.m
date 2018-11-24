function [directions,contrasts] = tfeQCMStimuliToDirectionsContrasts(stimuli)
% Compute directions and contrasts from stimuli
%
% Synopsis
%    [directions,contrasts] = tfeQCMStimuliToDirectionsContrasts(directions,stimuli)
%
% Description:
%    Convert stimuli to directions and contrasts.
%
%    There is some ambiguity for stimuli of length 0.  We return
%    a unit vector of passed dimension with 1 as the first entry,
%    and contrast of 0 in this case
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
% See also: tfeQCMDirectionsContrastsToStimuli
%

% History:
%   11/24/19  dhb    Wrote it for modularity.

%% Convert stimuli to directions/contrasts
nStimuli = size(stimuli,2);
directions = zeros(size(stimuli));
contrasts = zeros(1,nStimuli);
for ii = 1:nStimuli
    contrasts(ii) = norm(stimuli(:,ii));
    if (contrasts(ii) == 0)
        directions(:,ii) = zeros(size(stimuli(:,ii)));
        directions(1,ii) = 1;
    else
        directions(:,ii) = stimuli(:,ii)/contrasts(ii);
    end
end

