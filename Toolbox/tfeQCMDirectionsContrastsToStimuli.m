function [stimuli] = tfeQCMDirectionsContrastsToStimuli(directions,contrasts)
% Compute stimuli from directions and contrasts
%
% Synopsis
%    [stimuli] = tfeQCMDirectionsContrastsToStimuli(directions,contrasts)
%
% Description:
%    Convert directions and contrasts to stimuli.
%
% Inputs:
%    directions    - Stimuli, with stimulus contrasts in columns and each
%                    column vector having a vector length of 1.
%    contrasts     - Row vector of contrasts. Multiply directions by
%                    contrasts to get stimuli.  Can have a scalar here,
%                    in which case it is in common across direcitons.
% 
% Outputs:
%    stimuli       - Stimuli corresponding to direcitons and contrasts
%
% See also: tfeQCMStimuliToDirectionsContrasts
%

% History:
%   11/24/19  dhb    Wrote it for modularity.

% Examples:
%{
    stimuli = tfeQCMDirectionsContrastsToStimuli([[1 0]', [0 1]'],[0.5 2])
    [directions, contrasts] = tfeQCMStimuliToDirectionsContrasts(stimuli)
%}

%% Convert directions/contrasts to stimuli
if (length(contrasts(:)) == 1)
    stimuli = directions*contrasts(1);
else
    stimuli = directions*diag(contrasts);
end

