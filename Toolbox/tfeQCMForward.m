function [responses] = tfeQCMForward(params,stimuli)
% Compute QCM forward model
%
% Synopsis
%    [responses] = tfeQCMForward(params,stimuli)
%
% Description:
%    Take stimuli and compute QCM responses.
%
% Inputs:
%      params        - QCM model parameter struct
%      stimuli       - Stimuli, with stimulus contrasts in columns
% 
% Outputs:
%      responses     - Row vector of responses

% History:
%   11/24/19  dhb    Wrote it for modularity.
%             dhb    Enforce constraint on minor axis.

%% Check parameterization
dimension = size(stimuli,1);
if (dimension == 2)
    if (params.Qvec(1) > 1)
        error('Minor axis must have length less than or equal to 1');
    end
elseif (dimension == 3)
    if (params.Qvec(1) > 1)
        error('First minor axis must have length less than or equal to 1');
    end
        if (params.Qvec(2) > params.Qvec(2))
        error('Second minor axis must have length less than or equal to that of first minor axis');
    end
else
    error('Can only handle dimension of 2 or 3 for stimuli');
end

%% Get the ellipsoid parameters in cannonical form
[~,~,Q] = EllipsoidMatricesGenerate([1 params.Qvec]','dimension',dimension);

%% Find the length of the points after application of the quadratic
%
% This represents the quadaratic component of the neural response after
% application of the quadratic
theLengths = diag(sqrt(stimuli'*Q*stimuli))';

%% Push the quadratic response through a Naka-Rushton non-linearity
responses = ComputeNakaRushton([params.crfAmp,params.crfSemi,params.crfExponent],theLengths) + params.offset;


