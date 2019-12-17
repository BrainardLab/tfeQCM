function [responses,quadraticFactors] = tfeQCMForward(params,stimuli)
% Compute QCM forward model
%
% Synopsis
%    [responses,quadraticFactors] = tfeQCMForward(params,stimuli)
%
% Description:
%    Take stimuli and compute QCM responses.
%
%    Computation method changes depending on number of stimuli,
%    as the fastest way changes with that variable as well. 
%    Currently coded with break point at 500, but you can
%    change in the source by changing value of variable methodSwitchN.
%    Could recode so that this is an optional parameter with a default.
%
% Inputs:
%      params        - QCM model parameter struct
%      stimuli       - Stimuli, with stimulus contrasts in columns
% 
% Outputs:
%      responses     - Row vector of responses
%      quadraticFactors - Vector with same length as number of stimuli.
%                      Divide the semi saturation constant of the
%                      Naka-Rushton by this if you want Naka-Rushton
%                      parameters to apply directly to the stimuli. (The
%                      parameters passed to this routine are applied after
%                      the quadratic structure has been taken into
%                      account.)

% History:
%   11/24/19  dhb    Wrote it for modularity.
%             dhb    Enforce constraint on minor axis.

%% Method switch N
methodSwitchN = 500;

%% Check parameterization
if (~tfeQCMCheckParams(params))
    error('Bad values in QCM parameters structure');
end

%% Get the ellipsoid parameters in cannonical form
dimension = size(stimuli,1);
[~,~,Q] = EllipsoidMatricesGenerate([1 params.Qvec]','dimension',dimension);

%% Find the length of the points after application of the quadratic
%
% This represents the quadaratic component of the neural response after
% application of the quadratic
%
% Which way is faster depends on number of time points.  500
% (default value of methodSwitchN) is a reasonable guess as the
% right point to switch over to the loop.
if (size(stimuli,2) < methodSwitchN)
    equivalentContrasts = sqrt(diag(stimuli'*Q*stimuli))';
else
    equivalentContrasts = zeros(1,size(stimuli,2));
    for ii = 1:size(stimuli,2)
        equivalentContrasts(ii) = sqrt(stimuli(:,ii)'*Q*stimuli(:,ii));
    end
end

% Only do this if we need it for the return variable.  Computing norm
% is a little faster than sqrt of the explicit dot product in the
% stimulus expression, at least for the one case I tested.
if ((isfield(params,'fitting') & params.fitting) | nargout < 2)
    quadraticFactors = [];
else
    if (size(stimuli,2) < methodSwitchN)
        stimulusContrasts = sqrt(diag(stimuli'*stimuli))';
    else
        stimulusContrasts = zeros(1,size(stimuli,2));
        for ii = 1:size(stimuli,2)
            stimulusContrasts(ii) = norm(stimuli(:,ii));
        end
    end
    quadraticFactors = equivalentContrasts ./ stimulusContrasts;
end

%% Push the quadratic response through a Naka-Rushton non-linearity
responses = ComputeNakaRushton([params.crfAmp,params.crfSemi,params.crfExponent],equivalentContrasts) + params.crfOffset;

%% Conceptually
% responses = ComputeNakaRushton([params.crfAmp,params.crfSemi,params.crfExponent],quadraticFactors.*stimulusContrasts) + params.crfOffset;
% responses = ComputeNakaRushton([params.crfAmp,params.crfSemi ./ quadraticFactors,params.crfExponent],stimulusContrasts) + params.crfOffset;


