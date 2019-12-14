function [responses,quadraticFactors] = tfeQCMForward(params,stimuli)
% Compute QCM forward model
%
% Synopsis
%    [responses,quadraticFactors] = tfeQCMForward(params,stimuli)
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
% Although this is conceptually what we want to do, it is
% much slower than the loop version when stimulus size gets
% large, because the loop does not compute the off diagonal
% elements that we then immediately toss. For a stimulus
% with 7200 time points, the loop is over a factor 10 faster.
% equivalentContrasts = sqrt(diag(stimuli'*Q*stimuli))';
equivalentContrasts = zeros(1,size(stimuli,2));
for ii = 1:size(stimuli,2)
  equivalentContrasts(ii) = sqrt(stimuli(:,ii)'*Q*stimuli(:,ii));  
end   

% Only do this if we need it for the return variable.  Computing norm
% is a little faster than sqrt of the explicit dot product.
if (nargout >= 2)
    stimulusContrasts = zeros(1,size(stimuli,2));
    for ii = 1:size(stimuli,2)
        stimulusContrasts(ii) = norm(stimuli(:,ii));  
    end
    quadraticFactors = equivalentContrasts ./ stimulusContrasts;
end

%% Push the quadratic response through a Naka-Rushton non-linearity
responses = ComputeNakaRushton([params.crfAmp,params.crfSemi,params.crfExponent],equivalentContrasts) + params.crfOffset;

%% Conceptually
% responses = ComputeNakaRushton([params.crfAmp,params.crfSemi,params.crfExponent],quadraticFactors.*stimulusContrasts) + params.crfOffset;
% responses = ComputeNakaRushton([params.crfAmp,params.crfSemi ./ quadraticFactors,params.crfExponent],stimulusContrasts) + params.crfOffset;


