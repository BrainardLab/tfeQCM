function [indDirectionNRParams,indDirectionPredictions,indDirectionResponses,indDirectionDirections,indDirectionContrasts,indDirectionIndices,nIndDirections] = ...
    tfeQCMFitNakaRushtonDirectionsContrasts(responses,directions,contrasts,varargin)
% Fit Naka-Rushton functions to contrast-response functions measured in a set of directions
%
% Synopsis
%    [indDirectionNRParams,indDirectionPredictions,indDirectionResponses,indDirectionDirections,indDirectionContrasts,indDirectionIndices,nIndDirecitons] = ...
%        tfeQCMFitNakaRushtonDirectionsContrasts(responses,directions,contrasts)
%
% Description:
%    Take stimuli in form where we have the direction and contrast of each stimulus, 
%    and fit Naka-Rushton functions to the contrast response functions in
%    each direction.
% 
%    This parses the stimuli into a set of contrast response functions in
%    individual directions and then fits each of these individual
%    directions.
%
%    Can place constraints that keep some parameters fixed across
%    direcitons.
%
% Inputs:
%      responses               - Row vector of responses measured for each
%                                stimulus.
%      directions              - Matrix of directions for each stimulus, with
%                                directions in columns. Make sure columns
%                                match exactly for directions that are
%                                supposed to be the same.
%      contrasts               - Row vector of contrasts corresponding to columns
%                                of directions.
% 
% Outputs:
%      indDirectionNRParams    - Struct array. The Naka-Rushton parameters for each direction.
%      indDirectionPredictions - Cell array. The predictions of responses for each direction.
%      indDirectionResponses   - Cell array. The responses for each direction.
%      indDirectionDirections  - Array. The unique directions are in the
%                                columns of this array.
%      indDirectionContrasts   - Cell array. The contrasts for each direction.
%      indDirectionIndices     - Cell array. For each direction, the
%                                indices into the input data corresponding to that
%                                direction.
%      nIndDirections          - Number of independent directions
%
% Optional key/value pairs:
%      'lockOffsetToZero'      - Logical (default false). Force fits to go through 0 at 0 contrast
%      'commonAmplitude'       - Logical (default false). Force common amplitude across directions.
%      'commonSemi'            - Logical (default false). Force common semi-saturation across directions.
%      'commonExp'             - Logical (default false). Force common exponent across directions.
%      'commonOffset'          - Logical (default false). Force common offset across directions.
%
% See also:
%

% History:
%   11/24/18  dhb    Wrote it.

%% Parse key/value pairs
% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = false;
p.addParameter('lockOffsetToZero',false,@islogical);
p.addParameter('commonAmplitude',false,@islogical);
p.addParameter('commonSemi',false,@islogical);
p.addParameter('commonExp',false,@islogical);
p.addParameter('commonOffset',false,@islogical);
p.parse(varargin{:});

%% Parameters
unique_tolerance = 1e-6;

%% Parse the stimuli into indDirections, contrasts and responses for each individual direction
% 
% The commented out call to uniquetol is almost good, but scrambles the
% order of the unique outputs.  One could fix that with a little work, if
% tolerance turns out to be an issue.
[indDirectionDirections,indDirectionIndices] = tfeQCMParseDirections(directions);
nIndDirections = size(indDirectionDirections,2);
for ii = 1:nIndDirections
    indDirectionResponses{ii} = responses(indDirectionIndices{ii});
    indDirectionContrasts{ii} = contrasts(indDirectionIndices{ii});
end

%% Initialize parameters
indDirectionNRParams0 = tfeNRInitializeParams(nIndDirections);
indDirectionParamsvec0 = tfeNRParamsToVec(indDirectionNRParams0);

%% Set up search bounds
ampLowBound = 0; ampHighBound = 5;
semiLowBound = 0.01; semiHighBound = 10;
expLowBound = 0.01; expHighBound = 10;
expFalloffLowBound = indDirectionNRParams0.expFalloff;
expFalloffHighBound = indDirectionNRParams0.expFalloff;
noiseSdLowBound = indDirectionNRParams0.noiseSd;
noiseSdHighBound = indDirectionNRParams0.noiseSd;
if (p.Results.lockOffsetToZero)
    offsetLowBound = 0;
    offsetHighBound = 0;
else
    offsetLowBound = -ampHighBound;
    offsetHighBound = ampHighBound;
end

% Pack bounds into vector form of parameters.
for ii = 1:nIndDirections
    NRParamsLow(ii).crfAmp = ampLowBound;
    NRParamsLow(ii).crfSemi = semiLowBound;
    NRParamsLow(ii).crfExponent = expLowBound;
    NRParamsLow(ii).crfOffset = offsetLowBound;
    NRParamsLow(ii).expFalloff = expFalloffLowBound;
    NRParamsLow(ii).noiseSd = noiseSdLowBound;
end
vlb = tfeNRParamsToVec(NRParamsLow);
for ii = 1:nIndDirections
    NRParamsHigh(ii).crfAmp = ampHighBound;
    NRParamsHigh(ii).crfSemi = semiHighBound;
    NRParamsHigh(ii).crfExponent = expHighBound;
    NRParamsHigh(ii).crfOffset = offsetHighBound;
    NRParamsHigh(ii).expFalloff = expFalloffHighBound;
    NRParamsHigh(ii).noiseSd = noiseSdHighBound;
end
vub = tfeNRParamsToVec(NRParamsHigh);

%% Set up linear parameter constraints
%
% This is a little tricky. We need to know the order
% of the parameters in the parameter vector, and then
% we set 1,-1 pairs between the entries for the first
% function and those for each subsequent function.
%
% I don't see any easy way around knowing the order in
% which the parameters are packed into vectors here.
[~,nParams] = tfeNRVecToParams(vlb,nIndDirections);
if (nIndDirections > 1)
    % Build constraint matrix if there is more than one direction.
    %
    % Initialize
    Aeq = [];
    beq = [];
    eqRowIndex = 1;
    
    % Common amplitude constraints
    if (p.Results.commonAmplitude)
        paramIndex = 1;
        for ii = 2:nIndDirections
            Aeq(eqRowIndex,:) = zeros(size(vlb'));
            Aeq(eqRowIndex,paramIndex) = 1;
            Aeq(eqRowIndex,(ii-1)*nParams+paramIndex) = -1;
            beq(eqRowIndex) = 0;
            eqRowIndex = eqRowIndex+1;
        end
    end
    
    % Common semi-saturation constant constraints
    if (p.Results.commonSemi)
        paramIndex = 2;
        for ii = 2:nIndDirections
            Aeq(eqRowIndex,:) = zeros(size(vlb'));
            Aeq(eqRowIndex,paramIndex) = 1;
            Aeq(eqRowIndex,(ii-1)*nParams+paramIndex) = -1;
            beq(eqRowIndex) = 0;
            eqRowIndex = eqRowIndex+1;
        end
    end
    
    % Common exponent constraints
    if (p.Results.commonExp)
        paramIndex = 3;
        for ii = 2:nIndDirections
            Aeq(eqRowIndex,:) = zeros(size(vlb'));
            Aeq(eqRowIndex,paramIndex) = 1;
            Aeq(eqRowIndex,(ii-1)*nParams+paramIndex) = -1;
            beq(eqRowIndex) = 0;
            eqRowIndex = eqRowIndex+1;
        end
    end
    
    % Common offset constraints
    if (p.Results.commonOffset)
        paramIndex = 4;
        for ii = 2:nIndDirections
            Aeq(eqRowIndex,:) = zeros(size(vlb'));
            Aeq(eqRowIndex,paramIndex) = 1;
            Aeq(eqRowIndex,(ii-1)*nParams+paramIndex) = -1;
            beq(eqRowIndex) = 0;
            eqRowIndex = eqRowIndex+1;
        end
    end
else
    % If there is only one independent direction, then there is no
    % constraint across directions to be had.
    Aeq = [];
    beq = [];
end
    

%% Fit using fmincon
options = optimset;
options = optimset(options,'Diagnostics','off','Display','off');
options = optimset(options,'LargeScale','off');
indDirectionParamsvec = fmincon(@(x)FitIndNakaRushtonFun(x,indDirectionResponses,indDirectionContrasts,nIndDirections), ...
    indDirectionParamsvec0,[],[],Aeq,beq,vlb,vub,[],options); 

[~,indDirectionPredictions] = FitIndNakaRushtonFun(indDirectionParamsvec,indDirectionResponses,indDirectionContrasts,nIndDirections);
indDirectionNRParams = tfeNRVecToParams(indDirectionParamsvec,nIndDirections);


end

% This is matched to what the object based error functions do.
function [f,indDirectionPredictions] = FitIndNakaRushtonFun(paramsvec,indDirectionResponses,indDirectionContrasts,nIndDirections)

% Convert paramsvec to cell array
NRParams = tfeNRVecToParams(paramsvec,nIndDirections);

% Loop over directions
diff = [];
for ii = 1:nIndDirections
    indDirectionPredictions{ii} = ComputeNakaRushton([NRParams(ii).crfAmp,NRParams(ii).crfSemi,NRParams(ii).crfExponent],indDirectionContrasts{ii}) + NRParams(ii).crfOffset;
    diff = [diff indDirectionResponses{ii}-indDirectionPredictions{ii}];
end
f = 1000*sqrt(nanmean(diff.^2));
end



