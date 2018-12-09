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
%      indDirectionDirections  - Cell array. The direction for each direction.
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
[indDirectionsTemp,~,whichColumnsOut] = unique(directions','rows','stable');  % uniquetol(directions',unique_tolerance,'ByRows',true);
indDirectionsTemp = indDirectionsTemp';
nIndDirections = size(indDirectionsTemp,2);
for ii = 1:nIndDirections
    whichColumns = find(whichColumnsOut == ii);
    indDirectionIndices{ii} = whichColumns;
    indDirectionResponses{ii} = responses(whichColumns);
    indDirectionDirections{ii} = indDirectionsTemp(:,ii);
    indDirectionContrasts{ii} = contrasts(whichColumns);
end

%% Initialize parameters
indDirectionNRParams0 = InitializeNRParams(nIndDirections);
indDirectionParamsvec0 = NRParamsToVec(indDirectionNRParams0);

%% Set up search bounds
ampLowBound = 0; ampHighBound = 5;
semiLowBound = 0.01; semiHighBound = 10;
expLowBound = 0.01; expHighBound = 10;
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
end
vlb = NRParamsToVec(NRParamsLow);
for ii = 1:nIndDirections
    NRParamsHigh(ii).crfAmp = ampHighBound;
    NRParamsHigh(ii).crfSemi = semiHighBound;
    NRParamsHigh(ii).crfExponent = expHighBound;
    NRParamsHigh(ii).crfOffset = offsetHighBound;
end
vub = NRParamsToVec(NRParamsHigh);

%% Set up linear parameter constraints
%
% This is a little tricky. We need to know the order
% of the parameters in the parameter vector, and then
% we set 1,-1 pairs between the entries for the first
% function and those for each subsequent function.
%
% I don't see any easy way around knowing the order in
% which the parameters are packed into vectors here.
[~,nParams] = NRVecToParams(vlb,nIndDirections);
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
indDirectionParamsvec = fmincon(@(x)FitIndNakaRushtonFun(x,indDirectionResponses,indDirectionDirections,indDirectionContrasts,nIndDirections), ...
    indDirectionParamsvec0,[],[],Aeq,beq,vlb,vub,[],options); 

[~,indDirectionPredictions] = FitIndNakaRushtonFun(indDirectionParamsvec,indDirectionResponses,indDirectionDirections,indDirectionContrasts,nIndDirections);
indDirectionNRParams = NRVecToParams(indDirectionParamsvec,nIndDirections);


end

function [f,indDirectionPredictions] = FitIndNakaRushtonFun(paramsvec,indDirectionResponses,indDirectionDirections,indDirectionContrasts,nIndDirections)

% Initialize error message
f = 0;

% Convert paramsvec to cell array
NRParams = NRVecToParams(paramsvec,nIndDirections);

% Loop over directions
for ii = 1:nIndDirections
    indDirectionPredictions{ii} = ComputetfeQCMComputeNakaRushton([NRParams(ii).crfAmp,NRParams(ii).crfSemi,NRParams(ii).crfExponent],indDirectionContrasts{ii}) + NRParams(ii).crfOffset;
    diff = indDirectionResponses{ii}-indDirectionPredictions{ii};
    f = f + sum(diff.^2);
end

end

function paramsvec = NRParamsToVec(NRParams)
% Convert the NR parameter structure array to one long vector
%
% See also: NRVecToParams, InitializeNRParams
%

    nIndDirections = length(NRParams);
    paramsvec = [];
    for ii = 1:nIndDirections
        paramsvec = [paramsvec ; [NRParams(ii).crfAmp NRParams(ii).crfSemi NRParams(ii).crfExponent NRParams(ii).crfOffset]'];
    end
end

function [NRParams,nParams] = NRVecToParams(paramsvec,nIndDirections)
% Convert the NR parameter vector to struct array
%
% See also: NRParamsToVec, InitializeNRParams
%
    
    % Number of parameters.  You just have to know this.  Coordinate
    % any changes with the NRParamsToVec routine.
    nParams = 4;
    
    % Unpack vector into struct array.
    for ii = 1:nIndDirections        
       NRParams(ii).crfAmp =  paramsvec((ii-1)*nParams+1);
       NRParams(ii).crfSemi = paramsvec((ii-1)*nParams+2);
       NRParams(ii).crfExponent = paramsvec((ii-1)*nParams+3);
       NRParams(ii).crfOffset = paramsvec((ii-1)*nParams+4);
    end
end

function NRParams = InitializeNRParams(nIndDirections)
% Initialize the NR parameter struct array
%
% See also: NRParamsToVec, NRVecToParams
%

% Reasonable defaults
nParams = 4;
Rmax = 1;
sigma = 0.06;
n = 2;
offset = 0;

% Set up the struct array
for ii = 1:nIndDirections
    NRParams(ii).crfAmp = Rmax;
    NRParams(ii).crfSemi = sigma;
    NRParams(ii).crfExponent = n;
    NRParams(ii).crfOffset = offset;
end

end


