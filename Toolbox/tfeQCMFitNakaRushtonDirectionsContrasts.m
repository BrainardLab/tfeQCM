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
%      indDirectionNRParams    - Cell array. The Naka-Rushton parameters for each direction.
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
%
% See also:
%

% History:
%   11/24/19  dhb    Wrote it.

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

%% Set up bounds. Could be a little slicker, maybe
vlbRaw = [0   0.1  0.01  -5]';
vubRaw = [5   10   100    5]';
vlb = [];
vub = [];
for ii = 1:nIndDirections
    vlb = [vlb ; vlbRaw];
    vub = [vub ; vubRaw];
end

%% Fit using fmincon
options = optimset;
options = optimset(options,'Diagnostics','off','Display','iter');
options = optimset(options,'LargeScale','off');
indDirectionParamsvec = fmincon(@(x)FitIndNakaRushtonFun(x,indDirectionResponses,indDirectionDirections,indDirectionContrasts,nIndDirections), ...
    indDirectionParamsvec0,[],[],[],[],vlb,vub,[],options); 

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
    indDirectionPredictions{ii} = ComputeNakaRushton([NRParams{ii}(1),NRParams{ii}(2),NRParams{ii}(3)],indDirectionContrasts{ii}) + NRParams{ii}(4);
    diff = indDirectionResponses{ii}-indDirectionPredictions{ii};
    f = f + sum(diff.^2);
end

end

function paramsvec = NRParamsToVec(NRParams)
    nIndDirections = length(NRParams);
    paramsvec = [];
    for ii = 1:nIndDirections
        paramsvec = [paramsvec ; NRParams{ii}];
    end
end

function NRParams = NRVecToParams(paramsvec,nIndDirections)
    nParams = 4;
    NRParams = cell(1,nIndDirections);
    for ii = 1:nIndDirections
        NRParams{ii} = paramsvec((ii-1)*nParams+1:ii*nParams);
    end
end

function NRParams = InitializeNRParams(nIndDirections)

% Reasonable defaults
nParams = 4;
Rmax = 1;
sigma = 0.06;
n = 2;
offset = 0;

% Set up the cell array
for ii = 1:nIndDirections
    NRParams{ii} = zeros(nParams,1);
    NRParams{ii}(1) = Rmax;
    NRParams{ii}(2) = sigma;
    NRParams{ii}(3) = n;
    NRParams{ii}(4) = offset;
end

end


