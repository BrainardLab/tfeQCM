function [paramsFit,fVal,modelResponseStruct] = fitResponse(obj,thePacket,varargin)
% [paramsFit,fVal,modelResponseStruct] = fitResponse(obj,thePacket,varargin)
%
% Fit method for the tfeNakaRushtonDirection class.  This overrides the tfeQCM
% method, which we need to do because that is the direct parent class of
% tfeNakaRushtonDirection.
%
% Inputs:
%   thePacket          - A valid packet
%
% Outputs:
%   paramsFit          - Fit parameters
%   fVal               - Fit error.
%   predictedResponse  - Response predicted from fit
%
% Optional key/value pairs
%   See tfe.fitResponse for these.

%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('thePacket',@isstruct);
p.addParameter('defaultParamsInfo',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('defaultParams',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('searchMethod','fmincon',@ischar);
p.addParameter('DiffMinChange',[],@isnumeric);
p.addParameter('fminconAlgorithm','active-set',@ischar);
p.addParameter('errorType','rmse',@ischar);
p.parse(thePacket,varargin{:});

%% Check packet validity
if (~obj.isPacket(thePacket))
    error('The passed packet is not valid for this model');
else
    switch (obj.verbosity)
        case 'high'
            fprintf('valid\n');
    end
end

%% Set initial values and reasonable bounds on parameters
[indDirectionNRParams0,NRParamsLow,NRParamsHigh] = obj.defaultParams('defaultParamsInfo',p.Results.defaultParamsInfo,'defaultParams',p.Results.defaultParams,varargin{:});
paramsFitVec0 = obj.paramsToVec(indDirectionNRParams0);

%% Set up search bounds
ampLowBound = -5; ampHighBound = 5;
semiLowBound = 0.01; semiHighBound = 10;
expLowBound = 0.01; expHighBound = 10;
expFalloffLowBound = indDirectionNRParams0.expFalloff;
expFalloffHighBound = indDirectionNRParams0.expFalloff;
noiseSdLowBound = indDirectionNRParams0.noiseSd;
noiseSdHighBound = indDirectionNRParams0.noiseSd;
if (obj.lockOffsetToZero)
    offsetLowBound = 0;
    offsetHighBound = 0;
else
    offsetLowBound = -ampHighBound;
    offsetHighBound = ampHighBound;
end

% Pack bounds into vector form of parameters.
for ii = 1:obj.nDirections
    NRParamsLow(ii).crfAmp = ampLowBound;
    NRParamsLow(ii).crfSemi = semiLowBound;
    NRParamsLow(ii).crfExponent = expLowBound;
    NRParamsLow(ii).crfOffset = offsetLowBound;
    NRParamsLow(ii).expFalloff = expFalloffLowBound;
    NRParamsLow(ii).noiseSd = noiseSdLowBound;
end
vlb = obj.paramsToVec(NRParamsLow);
for ii = 1:obj.nDirections
    NRParamsHigh(ii).crfAmp = ampHighBound;
    NRParamsHigh(ii).crfSemi = semiHighBound;
    NRParamsHigh(ii).crfExponent = expHighBound;
    NRParamsHigh(ii).crfOffset = offsetHighBound;
    NRParamsHigh(ii).expFalloff = expFalloffHighBound;
    NRParamsHigh(ii).noiseSd = noiseSdHighBound;
end
vub = obj.paramsToVec(NRParamsHigh);

%% Set up linear parameter constraints
%
% This is a little tricky. We need to know the order
% of the parameters in the parameter vector, and then
% we set 1,-1 pairs between the entries for the first
% function and those for each subsequent function.
%
% I don't see any easy way around knowing the order in
% which the parameters are packed into vectors here.
tempvec = obj.paramsToVec(NRParamsLow(1));
nParams = length(tempvec);
if (obj.nDirections > 1)
    % Build constraint matrix if there is more than one direction.
    %
    % Initialize
    Aeq = [];
    beq = [];
    eqRowIndex = 1;
    
    % Common amplitude constraints
    if (obj.commonAmplitude)
        paramIndex = 1;
        for ii = 2:obj.nDirections
            Aeq(eqRowIndex,:) = zeros(1,nParams*length(NRParamsLow'));
            Aeq(eqRowIndex,paramIndex) = 1;
            Aeq(eqRowIndex,(ii-1)*nParams+paramIndex) = -1;
            beq(eqRowIndex) = 0;
            eqRowIndex = eqRowIndex+1;
        end
    end
    
    % Common semi-saturation constant constraints
    if (obj.commonSemi)
        paramIndex = 2;
        for ii = 2:obj.nDirections
            Aeq(eqRowIndex,:) = zeros(1,nParams*length(NRParamsLow'));
            Aeq(eqRowIndex,paramIndex) = 1;
            Aeq(eqRowIndex,(ii-1)*nParams+paramIndex) = -1;
            beq(eqRowIndex) = 0;
            eqRowIndex = eqRowIndex+1;
        end
    end
    
    % Common exponent constraints
    if (obj.commonExp)
        paramIndex = 3;
        for ii = 2:obj.nDirections
            Aeq(eqRowIndex,:) = zeros(1,nParams*length(NRParamsLow'));
            Aeq(eqRowIndex,paramIndex) = 1;
            Aeq(eqRowIndex,(ii-1)*nParams+paramIndex) = -1;
            beq(eqRowIndex) = 0;
            eqRowIndex = eqRowIndex+1;
        end
    end
    
    % Common offset constraints
    if (obj.commonOffset)
        paramIndex = 4;
        for ii = 2:obj.nDirections
            Aeq(eqRowIndex,:) = zeros(1,nParams*length(NRParamsLow'));
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

%% David sez: "Fit that sucker"
switch (obj.verbosity)
    case 'high'
        fprintf('Fitting.');
end

%% fmincon fit
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','iter','LargeScale','off','Algorithm',p.Results.fminconAlgorithm);
options = optimset(options,'TolCon',1e-3);
if ~isempty(p.Results.DiffMinChange)
    options = optimset(options,'DiffMinChange',p.Results.DiffMinChange);
end
paramsFitVec = fmincon(@(modelParamsVec)obj.fitError(modelParamsVec, ...
    thePacket),paramsFitVec0,[],[],Aeq,beq,vlb,vub,[],options);

% Get error and predicted response for final parameters
[fVal,modelResponseStruct] = obj.fitError(paramsFitVec,thePacket,'errorType',p.Results.errorType);

switch (obj.verbosity)
    case 'high'
        fprintf('\n');
        fprintf('Fit error value: %g', fVal);
        fprintf('\n');
end

% Convert fit parameters for return
paramsFit = obj.vecToParams(paramsFitVec);

end
        

