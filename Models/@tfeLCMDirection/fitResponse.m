function [paramsFit,fVal,modelResponseStruct] = fitResponse(obj,thePacket,varargin)
% Fit LCM model to a packet.
%
% Synopsis:
%    [paramsFit,fVal,modelResponseStruct] = fitResponse(obj,thePacket,varargin)
%
% Description:
%    Fit method for the tfeLCMDirection class.  This overrides the tfe method, allowing
%    a certain amount of customizaiton.
%
%    Also eventually will allow us to put on some parameter constraints that
%    are not generic.
%
%    Currently constrains thing by using vlb and vub passed to fmincon that
%    keep the first channel coefficient at 1, implemented by the default
%    parameters.  A better method would be to constrain the vector length of
%    the channel weights, but this is good enough for our current purposes.
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
%  Those understood by tfe.fitResponse with the following adjustments.
%
% 'fminconAlgorithm'     - String (default 'interior-point'). If set to a string,
%                           passed on as algorithm in options to fmincon.
%                           Can be empty or any algorithm string understood
%                           by fmincon.
%                              [] - Use fmincon's current default algorithm
%                              'active-set' - Active set algorithm
%                              'interior-point' - Interior point algorithm.
%  'fitErrorScalar'       - Computed fit error is multiplied by this before
%                           return.  Sometimes getting the objective
%                           function onto the right scale makes all the
%                           difference in fitting. Passed along as an
%                           option to the fitError method, but overrides
%                           the fitError's default value, Default here is
%                           1000.
%   'noNakaRushton'       - true/false (default false). If true, don't apply
%                           Naka-Rushton, just return equivalent contrast.
%
%

%% History
%    03/01/21  dhb  Wrote from QCM version.
%    04/03/21  dhb  New normalization, and enforce positivity on responses
%                   around angle.

%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('thePacket',@isstruct);
p.addParameter('initialParams',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('fitErrorScalar',1000,@isnumeric);
p.addParameter('maxIter',[],@isnumeric);
p.addParameter('maxFunEval',[],@isnumeric);
p.addParameter('noNakaRushton',false,@islogical);
p.addParameter('fminconAlgorithm','interior-point',@(x) (isempty(x) | ischar(x)));
p.parse(thePacket,varargin{:});

%% Initial parameters
[initialParams,vlbParams,vubParams] = obj.defaultParams;
if (~isempty(p.Results.initialParams))
    initialParams = p.Results.initialParams;
end

% Normalize the initial parameters to norm to 1
initialVec = obj.paramsToVec(initialParams);
initialVec(1:obj.nChannels/2) = initialVec(obj.nChannels/2)/norm(initialVec(1:obj.nChannels/2));
initialParams = obj.vecToParams(initialVec);

%% Locked Naka-Rushton parameters
if (~isempty(obj.lockedCrfAmp))
    initialParams.crfAmp = obj.lockedCrfAmp;
    vlbParams.crfAmp = obj.lockedCrfAmp;
    vubParams.crfAmp = obj.lockedCrfAmp;
end
if (~isempty(obj.lockedCrfExponent))
    initialParams.crfExponent = obj.lockedCrfExponent;
    vlbParams.crfExponent = obj.lockedCrfExponent;
    vubParams.crfExponent = obj.lockedCrfExponent;
end
if (~isempty(obj.lockedCrfSemi))
    initialParams.crfSemi = obj.lockedCrfSemi;
    vlbParams.crfSemi = obj.lockedCrfSemi;
    vubParams.crfSemi = obj.lockedCrfSemi;
end
if (~isempty(obj.lockedCrfOffset))
    initialParams.crfOffset = obj.lockedCrfOffset;
    vlbParams.crfOffset = obj.lockedCrfOffset;
    vubParams.crfOffset= obj.lockedCrfOffset;
end

% Setting the fitting flag allows the eventually called
% computeResponse routine to take some shortcuts that we
% don't want taken in general.  But we care enough about execution
% speed to do this.  Be sure to set stimuli back to empty when setting
% flag back to false.  The computeResponse routine will compute it when
% the fitting flag is true and it is empty.
obj.fitting = true; obj.angles = [];
[paramsFit1,fVal1,modelResponseStruct1] = fitResponse@tfe(obj,thePacket,...
    'initialParams',initialParams,'vlbParams',vlbParams,'vubParams',vubParams,...
    'fitErrorScalar',p.Results.fitErrorScalar,...
    'noNakaRushton',p.Results.noNakaRushton, ...
    'MaxIter',p.Results.maxIter, ...
    'MaxFunEval',p.Results.maxFunEval, ...
    'nlcon',@(x)fitNlCon(x,obj));
obj.fitting = false;
obj.angles = [];

% Set return
paramsFit = paramsFit1;
fVal = fVal1;
modelResponseStruct = modelResponseStruct1;

% Sanity check. This is not needed with new normalization procedure.
% for ii = 2:length(paramsFit.channelWeightsPos)/2
%     checkVal = paramsFit.channelWeightsPos(ii)/paramsFit.channelWeightsPos(1);
%     if (checkVal < 1e-4 | checkVal > 1e4)
%         error('Large asymmetry in fit channel weights. May be a symptom that locking first weight is not appropriate for this data set');
%     end
% end

% Use this to check error value as it sits here
fValCheck = obj.fitError(obj.paramsToVec(paramsFit),thePacket,varargin{:},'fitErrorScalar',p.Results.fitErrorScalar,'noNakaRushton',p.Results.noNakaRushton);
if (fValCheck ~= fVal)
    error('Cannot compute the same fit error twice the same way. Check.');
end


end

function [C,Ceq] = fitNlCon(x,obj)

% Constraint on vector length
Ceq = norm(x(1:obj.nChannels/2))-1;

% Keep responses around the circle positive
if (obj.dimension ~= 2)
    error('This only works if stimlus dimension is 2');
end
stimulusStruct.timebase = 1:length(obj.angleSupport);
stimulusStruct.values(1,:) = cosd(obj.angleSupport);
stimulusStruct.values(2,:) = sind(obj.angleSupport);
stimulusStruct.values(3,:) = ones(size(obj.angleSupport)); 
unitContrastResponseStruct = obj.computeResponse(obj.vecToParams(x),stimulusStruct,[],'noNakaRushton',true);
C = -min(unitContrastResponseStruct.values);

end





