function [isoContrast,unitEquivContrastResponse,angleSupport] = getIsoContrast(obj,params,varargin)
% [isoResponse,unitResponse,angleSupport] = getIsoContrast(obj,params,varargin)
%
% Synopsis:
%   Get isoresponse contour for the LCM (linear channel model).
%
% Inputs:
%   obj        - the tfeLCMDirection object
%   params     - parameter structure for channel being evaluated
%
% Outputs:
%   isoContrast - Contrast as a function of angle support that produces
%                 criterion response (obj.criterionResp)
%   unitContrastResponse - Reponse as a function of angle support to unit
%                 contrast input.
%   angleSupport - the angle support in degrees (almost surely 1:360 in 1 
%                 degree increments.)
%
% Optional key/value pairs
%   'None.

%% History
%    03/01/21  dhb  Wrote.

%% Parse input
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true;
p.addRequired('params',@isstruct);
p.parse(params,varargin{:});
params = p.Results.params;

% Create a single linear channel as a weighted combination
% of the underlying channels. 
%
% Since we use symmetric modulations in our experiments, we'll keep
% this symmetric.  We specify the first three weights and then duplicat
% for the second three.
theChannelWeights = [[params.channelWeightsPos]' ; [params.channelWeightsPos]'];
theChannel = (obj.underlyingChannels'*theChannelWeights)';

% Response around circle to unit contrast is just the channel sensitivity
unitEquivContrastResponse1 = theChannel;

% Compute linear response (no NR) using compute method
if (obj.dimension ~= 2)
    error('This only works if stimlus dimension is 2');
end
directionStimulusStruct.timebase = 1:length(obj.angleSupport);
directionStimulusStruct.values(1,:) = cosd(obj.angleSupport);
directionStimulusStruct.values(2,:) = sind(obj.angleSupport);
directionStimulusStruct.values(3,:) = ones(size(obj.angleSupport)); 
unitEquivContrastResponseStruct = obj.computeResponse(params,directionStimulusStruct,[],'noNakaRushton',true);
unitEquivContrastResponse2 = unitEquivContrastResponseStruct.values;

% Check
if (obj.summationExponent == 1)
    if (max(abs(unitEquivContrastResponse1-unitEquivContrastResponse2)) > 1e-10)
        error('Do not get same answer twice');
    end
    unitEquivContrastResponse = unitEquivContrastResponse1;
else
    unitEquivContrastResponse = unitEquivContrastResponse2;
end

% Deal with small value problem
unitEquivContrastResponse(unitEquivContrastResponse <= 1e-6) = 1e-6;

% Invert Naka-Rushton to get the equivalent contrast corresponding to
% the criterion response
criterionEquivContrastCell = tfeNRInvert(params,{obj.criterionResp});
criterionEquivContrast = criterionEquivContrastCell{1};

exponentiatedContrast = criterionEquivContrast./unitEquivContrastResponse;
isoContrast = exponentiatedContrast.^(1/obj.summationExponent);

% Get isocontrast around unit circle to produce criterion response.
%isoContrast = obj.criterionResp./(criterionEquivContrast.^(1/obj.summationExponent));

% Check that isoContrast really does produce an isoresponse.
directionStimulusStruct.values(3,:) = isoContrast;
checkContrastResponseStruct = obj.computeResponse(params,directionStimulusStruct,[],'noNakaRushton',false);
checkContrastResponse = checkContrastResponseStruct.values;
if (max(abs(checkContrastResponse - obj.criterionResp)) > 1e-8)
    error('iso contrast contour does not actually produce iso response');
end

% Set angle support for return
angleSupport = obj.angleSupport;
if (any(1:1:360 ~= angleSupport))
    error('Someone changed angle support property, make sure this was done carefully')
end
    

end

