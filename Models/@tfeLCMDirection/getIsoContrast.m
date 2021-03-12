function [isoContrast,unitContrastResponse,angleSupport] = getIsoContrast(obj,params,varargin)
% [isoResponse,unitResponse,angleSupport] = getIsoContrast(obj,params,varargin)
%
% Get isoresponse contour for the LCM (linear channel model).
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
unitContrastResponse = theChannel;

% Deal with small value problem
unitContrastResponse(unitContrastResponse <= 1e-6) = 1e-6;

% Get isocontrast around unit circle to produce criterion response.
isoContrast = obj.criterionResp./unitContrastResponse;

% Set angle support for return
angleSupport = obj.angleSupport;
if (any(1:1:360 ~= angleSupport))
    error('Someone changed angle support property, make sure this was done carefully')
end
    

end

