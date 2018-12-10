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

if (isempty(p.Results.defaultParams))
    defaultParamsVals = defaultParams(obj);
else
    defaultParamsVals = p.Results.defaultParams;
end

%% Call the tfe version
[paramsFit,fVal,modelResponseStruct] = fitResponse@tfe(obj,thePacket,varargin{:},'defaultParams',defaultParamsVals);

