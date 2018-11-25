function [paramsFit,fVal,modelResponseStruct] = fitResponse(obj,thePacket,varargin)
% [paramsFit,fVal,modelResponseStruct] = fitResponse(obj,thePacket,varargin)
%
% Fit method for the tfeQCM class.  This overrides the tfe method, for the
% main purpose of allowing us to have multiple starting points chosen
% sensibly for the QCM model.
%
% Inputs:
%   thePacket: a valid packet
%
% Optional key/value pairs
%   See tfe.fitResponse for these.
%
% Outputs:
%   paramsFit: fit parameters
%   fVal: mean value of fit error, mean taken over runs.
%   predictedResponse: big vector containing the fit response

%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
% p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
% p.addRequired('thePacket',@isstruct);
% p.addParameter('defaultParamsInfo',[],@(x)(isempty(x) | isstruct(x)));
% p.addParameter('defaultParams',[],@(x)(isempty(x) | isstruct(x)));
% p.addParameter('searchMethod','fmincon',@ischar);
% p.addParameter('DiffMinChange',[],@isnumeric);
% p.addParameter('fminconAlgorithm','active-set',@ischar);
% p.addParameter('errorType','rmse',@ischar);
% p.parse(thePacket,varargin{:});


% Pass on through to parent class method
if (obj.dimension == 2)
    defaultParams = defaultParams(obj);
    [paramsFit,fVal,modelResponseStruct] = fitResponse@tfe(obj,thePacket,varargin{:});
else
     [paramsFit,fVal,modelResponseStruct] = fitResponse@tfe(obj,thePacket,varargin{:});
end   

end



