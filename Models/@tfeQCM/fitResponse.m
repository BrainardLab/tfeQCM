function [paramsFit,fVal,modelResponseStruct] = fitResponse(obj,thePacket,varargin)
% [paramsFit,fVal,modelResponseStruct] = fitResponse(obj,thePacket,varargin)
%
% Fit method for the tfeQCM class.  This overrides the tfe method, for the
% main purpose of allowing us to have multiple starting points chosen
% sensibly for the QCM model.
%
% Also eventually will allow us to put on some parameter constraints that
% are not generic.
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
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('thePacket',@isstruct);
p.addParameter('defaultParamsInfo',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('defaultParams',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('searchMethod','fmincon',@ischar);
p.addParameter('DiffMinChange',[],@isnumeric);
p.addParameter('fminconAlgorithm','active-set',@ischar);
p.addParameter('errorType','rmse',@ischar);
p.parse(thePacket,varargin{:});

% Some custom fitting
if (obj.dimension == 2)
    % Fit with standard deafult parameters
    if (isempty(p.Results.defaultParams))
        defaultParamsVals = defaultParams(obj);
    else
        defaultParamsVals = p.Results.defaultParams;
    end
    [paramsFit1,fVal1,modelResponseStruct1] = fitResponse@tfe(obj,thePacket,varargin{:},'defaultParams',defaultParamsVals);
    
    % Perturb angle by 90 degrees and fit again
    defaultParamsVals.Qvec(2) = defaultParamsVals.Qvec(2)-90;
    [paramsFit2,fVal2,modelResponseStruct2] = fitResponse@tfe(obj,thePacket,varargin{:},'defaultParams',defaultParamsVals);
    
    % Pick the winner
    if (fVal1 <= fVal2)
        paramsFit = paramsFit1;
        fVal = fVal1;
        modelResponseStruct = modelResponseStruct1;
    else
        paramsFit = paramsFit2;
        fVal = fVal2;
        modelResponseStruct = modelResponseStruct2;
    end
else
     [paramsFit,fVal,modelResponseStruct] = fitResponse@tfe(obj,thePacket,varargin{:});
end   

end



