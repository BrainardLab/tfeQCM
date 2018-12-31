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

% Some custom fitting
if (obj.dimension == 2)
    % Get initial parameters
    if (isempty(p.Results.defaultParams))
        initialParamsVals = obj.defaultParams;
    else
        initialParamsVals = p.Results.defaultParams;
    end
    [paramsFit1,fVal1,modelResponseStruct1] = fitResponse@tfe(obj,thePacket,varargin{:},'defaultParams',initialParamsVals);
    
    % Perturb angle by 90 degrees and fit again
    initialParamsVals.Qvec(2) = initialParamsVals.Qvec(2)-90;
    [paramsFit2,fVal2,modelResponseStruct2] = fitResponse@tfe(obj,thePacket,varargin{:},'defaultParams',initialParamsVals);
    
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
    
    % Put angle into canonical range
    if  (paramsFit.Qvec(2) < -90)
        paramsFit.Qvec(2) = paramsFit.Qvec(2) + 180;
    elseif (paramsFit.Qvec(2) > 90)
        paramsFit.Qvec(2) = paramsFit.Qvec(2) - 180;
    end
else
     [paramsFit,fVal,modelResponseStruct] = fitResponse@tfe(obj,thePacket,varargin{:});
end   

end



