function newCriterionResp = scaleToFitIsoContrast(obj,params,isoContrast,varargin)
% Find criterion response such that LCM isocontrast contour best fits a given one.
%
% Syntax:
%    newCriterionResp = scaleIsoContrast(obj,params,varargin)
%
% Description:
%    Find the criterion response so as to change the size of the isoresponse
%    contour to fit a passed isoresponse contour.
%
%    Once you have the return value from this method, set the criterion
%    response to that value. This function does not change the object.
%    So usage would be:
%       newCriterionResp = theObj.scaleToFitIsoContrast(params,isoContrast)
%       saveCriterionResp = theObj.criterionResp;
%       theObj.criterionResp = newCriterionResp;
%       isoContrastFit = theObj.getIsoContrast(params);
%       theObj.criterionResp = saveCriterionResp;
%
% Optional key/value pairs
%   None.

% History
%   03/12/21  dhb  Wrote it.%
%   04/03/21  dhb  Rewrote based on fmincon 

% Parse input.
p = inputParser;
p.addRequired('params',@isstruct);
p.addRequired('isoContrast',@isnumeric);
p.parse(params,isoContrast,varargin{:});

% Save criterion response
saveCriterionResp = obj.criterionResp;

% fmincon options
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','iter','LargeScale','off');
vlbVec = params.crfOffset;
vubVec = params.crfAmp+params.crfOffset;

% Do the fit
newCriterionResp = fmincon(@(x)FitFunction(x,obj,params,isoContrast), ...
    saveCriterionResp,[],[],[],[],vlbVec,vubVec,[],options);

% Put criterionResp back
obj.criterionResp = saveCriterionResp;

end

function f = FitFunction(x,obj,params,isoContrast)

% Set criterion
obj.criterionResp = x;

% Get isocontrast
isoContrastFit = obj.getIsoContrast(params);

% Compare with target
isoDifference = isoContrast-isoContrastFit;

% Return RMSE
f = sqrt(mean(isoDifference.^2));

end
