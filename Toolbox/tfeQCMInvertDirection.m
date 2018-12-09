function [contrasts,stimuli] = tfeQCMInvertDirection(params,directions,response)
% Invert QCM model for a particular direction
%
% Synopsis
%    [contrasts,stimuli] = tfeQCMInvertDirection(params,directions,response)
%
% Description:
%    Take a unit vector specifying a direction in stimulus space, the
%    parameters of the QCM model, and a desirec response, and find the
%    contrast that produces that response.  Also returns the corresponding
%    stimulus.
%
% Inputs:
%      params        - QCM model parameter struct
%      directions    - Directions, each column should have unit length
%      response      - Desired response
% 
% Outputs:
%      contrasts     - Contrasts that produces desired response, in a row
%                      vector. NaN if passed response not attainable.
%      stimuli       - Stimuli corresponding to contrasts, given by
%                      contrasts(ii)*directions(:,ii) for each column direction of
%                      passed matrix directions.

% History:
%   11/20/18  dhb, mab  Wrote it.
%   11/24/18  dhb       Work on an input matrix, change input row/col
%                       convention to match QCM. Works with multiple
%                       inputs.
%   11/24/18  dhb       Better name.

%% Paramters check
if (~tfeQCMCheckParams(params))
    error('Bad values in QCM parameters structure');
end

%% Handle the pesky offset and check that there is an in-range contrast.
% Return NaN if not.
offsetResponse = response-params.crfOffset;
if (offsetResponse > params.crfAmp | offsetResponse < 0)
    contrasts = NaN;
    stimuli = NaN*ones(size(directions));
    return;
end

%% Invert Naka rushton to find desired output of quadratic computation.
theDimension = size(directions,1);
desiredEqContrast = InverttfeQCMComputeNakaRushton([params.crfAmp,params.crfSemi,params.crfExponent],offsetResponse);

%% Loop over directions and invert each one
%
% Find what comes out of quadratic for the passed direction.
nDirections = size(directions,2);
contrasts = zeros(1,nDirections);
stimuli = zeros(size(directions));
for ii = 1:nDirections
    direction = directions(:,ii);
    if (abs(norm(direction)-1) > 1e-5)
        error('A passed direciton does not have unit vector length');
    end
    
    % Do the work of inverting the QCM model for one direction
    [~,Ainv,Q] = EllipsoidMatricesGenerate([1 params.Qvec],'dimension',theDimension);
    directionEqContrast = diag(sqrt(direction'*Q*direction));
    contrasts(ii) = desiredEqContrast/directionEqContrast;
    stimuli(:,ii) = contrasts(ii)*direction;
end

