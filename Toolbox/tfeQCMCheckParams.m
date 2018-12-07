function [status] = tfeQCMCheckParams(params)
% Check that tfeQCM parameters structure is legitimate
%
% Synopsis
%    [status] = tfeQCMCheckParams(params)
%
% Description:
%    Take stimuli and compute QCM responses.
%
% Inputs:
%      params        - QCM model parameter struct
% 
% Outputs:
%      status        - True if parameters OK, false otherwise.

% History:
%   11/24/19  dhb    Wrote it for modularity.

%% Infer dimension
dimension = 0;
if (length(params.Qvec) == 2)
    dimension = 2;
elseif (length(params.Qvec) == 5)
    dimension = 3;
end

%% Check parameterization
status = true;
tolerance = 1e-3;
if (dimension == 2)
    if (params.Qvec(1) > 10+tolerance)
        status = false;
    end
elseif (dimension == 3)
    if (params.Qvec(1) > 10+tolerance)
        status = false;
    end
    if (params.Qvec(2) > params.Qvec(1)+tolerance)
        status = false;
    end
else
    error('Can only handle dimension of 2 or 3 for stimuli');
end

