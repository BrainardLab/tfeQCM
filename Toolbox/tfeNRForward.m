function responses = tfeNRForward(NRParams,contrasts)
% Compute Naka-Rushton over separate directions, from passed contrasts.
%
% Syntax:
%    responses = tfeNRForward(NRParams,contrasts)   
%
% Description:
%    Take cell array of contrasts in each of a series of directions,
%    corresponding parameters in a struct array, and compute cell array of 
%    responses.
%
% Inputs:
%    NRParams       - Struct array with each entry containing Naka-Rushton
%                     parameteres for one direction.  See
%                     tfeNakaRushton.defaultParams for form of parameter
%                     structure.
%    contrasts      - Cell array of contrasts. Contrasts are in row
%                     vectors of each cell entry.
%
% Outputs:
%    responses      - Cell array of responses. Responses are in row vectors
%                     for each cell entry.
%
% See also:
%

% History:
%   12/10/18  dhb  Wrote it.

% Loop over directions
nDirections = length(NRParams);
for ii = 1:nDirections
    responses{ii} = tfeQCMComputeNakaRushton(contrasts{ii},NRParams(ii).crfSemi,NRParams(ii).crfExponent,NRParams(ii).crfAmp,NRParams(ii).crfOffset);
end