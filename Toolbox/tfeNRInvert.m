function contrasts = tfeNRInvert(NRParams,responses)
% Invert Naka-Rushton over separate directions, from passed responses.
%
% Syntax:
%    contrasts = tfeNRInvert(NRParams,responses)   
%
% Description:
%    Take cell array of responses in each of a series of directions,
%    corresponding parameters in a struct array, and compute cell array of 
%    contrasts.
%
% Inputs:
%    NRParams       - Struct array with each entry containing Naka-Rushton
%                     parameteres for one direction.  See
%                     tfeNakaRushton.defaultParams for form of parameter
%                     structure.
%    responses      - Cell array of responses. Responses are in row vectors
%                     for each cell entry.
%
% Outputs:
%
%    contrasts      - Cell array of contrasts. Contrasts are in row
%                     vectors of each cell entry.
%
% See also: tfeNakaRushtonDirection, tfeNRForward
%

% History:
%   12/10/18  dhb  Wrote it.

% Loop over directions
nDirections = length(NRParams);
for ii = 1:nDirections
    contrasts{ii} = InvertNakaRushton([NRParams(ii).crfAmp,NRParams(ii).crfSemi,NRParams(ii).crfExponent],responses{ii}-NRParams(ii).crfOffset);
end