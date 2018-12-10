function paramsvec = tfeNRParamsToVec(NRParams)
% Convert the NR parameter structure array to one long vector
%
% Syntax:
%    paramsvec = tfeNRParamsToVec(NRParams)
%
% Description:
%    Convert vector form of Naka-Rushton parameters for multiple color
%    directions into struct array form.
%
% Inputs:
%    NRParams        - Struct array of parameters.  See
%                      tfeNRInitializeParams for the fields.
%
% Outputs:
%    paramsvec       - Vector form of parameters.
%
% Optional key/value pairs:
%    None.
%
% See also: tfeNRVecToParams, tfeInitializeNRParams
%


% History:
%   12/10/18  dhb  Header comments.

    nDirections = length(NRParams);
    paramsvec = [];
    for ii = 1:nDirections
        paramsvec = [paramsvec ; [NRParams(ii).crfAmp NRParams(ii).crfSemi NRParams(ii).crfExponent NRParams(ii).crfOffset NRParams(ii).expFalloff NRParams(ii).noiseSd]'];
    end
end