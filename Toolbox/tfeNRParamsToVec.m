function paramsvec = tfeNRParamsToVec(NRParams)
% Convert the NR parameter structure array to one long vector
%
% See also: tfeNRVecToParams, InitializeNRParams
%

    nIndDirections = length(NRParams);
    paramsvec = [];
    for ii = 1:nIndDirections
        paramsvec = [paramsvec ; [NRParams(ii).crfAmp NRParams(ii).crfSemi NRParams(ii).crfExponent NRParams(ii).crfOffset NRParams(ii).expFalloff NRParams(ii).noiseSd]'];
    end
end