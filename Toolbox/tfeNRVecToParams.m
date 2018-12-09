function [NRParams,nParams] = tfeNRVecToParams(paramsvec,nIndDirections)
% Convert the NR parameter vector to struct array
%
% Syntax:
%    [NRParams,nParams] = tfeNRVecToParams(paramsvec,nIndDirections)
%
% See also: tfeNRParamsToVec, InitializeNRParams
%
    
    % Number of parameters.  You just have to know this.  Coordinate
    % any changes with the tfeNRParamsToVec routine.
    nParams = 6;
    
    % Unpack vector into struct array.
    for ii = 1:nIndDirections        
        NRParams(ii).crfAmp =  paramsvec((ii-1)*nParams+1);
        NRParams(ii).crfSemi = paramsvec((ii-1)*nParams+2);
        NRParams(ii).crfExponent = paramsvec((ii-1)*nParams+3);
        NRParams(ii).crfOffset = paramsvec((ii-1)*nParams+4);
        NRParams(ii).expFalloff = paramsvec((ii-1)*nParams+5);
        NRParams(ii).noiseSd = paramsvec((ii-1)*nParams+6);
    end
end