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
% See also: tfeNakaRushtonDirection, tfeNRInvert
%

% History:
%   12/10/18  dhb  Wrote it.

% Examples:
%{
    NRParams(1).crfAmp = 1.5;
    NRParams(1).crfExponent = 2.2;
    NRParams(1).crfSemi = 0.1;
    NRParams(1).crfOffset = -0.1;
    NRParams(2).crfAmp = 0.5;
    NRParams(2).crfExponent = 1.6;
    NRParams(2).crfSemi = 0.2;
    NRParams(2).crfOffset = 0.2;
    contrasts{1} = linspace(0,0.5,10);
    contrasts{2} = linspace(0,1,100);
    responses = tfeNRForward(NRParams,contrasts);
    figure; clf; hold on
    plot(contrasts{1},responses{1},'r','LineWidth',3);
    plot(contrasts{2},responses{2},'g','LineWidth',3);
    xlabel('Contrast');
    ylabel('Response');

    contrastsChk = tfeNRInvert(NRParams,responses);
    plot(contrastsChk{1},responses{1},'k','LineWidth',2);
    plot(contrastsChk{2},responses{2},'k','LineWidth',2);

    if (max(abs(contrasts{1}-contrastsChk{1})) > 1e-8)
       error('NR function does not invert right');
    end
    if (max(abs(contrasts{2}-contrastsChk{2})) > 1e-8)
       error('NR function does not invert right');
    end
%}

% Loop over directions
nDirections = length(NRParams);
for ii = 1:nDirections
    responses{ii} = tfeQCMComputeNakaRushton(contrasts{ii},NRParams(ii).crfSemi,NRParams(ii).crfExponent,NRParams(ii).crfAmp,NRParams(ii).crfOffset);
end