function status = tfeQCMRunExamplesAll
%% Run all the examples in the tfeQCM tree
%
% Syntax:
%     tfeQCMRunExamplesAll
%
% Description:
%     Run all the examples in the tfeQCM tree,
%     excepthose that contain a line of the form
%     "% ETTBSkip"
%
% Inputs:
%    None.
%
% Outputs:
%    status    - 1 if all examples run OK, 0 otherwise.
%
% Optional key/value pairs:
%    None.
%
% See also:
%   ieValidateFullAll, ieRunTutorialsAll

% History:
%   01/17/18  dhb  Wrote it.

[~, functionStatus] = ExecuteExamplesInDirectory(tbLocateToolbox('tfeQCM'),'verbose',false);
status = all(functionStatus ~= -1);