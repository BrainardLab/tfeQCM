function tfeQCMValidateFullAll(varargin)
% tfeQCValidateFullAll(varargin)
%
% Full data check (no figures, no publish) of all validation functions
%
% Optional key/value pairs
%    'asAssertion' - true/false (default false).  Run as an assertion? (for build integration).


%% Parse input and set settable prefs
p = inputParser; p.PartialMatching = false;
p.addParameter('asAssertion',false,@islogical);
p.parse(varargin{:});

%% Close figs
close all;

%% Now check tutorials
tutorialStatus = tfeQCMRunTutorialsAll;
if (p.Results.asAssertion)
    assert(tutorialStatus, 'One or more validations failed.');
end

%% And examples
exampleStatus = tfeQCMRunExamplesAll;
if (p.Results.asAssertion)
    assert(exampleStatus, 'One or more examples failed.');
end

%% Report
if (tutorialStatus & exampleStatus)
    fprintf('\n*** ISET3D validations PASS ***\n');
else
    fprintf('\n*** ISET3D validations FAIL ***\n');
    if (~tutorialStatus)
        fprintf('\tOne or more tutorials failed\n');
    end
    if (~exampleStatus)
        fprintf('\tOne or more examples failed\n');
    end
end
        
        

end