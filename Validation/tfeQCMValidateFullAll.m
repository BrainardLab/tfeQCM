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
        

end