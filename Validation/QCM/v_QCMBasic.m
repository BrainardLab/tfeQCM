function varargout = v_QCMBasic(varargin)
% varargout = v_QCMBasic(varargin)
%
% Works by running t_QCMBasic with various arguments and comparing
% results with those stored.
%
% Validate applyKernel method of tfe parent class.

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)
    
    %% Freeze RNG so validations work
    rng(1);
    
    %% Basic validation
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_QCMBasic *****');
    validationData1 = t_QCMBasic('generatePlots',runTimeParams.generatePlots);
    UnitTest.validationData('validationData1',validationData1);
    
end



