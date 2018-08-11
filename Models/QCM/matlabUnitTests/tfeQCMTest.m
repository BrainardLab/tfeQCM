classdef tfeQCMTest < matlab.unittest.TestCase
    % Basic unit tests for QCM model code, using Matlab's unit test
    % framework.  Change to directory that contains this program and 
    % execute "runtests" at the Matlab command line.
    
    properties
        testFigure;
    end
    
    % Initial setup for tests
    methods (TestClassSetup)
        function initTests(obj)
            obj.testFigure = figure; clf; hold on
        end
    end
    
    % Teardown after tests
    methods (TestClassTeardown)
        function closeTests(obj)
            close(obj.testFigure);
        end
    end
    
    methods (Test)
        % Test that instantiating the class doesn't crash
        function tfeQCMConstructorTest(obj)
            tfe = tfeQCM;
        end
        
        % Test that we can set verbosity
        function tfeQCMVerbosityTest(obj)
            tfe = tfeQCM('verbosity','high');
            if (~strcmp(tfe.verbosity,'high'))
                error('Failed to set verbosity to high on create');
            end
        end
        
        % Test that paramsToVec and vecToParams invert
        % 
        % This also tests that defaultParams returns a parameter struct and
        % that print prints it out.
        function paramsToVecTest(obj)
            tfe = tfeQCM;
            params0 = tfe.defaultParams;
            tfe.paramPrint(params0);
            
            x0 = tfe.paramsToVec(params0);
            x1 = x0;
            x1(1) = 2;
            x1(2) = 0.5;
            x1(3) = pi/2;
            x1(7) = 3;
            params1 = tfe.vecToParams(x1);
            x2 = tfe.paramsToVec(params1);
            obj.assertEqual(x1,x2);
        end
        
        % Test that we can simulate a neural response
        function neuralResponseTest(obj)
            % Construct the model object
            tfe = tfeQCM;
            
            % Set parameters
            params0 = tfe.defaultParams;

            % Set the timebase we want to compute on
            deltaT = 1;
            totalTime = 1000;
            stimulusStruct.timebase = 0:deltaT:totalTime;
            
            % Specify the stimulus.
            nTimeSamples = size(stimulusStruct.timebase,2);
            filter = fspecial('gaussian',[1 nTimeSamples],6);
            stimulusStruct.values = rand(3,nTimeSamples);
            for i = 1:3
                stimulusStruct.values(i,:) = ifft(fft(stimulusStruct.values(i,:)) .* fft(filter));
            end
            
            % Generate response and plot. No convolution kernal passed.
            modelResponseStruct = tfe.computeResponse(params0,stimulusStruct,[],'AddNoise',true);
            figure(obj.testFigure);
            tfe.plot(modelResponseStruct,'NewWindow',false);
        end
    end
    
end