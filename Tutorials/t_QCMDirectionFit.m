% Test the QCM fitting when the design involves explict directions

% Description:
%   This script synthesizes data for the QCM model and then fits it,
%   to ensure that we can get back what we put in.  Etc. It also
%   demonstrates fitting contrast response functions in separate directions
%   with variously constrained Naka-Rushton functions.
% 
%   Note that this is developed for the two-dimensional (ellipse) version.
% 
%   NOTE: When we get to ellipses, we need to remember to put a constraint into the
%   fitting that keeps the length of the third axis smaller than the second.
%

% History:
%   11/20/18  dhb, mab  Tuned this up; things are making sense.
%   1//24/18  dhb       A bit more tuning.  Fix row/col convention.
%   11/25/18  dhb       Demonstrate code that fits NR functions, and allows
%                       locking of parameters across directions.
%   12/8/18   dhb       Start so can make tfe objects for the various NR
%                       fits

%% Initialize
clear; close all

%% Set up params
%
% Naka-Rushton Params
FIT_NAKARUSHTON = true;
NOOFFSET = false;
theDimension = 2;
Rmax   = 0.9;
sigma  = 0.1;
n      = 2.1;
if (NOOFFSET)
    offset = 0;
else
    offset = -0.1;
end

% Ellipse parameters
minorAxis = 0.3;
rotdeg    = 45;
ellParams = [1 minorAxis rotdeg];

%% Set up QCM matching parameters above
%
% Keep noise very small for testing
QCMObj = tfeQCM('verbosity','none','dimension',theDimension);
paramsQCM = QCMObj.defaultParams;
paramsQCM.Qvec = [minorAxis rotdeg];
paramsQCM.crfAmp = Rmax;
paramsQCM.crfSemi = sigma;
paramsQCM.crfExponent = n;
paramsQCM.noiseSd = 0.01;
paramsQCM.crfOffset = offset;
paramsQCM.expFalloff = 0.3;
paramsQCM.noiseSd = 0.003;
fprintf('\nSimulated QCM parameters:\n');
QCMObj.paramPrint(paramsQCM);

%% Generate stimulus
%
% Random stimuli bounded between 1 and -1
RANDOM_STIMULI = false;
if (RANDOM_STIMULI)
    numStim = 300;
    stimuli = (2*rand(2,numStim) -1);
    
% Stimuli that vary along specified directions
else
    % This matches how we tend to think about our stimuli.  For fitting the
    % Naka_Rushton, it's good to build things up this way, because then the
    % directions that come out of the stimuli match up exactly with the
    % unique directions that we started with.
    %
    % Speicfy stimulus contrasts to use in each direction.
    maxContrast = 0.6;
    contrastsInEachDirection = maxContrast*[0.0625, 0.125, 0.25, 0.5, 1];
    nContrastsPerDirection = length(contrastsInEachDirection);
    
    % Specify the directions, which must each have unit length.
    indDirections = [ [1 0]',  [1 1]',  [0 1]', [1 -1]' ];
    for ii = 1:size(indDirections,2)
        indDirections(:,ii) = indDirections(:,ii)/norm(indDirections(:,ii));
    end
    nUniqueDirections = size(indDirections,2);
    
    % Construct the contrasts crossed with directions set of stimuli
    stimuli = kron(indDirections',contrastsInEachDirection')';
    numStim = size(stimuli,2);
    
    % Get directions/contrasts format from stimuli, for later use
    [stimDirections,stimContrasts] = tfeQCMStimuliToDirectionsContrasts(stimuli);
end

%% Generate response by hand
% This does by hand what our compute response routine does, and is here
% just to expose that calculation and make sure it does what we think it
% should.
%
% Set up scale matrix
S = diag(ellParams(1:2));

% Set up rotation matrix
%
% Whether or not we tag on a transpose here defines 
% the convention for angle. Currently matched to the
% tfeQCM implementation and seems like the right convention
% for the two-dimensional case at least.
V = deg2rotm(ellParams(3))';  

% Get the Q matrix that takes stimulus and get a radius. 
A = S*V';
Q = A'*A;

% Get the radius
radius = diag(sqrt(stimuli'*Q*stimuli))';

% Get the neural response values using the Naka-Rushton function
QCMResponsesByHand = tfeQCMComputeNakaRushton(radius,sigma,n,Rmax,offset);

%% Let's check that the QCM forward model gives the same responses.
%
% Construct the model object 
stimulusStruct.values = stimuli;
stimulusStruct.timebase = 1:numStim;

% Set parameters and simulate responses
QCMResponseStruct = QCMObj.computeResponse(paramsQCM,stimulusStruct,[],'addNoise',false);
if (max(abs(QCMResponsesByHand-QCMResponseStruct.values)/max(QCMResponseStruct.values(:))) > 1e-6)
    error('Hand computation of QCM model does not match tfeQCM forward model');
end

%%  Use the tfeQCM to fit the stim/resp:
%
% Get the tfeQCM object
QCMNoisyResponseStruct = QCMObj.computeResponse(paramsQCM,stimulusStruct,[],'addNoise',true);

% Construct a packet for the QCM to fit.
thePacket.stimulus = stimulusStruct;
thePacket.response = QCMNoisyResponseStruct;
thePacket.kernel = [];
thePacket.metaData = [];

% Fit the packet
if (NOOFFSET)
    defaultParamsInfo.noOffset = true;
else
    defaultParamsInfo.noOffset = false;
end
[fitQCMParams,fVal,fitResponseStructQCM] = QCMObj.fitResponse(thePacket,'defaultParamsInfo',defaultParamsInfo);
fprintf('\nQCM parameters from fit:\n');
QCMObj.paramPrint(fitQCMParams)

%%  Check that the fit recovers the responses we put in to reasonable approximation
% This will break if we simulate too much noise
fitQCMResponseStruct = QCMObj.computeResponse(fitQCMParams,stimulusStruct,[],'addNoise',false);
if (max(abs(fitQCMResponseStruct.values-QCMResponseStruct.values)/max(QCMResponseStruct.values(:))) > 1e-2)
    error('Fit does not do a good job of recovering responses');
end

%% Verify QCMDirection fit 
if (~RANDOM_STIMULI)
    % Create direction QCM object
    QCMDirectionObj = tfeQCMDirection('verbosity','none','dimension',theDimension);

    % Use direction object to compute response. Should get same answer as
    % with regular QCM above.
    directionStimulusStruct = stimulusStruct;
    directionStimulusStruct.values(1:theDimension,:) = stimDirections;
    directionStimulusStruct.values(theDimension+1,:) = stimContrasts;
    QCMDirectionResponseStruct = QCMDirectionObj.computeResponse(paramsQCM,directionStimulusStruct,[],'AddNoise',false);
    if (max(abs(QCMDirectionResponseStruct.values-QCMResponseStruct.values)/max(QCMResponseStruct.values(:))) > 1e-8)
        error('Direction and direct QCM objects do not compute same responses');
    end
    
    % Now fit with direction object and make sure that works the same.
    theDirectionPacket = thePacket;
    theDirectionPacket.stimulus = directionStimulusStruct;
    [fitQCMDirectionParams,fVal,fitQCMDirectionResponseStruct] = QCMDirectionObj.fitResponse(theDirectionPacket,'defaultParamsInfo',defaultParamsInfo);
    fprintf('\nQCM parameters from direction fit:\n');
    QCMDirectionObj.paramPrint(fitQCMDirectionParams)
    if (max(abs(fitQCMDirectionResponseStruct.values-fitQCMResponseStruct.values)/max(fitQCMResponseStruct.values(:))) > 1e-8)
        error('Direction and direct QCM objects do not fit the same');
    end
end

%%  Check that Naka-Rushton funciton inverts
thresholdResponse = fitQCMParams.crfAmp/3;
eqContrast = InvertNakaRushton([fitQCMParams.crfAmp,fitQCMParams.crfSemi,fitQCMParams.crfExponent],thresholdResponse-fitQCMParams.crfOffset);
circlePoints = eqContrast*UnitCircleGenerate(numStim);
[~,Ainv,Q] = EllipsoidMatricesGenerate([1 fitQCMParams.Qvec],'dimension',2);
ellipsePoints = Ainv*circlePoints;
checkThresh = tfeQCMComputeNakaRushton(diag(sqrt(ellipsePoints'*Q*ellipsePoints)),fitQCMParams.crfSemi,fitQCMParams.crfExponent,fitQCMParams.crfAmp,fitQCMParams.crfOffset);
checkThresh1 = ComputeNakaRushton([fitQCMParams.crfAmp,fitQCMParams.crfSemi,fitQCMParams.crfExponent],diag(sqrt(ellipsePoints'*Q*ellipsePoints))) + fitQCMParams.crfOffset;
if (any(abs(checkThresh-thresholdResponse) > 1e-10))
    error('Did not invert QCM model correctly');
end
if (any(abs(checkThresh-checkThresh1) > 1e-10))
    error('Naka-Rushton issue in two ways of computing it');
end

%%  Find contrast in given direction that produces desired response
if (~RANDOM_STIMULI)
    maxResponseFactor = 3;
    whichDirection = 2;
    
    % Pull out this direction and simulated responses
    theDirection = indDirections(:,whichDirection);
    directionResponses = QCMResponsesByHand(nContrastsPerDirection*(whichDirection-1)+1:nContrastsPerDirection*whichDirection);
    maxResponse = max(directionResponses);
    
    % Invert model for chosen direction
    [contrastFromSim,stimulusFromSim] = tfeQCMInvertDirection(paramsQCM,theDirection,maxResponse/maxResponseFactor);
    [contrastFromFit,stimulusFromFit] = tfeQCMInvertDirection(fitQCMParams,theDirection,maxResponse/maxResponseFactor);

    % Plot simulated CRF and inverted points
    figure; hold on
    plot(contrastsInEachDirection,directionResponses,'ro','MarkerFaceColor','r','MarkerSize',12);
    plot(contrastFromSim,maxResponse/maxResponseFactor,'bo','MarkerFaceColor','b','MarkerSize',8);
    plot(contrastFromFit,maxResponse/maxResponseFactor,'gx','MarkerSize',14); 
    xlabel('Contrast'); ylabel('Response');
    
    % Plot an isoresponse contour of the simualted and fit model
    nTheta = 100;
    circleDirections = UnitCircleGenerate(nTheta);
    [contrasts1,stimuli1] = tfeQCMInvertDirection(paramsQCM,circleDirections,paramsQCM.crfAmp/3);
    figure; hold on
    plot(stimuli1(1,:),stimuli1(2,:),'r','LineWidth',3);
    [contrastsFit,stimuliFit] = tfeQCMInvertDirection(fitQCMParams,circleDirections,paramsQCM.crfAmp/3);
    plot(stimuliFit(1,:),stimuliFit(2,:),'b','LineWidth',2);
end

%% Fit Naka-Rushton function to individual directions
if (~RANDOM_STIMULI & FIT_NAKARUSHTON)
    [indDirectionNRParams,indDirectionPredictions,indDirectionResponses,indDirectionDirections,indDirectionContrasts,indDirectionIndices,nIndDirections] = ...
        tfeQCMFitNakaRushtonDirectionsContrasts(QCMResponsesByHand,stimDirections,stimContrasts);
    
    % Check that directions, contrasts and responses came out the way they
    % went in. If this isn't right, then the fit parameters we get back are
    % going to be nonsense.  The biggest worry is that unique() will decide
    % two directions are different because of some numerical tolerance
    % issue.
    inputCounter = 1;
    if (nIndDirections ~= nUniqueDirections)
        error('Did not recover correct number of unique directions from stimulus description');
    end
    for ii = 1:nIndDirections
        for jj = 1:length(indDirectionIndices{ii})
            if (max(abs(stimDirections(:,inputCounter)-indDirectionDirections{ii})) > 1e-7)
                error('Did not properly recover stimulus directions from stimulus description');
            end
            if (QCMResponsesByHand(inputCounter) ~= indDirectionResponses{ii}(jj))
                error('Did not properly recover responses for independent directions');
            end
            if (max(abs(contrastsInEachDirection(jj) - indDirectionContrasts{ii}(jj))) > 1e-7)
                error('Did not properly recover contrasts in each direction for independent directions');
            end
            if (max(abs(indDirectionContrasts{ii}(jj)-stimContrasts(jj))) > 1e-7)
                error('Did not properly recover stimulus contrasts from stimulus description');
            end
            inputCounter = inputCounter+1;
        end
    end
    
    % Fit with things common across directions
    [indDirectionNRParamsCommon] = ...
        tfeQCMFitNakaRushtonDirectionsContrasts(QCMResponsesByHand,stimDirections,stimContrasts,...
        'lockOffsetToZero',false,'commonAmp',true,'commonSemi',false,'commonExp',true,'commonOffset',true);
    
    % Make plot of the individual contrast-response functions and fits
    figure; clf;
    for ii = 1:nIndDirections
        subplot(nIndDirections,1,ii); hold on;
        
        % Dump parameters from independent fits
        fprintf('Parameters for independent direction %d, independent fit\n',ii);
        indDirectionNRParams(ii)
        fprintf('\n');
        
        % Plot simulated data
        plot(indDirectionContrasts{ii},indDirectionResponses{ii},'ro','MarkerFaceColor','r','MarkerSize',12);
        
        % Compute and plot predicted functions
        plotContrasts = linspace(0,max(indDirectionContrasts{ii}),100);
        plotPredictions = ComputeNakaRushton([indDirectionNRParams(ii).crfAmp,indDirectionNRParams(ii).crfSemi,indDirectionNRParams(ii).crfExponent],plotContrasts) + indDirectionNRParams(ii).crfOffset;
        plot(plotContrasts,plotPredictions,'b','LineWidth',4);
        
        % Compute and plot predicted functions, common amplitude
        plotContrasts = linspace(0,max(indDirectionContrasts{ii}),100);
        plotPredictionsCommon = ComputeNakaRushton([indDirectionNRParamsCommon(ii).crfAmp,indDirectionNRParamsCommon(ii).crfSemi,indDirectionNRParamsCommon(ii).crfExponent],plotContrasts) + indDirectionNRParamsCommon(ii).crfOffset;
        plot(plotContrasts,plotPredictionsCommon,'g','LineWidth',2);
    end
    
    % Now try with the tfeNakeRushtonDirection object.  Because we fit the
    % NR to noise free responses that were in turn generated with an NR in
    % each direction, we should get pretty much the same responses in
    % order, as we put in.
    for ii = 1:length(indDirectionDirections)     
        indDirections(:,ii) = indDirectionDirections{ii};
    end
    stimulusStruct.values = [stimDirections ; stimContrasts];
    stimulusStruct.timebase = 1:size(stimulusStruct.values,2);
    NRObj = tfeNakaRushtonDirection(indDirections);
    objResponses = NRObj.computeResponse(indDirectionNRParamsCommon,stimulusStruct,[]);
    if (max(abs(QCMResponsesByHand-objResponses.values)/max(QCMResponsesByHand)) > 1e-6)
        error('tfeNakaRushtonDirection object computeResponse method does not give right answer');
    end
    
end

