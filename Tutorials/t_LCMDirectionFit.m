% Test the LCM fitting when the design involves explict directions
%
% Description:
%   This script synthesizes data for the LCM model and then fits it,
%   to ensure that we can get back what we put in.  Etc.
% 
%   The parameters used to synthesize the data were just pulled out of a
%   hat, so there is no particular expected shape the resulting specified
%   LCM isoresponse contours in Figure 1 and 2.  The models are then fit to
%   match an ellipse in Figure 3, so there something close to an ellipse is
%   expected.
%
%   There are some local minimum issues with the fitting parameter search
%   in some cases.  The default example parameters are set so these are
%   avoided in both 2019a and 2021a, and the search starts at multiple
%   points as an additional attempt to avoid this issue.  But if you play
%   around with parameters you may find the search gets stuck which will
%   cause some error checks to throw.
%

% History:
%   03/01/21  dhb       Wrote from QCM test version.
%   04/04/21  dhb       Lots of changes.
%   04/10/21  dhb       Illustrate regression scaling of isocontrast
%                       contours.

%% Initialize
clear; close all
rng(0);

%% Test case specific parameter choices
%
% See switch statement below for options and explanations.
whichTest = 'qcmFiveChannel';
switch (whichTest)
    case 'brouwerHeegerBasic'
        % This implements our version of the Brouwer and Heeger model.
        % Six channels with cos^2 sensitivity, linear combination.
        %
        % Channel properties
        nChannels = 6;
        channelExponent = 2;
        summationExponent = 1;
        startCenter = 0;
        channelWeightsPos = [0.5 0.1 0.2];

    case 'kimEtAl'
        % This implements the Kim et al. variant of
        % the Brouwer and Heeger model.
        % Eight channels with cos^4 sensitivity, linear combination.
        %
        % Channel properties
        nChannels = 8;
        channelExponent = 6;
        summationExponent = 1;
        startCenter = 0;
        channelWeightsPos = [1 0.6 0.2 0.4];

    case 'qcmTwoChannel'
        % This implements a QCM with L+M and L-M underlying mechanisms.
        %
        % Channel properties
        nChannels = 4;
        channelExponent = 1;
        summationExponent = 2;
        startCenter = 45;
        channelWeightsPos = 0.8*[1 1/(0.4.^2)];

    case 'qcmFiveChannel'
        % This implements a QCM with five of underlying
        % mechanisms. Note that the simulated data underdetermine
        % this model, so that although the fit to the simulated data
        % reproduce the fit responses, they don't recover the parameters we
        % put in.  An example of how things can be ambiguous in this
        % general space.
        %
        % Channel properties
        nChannels = 10;
        channelExponent = 1;
        summationExponent = 2;
        startCenter = 0;
        channelWeightsPos = [1 0.6 0.2 0.4 0.1];

    otherwise
        error('Unknown test case specified');
end

%% Ellipse to fit
ellipseAngle = 45;          % Angle
ellipseAspectRatio = 0.4;   % Minor axis aspect ratio

%% Set up params
%
% Naka-Rushton Params
NOOFFSET = false;
LOCKEDOFFSET = true;

% Other parameters.  Keep Rmax larger than criterionResp-offset
theDimension = 2;
Rmax   = 2;
sigma  = 0.1;
n      = 2.1;
if (NOOFFSET)
    offset = 0;
else
    offset = -0.1;
end

% Fit error scalar.  Big value seems
% to work better here.
fitErrorScalar = 10000;
maxIter = 1000;

% Follow normalization convention with channel weights.
channelWeightsPos = channelWeightsPos/norm(channelWeightsPos);

% Criterion response for level of isoresponse contours
criterionResponse = 1;

%% Set up LCM matching parameters above
%
% Keep noise very small for testing
if (LOCKEDOFFSET)
    LCMObj = tfeLCMDirection('verbosity','none','dimension',theDimension,'lockedCrfOffset',offset,'criterionResp',criterionResponse, ...
        'nChannels',nChannels,'channelExponent',channelExponent,'summationExponent',summationExponent,'startCenter',startCenter);
else
    LCMObj = tfeLCMDirection('verbosity','none','dimension',theDimension,'criterionResp',criterionResponse, ...
        'nChannels',nChannels,'channelExponent',channelExponent,'summationExponent',summationExponent,'startCenter',startCenter);
end
paramsLCM = LCMObj.defaultParams;
paramsLCM.channelWeightsPos = channelWeightsPos;
paramsLCM.crfAmp = Rmax;
paramsLCM.crfSemi = sigma;
paramsLCM.crfExponent = n;
paramsLCM.noiseSd = 0.01;
paramsLCM.crfOffset = offset;
paramsLCM.expFalloff = 0.3;
paramsLCM.noiseSd = 0.005;
fprintf('\nSimulated LCM parameters:\n');
LCMObj.paramPrint(paramsLCM);

%% Generate stimulus
%
% Stimuli that vary along specified directions

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
nUniqueDirections = 10;
indDirectionAngles = linspace(1,360,nUniqueDirections);
indDirectionDirections(1,:) = cosd(indDirectionAngles);
indDirectionDirections(2,:) = sind(indDirectionAngles);
for ii = 1:size(indDirectionDirections,2)
    indDirectionDirections(:,ii) = indDirectionDirections(:,ii)/norm(indDirectionDirections(:,ii));
end

% Construct the contrasts crossed with directions set of stimuli
stimuli = kron(indDirectionDirections',contrastsInEachDirection')';
numStim = size(stimuli,2);

% Get directions/contrasts format from stimuli, for later use
[stimDirections,stimContrasts] = tfeQCMStimuliToDirectionsContrasts(stimuli);

% Set up stimulus struct
stimulusStruct.values = stimuli;
stimulusStruct.timebase = 1:numStim;

% Put into direction format
directionStimulusStruct = stimulusStruct;
directionStimulusStruct.values(1:theDimension,:) = stimDirections;
directionStimulusStruct.values(theDimension+1,:) = stimContrasts;

%% Generate response and isocontrast contour
LCMResponseStruct = LCMObj.computeResponse(paramsLCM,directionStimulusStruct,[],'addNoise',false);
isoContrastLCM = LCMObj.getIsoContrast(paramsLCM);
angleSupport = LCMObj.angleSupport;

%%  Use the tfeLCMDirection object to fit noisy LCM responses
LCMNoisyResponseStruct = LCMObj.computeResponse(paramsLCM,directionStimulusStruct,[],'addNoise',true);

% Construct a packet for the LCM to fit.
thePacket.stimulus = directionStimulusStruct;
thePacket.response = LCMNoisyResponseStruct;
thePacket.kernel = [];
thePacket.metaData = [];

% Fit the packet
paramsLCMFit = LCMObj.fitResponse(thePacket,'fitErrorScalar',fitErrorScalar,...
    'maxIter',maxIter,'maxFunEval',maxIter);
fprintf('\nLCM parameters from fit to LCM:\n');
LCMObj.paramPrint(paramsLCMFit)

% Check that the fit recovers the responses we put in to reasonable approximation
%
% This will break if we simulate too much noise
fitLCMResponseStruct = LCMObj.computeResponse(paramsLCMFit,directionStimulusStruct,[],'addNoise',false);
if (max(abs(fitLCMResponseStruct.values-LCMResponseStruct.values)/max(LCMResponseStruct.values(:))) > 5e-2)
    error('Fit does not do a good job of recovering responses');
end
isoContrastLCMFit = LCMObj.getIsoContrast(paramsLCMFit);

%% Scale isocontrast contour.
%
% By changing the crition response in the LCM object, we can
% get an iso response contour of a different scale.  Be
% sure not to make the criterion response larger than
% the maximum response in the Naka-Rushton.  Also note
% that becaues the Naka-Rushton is non-linear, the scaling
% of the isoresponse contour will not match the change in
% criterion response.
%
% You probably wouldn't do this for a practical purpose, but
% it's not bad to understand the the ideas.  See below for a 
% regression method for scaling one isocontrast contour to another
% for plotting purposes.
critResponseScaleFactor = 0.5;
saveCriterionResp = LCMObj.criterionResp;
LCMObj.criterionResp = critResponseScaleFactor*LCMObj.criterionResp;
[isoContrastScale] = LCMObj.getIsoContrast(paramsLCM);
LCMObj.criterionResp = saveCriterionResp;

%% Plot isoresponse contours
figure; hold on
plot(isoContrastLCM.*cosd(angleSupport),isoContrastLCM.*sind(angleSupport),'k','LineWidth',5);
plot(isoContrastLCMFit.*cosd(angleSupport),isoContrastLCMFit.*sind(angleSupport),'r','LineWidth',3);
plot(isoContrastScale.*cosd(angleSupport),isoContrastScale.*sind(angleSupport),'b','LineWidth',2);
axis('square');
xlim([-2 2]); ylim([-2 2]);
xlabel('Cone 1 Contrast');
ylabel('Cone 2 Contrast');
title('IsoContrast');
legend({'LCM', 'LCM Fit', sprintf('LCM Scaled by %g in crit resp',critResponseScaleFactor)},'Location','NortheastOutside');

%% Generate an ellipsoidal isoresponse contour
ellipticalIsoContrast = tfeEllipticalIsoContrast(ellipseAngle,ellipseAspectRatio,LCMObj.angleSupport,LCMObj.criterionResp);

%% Scale an isoresponse contour to fit another one.
%
% This method doesn't adjust or affect the model, it just shows how
% you can use regression to adjust size of one isocontrast to fit
% another, for plotting purposes.  Allows us to compare shapes in
% cases where different models normalize their basic isocontrast
% contours differently.
scaleFactor = (isoContrastLCM')\(ellipticalIsoContrast');
scaledToEllipseIsoContrast = scaleFactor*isoContrastLCMFit; 

% Plot
figure; hold on
plot(ellipticalIsoContrast.*cosd(angleSupport),ellipticalIsoContrast.*sind(angleSupport),'k','LineWidth',5);
plot(isoContrastLCM.*cosd(angleSupport),isoContrastLCM.*sind(angleSupport),'r','LineWidth',3);
plot(scaledToEllipseIsoContrast.*cosd(angleSupport),scaledToEllipseIsoContrast.*sind(angleSupport),'b','LineWidth',2);
axis('square');
xlim([-2 2]); ylim([-2 2]);
xlabel('Cone 1 Contrast');
ylabel('Cone 2 Contrast');
title('IsoContrast');
legend({'Ellipse', 'LCM' 'LCM Scaled To Fit'},'Location','NortheastOutside');

%% Now fit the model to the ellipsoidal isoresponse contour
%
% This actually tries to fit the ellipse isoresopnse contour,
% by fitting data that correspond to the iso response contrasts.
%
% Set up structure describing stimuli around the circle.  We generate
% responses for several different contrasts so that we actually constrain
% the model.
fitStimulusStruct.timebase = 1:3*length(LCMObj.angleSupport);
fitStimulusStruct.values(1,:) = [cosd(LCMObj.angleSupport) cosd(LCMObj.angleSupport) cosd(LCMObj.angleSupport)];
fitStimulusStruct.values(2,:) = [sind(LCMObj.angleSupport) sind(LCMObj.angleSupport) sind(LCMObj.angleSupport)];
fitStimulusStruct.values(3,:) = [ellipticalIsoContrast 0.5*ellipticalIsoContrast 2*ellipticalIsoContrast];

% Want to produce criterion response the ellipticalIsoContrast stimuli, and
% scaled versions for the scale versions of that.  We make sure some of the
% responses exceed LCMObj.criterionResp so that the fit Rmax of the
% Naka-Rushton exceeds the criterion response.  That avoids problems when
% we get the isoresponse contour corresponding to the criterion response
% below.
fitResponseStruct.timebase = 1:3*length(LCMObj.angleSupport);
fitResponseStruct.values = [LCMObj.criterionResp*ones(size(LCMObj.angleSupport)) ...
    0.25*LCMObj.criterionResp*ones(size(LCMObj.angleSupport)) ...
    1.5*LCMObj.criterionResp*ones(size(LCMObj.angleSupport))];

% Construct the packet for the fit
thePacket.stimulus = fitStimulusStruct;
thePacket.response = fitResponseStruct;
thePacket.kernel = [];
thePacket.metaData = [];

% Do the fit.
[LCMFitToEllipseParams,fVal,fitToEllpiseStructLCM] = LCMObj.fitResponse(thePacket,'fitErrorScalar',fitErrorScalar,'noNakaRushton',false, ...
        'maxIter',maxIter,'maxFunEval',maxIter); 
    
% Get isoresponse contour corresponding to fit parameters.
[isoContrastFromLCMEllipseFit] = LCMObj.getIsoContrast(LCMFitToEllipseParams);
fprintf('\nLCM parameters from fit to ellipse:\n');
LCMObj.paramPrint(LCMFitToEllipseParams)

% Plot
figure; hold on
plot(ellipticalIsoContrast.*cosd(angleSupport),ellipticalIsoContrast.*sind(angleSupport),'k','LineWidth',5);
plot(isoContrastFromLCMEllipseFit.*cosd(angleSupport),isoContrastFromLCMEllipseFit.*sind(angleSupport),'r','LineWidth',3);
axis('square');
xlim([-2 2]); ylim([-2 2]);
xlabel('Cone 1 Contrast');
ylabel('Cone 2 Contrast');
title('IsoContrast');
legend({'Ellipse', 'LCM Fit to Ellipse'},'Location','NortheastOutside');


