% Test the LCM fitting when the design involves explict directions
%
% Description:
%   This script synthesizes data for the LCM model and then fits it,
%   to ensure that we can get back what we put in.  Etc.
% 
%   The parameters used to synthesize the data were just pulled out of a
%   hat, so there is no particular expected shape the resulting isoresponse
%   contour.

% History:
%   03/01/21  dhb       Wrote from QCM test version.

%% Initialize
clear; close all
rng(0);

%% Test case specific parameter choices
whichTest = 'qcmTwoChannel';
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
        channelWeightsPos = [1 0.6 0.2];
        
        % Ellipse to fit
        ellipseAngle = 105;         % Angle
        ellipseAspectRatio = 0.4;   % Minor axis aspect ratio

    case 'kimEtAl'
        % This implements the Kim et al. variant of
        % the Brouwer and Heeger model.
        % Eight channels with cos^4sensitivity, linear combination.
        %
        % Channel properties
        nChannels = 8;
        channelExponent = 4;
        summationExponent = 1;
        startCenter = 0;
        channelWeightsPos = [1 0.6 0.2 0.4];
        
        % Ellipse to fit
        ellipseAngle = 105;         % Angle
        ellipseAspectRatio = 0.4;   % Minor axis aspect ratio

    case 'qcmTwoChannel'
        % This implements a QCM with L+M and L-M underlying mechanisms.
        %
        % Channel properties
        nChannels = 4;
        channelExponent = 1;
        summationExponent = 2;
        startCenter = 45;
        channelWeightsPos = 0.8*[1 1/(0.4.^2)];
        
        % Ellipse to fit
        ellipseAngle = 45;          % Angle
        ellipseAspectRatio = 0.4;   % Minor axis aspect ratio
    case 'qcmFiveChannel'
    otherwise
        error('Unknown test case specified');
end

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
paramsLCMFit = LCMObj.fitResponse(thePacket,'fitErrorScalar',fitErrorScalar);
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

%% Scale isocontrast contour by a specified amount.
%
% The amount is the change in criterion response, so with
% the Naka-Rushton it will not in general be the factor by
% which the contour itself grows or shrinks.
%
% You probably wouldn't do this, but understanding the idea
% is the basis of the code that scales an LCM isocontrast
% contour to fit another one for plotting.
contourScaleFactor = 0.25;
saveCriterionResp = LCMObj.criterionResp;
LCMObj.criterionResp = 0.5*LCMObj.criterionResp;
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
legend({'LCM', 'LCM Fit', sprintf('LCM Scaled by %g in crit resp',contourScaleFactor)},'Location','NortheastOutside');

%% Compute responses with original and scaled parameters
% fitResponseStruct = LCMObj.computeResponse(paramsLCM,directionStimulusStruct,[]);
% fitResponseStructScale = LCMObj.computeResponse(paramsLCMScale,directionStimulusStruct,[]);
% if (max(abs(fitResponseStruct.values-fitResponseStructScale.values)./fitResponseStruct.values) > 1e-9)
%     error('Scaling does not preserve response');
% end

%% Generate an ellipsoidal isoresponse contour
ellipticalIsoContrast = tfeEllipticalIsoContrast(ellipseAngle,ellipseAspectRatio,LCMObj.angleSupport,LCMObj.criterionResp);

%% Get parameters with model parameters scaled to produce best fit to elliptical contour
%
% This shows how to scale an isoresponse contour, but doesn't try to adjust
% relative channel weights to best fit the ellipse.
newCriterionResp = LCMObj.scaleToFitIsoContrast(paramsLCM,ellipticalIsoContrast);
saveCriterionResp = LCMObj.criterionResp;
LCMObj.criterionResp = newCriterionResp;
scaledToEllipseIsoContrast = LCMObj.getIsoContrast(paramsLCM);
LCMObj.criterionResp = saveCriterionResp;

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
% Note use of 'noNakeRushton' flag so that we just work with the
% linear response, since what we care about here is the shape. 
%
% Set up structure describing stimuli around the circle
fitStimulusStruct.timebase = 1:length(LCMObj.angleSupport);
fitStimulusStruct.values(1,:) = cosd(LCMObj.angleSupport);
fitStimulusStruct.values(2,:) = sind(LCMObj.angleSupport);
fitStimulusStruct.values(3,:) = ellipticalIsoContrast;

% Want to produce criterion response for those stimuli
fitResponseStruct.timebase = fitStimulusStruct.timebase;
fitResponseStruct.values = LCMObj.criterionResp*ones(size(fitResponseStruct.timebase));

% Construct the packet for the fit
thePacket.stimulus = fitStimulusStruct;
thePacket.response = fitResponseStruct;
thePacket.kernel = [];
thePacket.metaData = [];

% Do the fit.
[fitToEllipseParams,fVal,fitToEllpiseStructLCM] = LCMObj.fitResponse(thePacket,'fitErrorScalar',fitErrorScalar,'noNakaRushton',false);
[isoContrastFromEllipseFit] = LCMObj.getIsoContrast(fitToEllipseParams);
paramsLCMFit = LCMObj.fitResponse(thePacket,'fitErrorScalar',fitErrorScalar);
fprintf('\nLCM parameters from fit to ellipse:\n');
LCMObj.paramPrint(paramsLCMFit)

% Plot
figure; hold on
plot(ellipticalIsoContrast.*cosd(angleSupport),ellipticalIsoContrast.*sind(angleSupport),'k','LineWidth',5);
plot(isoContrastFromEllipseFit.*cosd(angleSupport),isoContrastFromEllipseFit.*sind(angleSupport),'r','LineWidth',3);
axis('square');
xlim([-2 2]); ylim([-2 2]);
xlabel('Cone 1 Contrast');
ylabel('Cone 2 Contrast');
title('IsoContrast');
legend({'Ellipse', 'LCM Fit to Ellipse'},'Location','NortheastOutside');


