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

%% Set up params
%
% Naka-Rushton Params
NOOFFSET = false;
LOCKEDOFFSET = true;
theDimension = 2;
Rmax   = 0.9;
sigma  = 0.1;
n      = 2.1;
if (NOOFFSET)
    offset = 0;
else
    offset = -0.1;
end

% Weight parameters
channelWeightsPos = [1 0.6 0.2];

% Ellipse parameters.  Parameters here picked by hand so that
% the scaled LCM isoresponse contour for the weights above
% comes out pretty close.
ellipseAngle = 105;         % Angle
ellipseAspectRatio = 0.4;   % Minor axis aspect ratio

% Criterion response for level of isoresponse contours
criterionResponse = 1;

%% Set up LCM matching parameters above
%
% Keep noise very small for testing
if (LOCKEDOFFSET)
    LCMObj = tfeLCMDirection('verbosity','none','dimension',theDimension,'lockedCrfOffset',offset,'criterionResp',criterionResponse);
else
    LCMObj = tfeLCMDirection('verbosity','none','dimension',theDimension,'criterionResp',criterionResponse);
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
indDirectionDirections = [ [1 0]',  [1 1]',  [0 1]', [1 -1]' ];
for ii = 1:size(indDirectionDirections,2)
    indDirectionDirections(:,ii) = indDirectionDirections(:,ii)/norm(indDirectionDirections(:,ii));
end
nUniqueDirections = size(indDirectionDirections,2);

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

%% Generate response
LCMResponseStruct = LCMObj.computeResponse(paramsLCM,directionStimulusStruct,[],'AddNoise',false);

%%  Use the tfeLCMDirection object to fit the stim/resp:
LCMNoisyResponseStruct = LCMObj.computeResponse(paramsLCM,directionStimulusStruct,[],'addNoise',true);

% Construct a packet for the LCM to fit.
thePacket.stimulus = directionStimulusStruct;
thePacket.response = LCMNoisyResponseStruct;
thePacket.kernel = [];
thePacket.metaData = [];

% Fit the packet
[fitLCMParams,fVal,fitResponseStructLCM] = LCMObj.fitResponse(thePacket);
fprintf('\nLCM parameters from fit:\n');
LCMObj.paramPrint(fitLCMParams)

%% Check that the fit recovers the responses we put in to reasonable approximation
%
% This will break if we simulate too much noise
fitLCMResponseStruct = LCMObj.computeResponse(fitLCMParams,directionStimulusStruct,[],'addNoise',false);
if (max(abs(fitLCMResponseStruct.values-LCMResponseStruct.values)/max(LCMResponseStruct.values(:))) > 5e-2)
    error('Fit does not do a good job of recovering responses');
end

%% Get and plot isoresponse contour
figure; hold on
[isoContrastFit,unitContrastResponseFit,angleSupport] = LCMObj.getIsoContrast(fitLCMParams);
plot(isoContrastFit.*cosd(angleSupport),isoContrastFit.*sind(angleSupport),'r','LineWidth',2);
axis('square');
xlim([-2 2]); ylim([-2 2]);
xlabel('Cone 1 Contrast');
ylabel('Cone 2 Contrast');
title('IsoContrast');

%% Scale isocontrast contour and replot
fitLCMParamsScale = LCMObj.scaleIsoContrast(fitLCMParams,0.5);
[isoContrastScale] = LCMObj.getIsoContrast(fitLCMParamsScale);
plot(isoContrastScale.*cosd(angleSupport),isoContrastScale.*sind(angleSupport),'b','LineWidth',2);

%% Compute responses with original and scaled parameters
fitResponseStruct = LCMObj.computeResponse(fitLCMParams,directionStimulusStruct,[]);
fitResponseStructScale = LCMObj.computeResponse(fitLCMParamsScale,directionStimulusStruct,[]);
if (max(abs(fitResponseStruct.values-fitResponseStructScale.values)./fitResponseStruct.values) > 1e-9)
    error('Scaling does not preserve response');
end

%% Generate an ellipsoidal isoresponse contour
ellipticalIsoContrast = tfeEllipticalIsoContrast(ellipseAngle,ellipseAspectRatio,LCMObj.angleSupport,LCMObj.criterionResp);

%% Get parameters with model parameters scaled to produce best fit to elliptical contour
fitLCMParamsScaledTOEllipse = LCMObj.scaleToFitIsoContrast(fitLCMParams,ellipticalIsoContrast);
scaledToEllipseIsoContrast = LCMObj.getIsoContrast(fitLCMParamsScaledTOEllipse);

% Plot
figure; hold on
plot(ellipticalIsoContrast.*cosd(angleSupport),ellipticalIsoContrast.*sind(angleSupport),'k','LineWidth',2);
plot(isoContrastFit.*cosd(angleSupport),isoContrastFit.*sind(angleSupport),'r','LineWidth',2);
plot(scaledToEllipseIsoContrast.*cosd(angleSupport),scaledToEllipseIsoContrast.*sind(angleSupport),'b','LineWidth',2);
axis('square');
xlim([-2 2]); ylim([-2 2]);
xlabel('Cone 1 Contrast');
ylabel('Cone 2 Contrast');
title('IsoContrast');


