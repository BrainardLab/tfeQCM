% Test the QCM fitting on an example test packet
%
% Description:
%   This script tests the QCM and NR models on a test packet that
%   contains CRFs for various directions in color space.,
%   to ensure that we can get back what we put in.  Etc. It also
%   demonstrates fitting contrast response functions in separate directions
%   with variously constrained Naka-Rushton functions.
%

% History:
%   12/30/18  dhb  Wrote it.

%% Initialize
clear; close all
rng(0);

%% Load the packet
load theTestPacket

%% Get directions/contrasts
[stimulusDirections,stimulusContrasts] = tfeQCMStimuliToDirectionsContrasts(thePacket.stimulus.values);
[uniqueDirections,directionIndices] = tfeQCMParseDirections(stimulusDirections);
nUniqueDirections = size(uniqueDirections,2);
theDimension = size(uniqueDirections,1);
theDirectionPacket = thePacket;
theDirectionPacket.stimulus.values = [stimulusDirections ; stimulusContrasts];

%% Plot the measured contrast response functions
%
% This is hard coded to the structure of our test packet. May need
% a little tweaking if a different packet is used.
subCols = 4;
subRows = 2;
directionIndex = 1;
minPlotContrast = 0;
maxPlotContrast = 0.6;
minPlotResponse = -0.5;
maxPlotResponse = 2;
crfFig = figure; clf; set(gcf,'Position',[50 550 1200 720]);
for ii = 1:subRows
    for jj = 1:subCols
        subplot(subRows,subCols,directionIndex); hold on;
        directionContrasts = stimulusContrasts(directionIndices{directionIndex});
        directionResponses = thePacket.response.values(directionIndices{directionIndex});
        plot(directionContrasts,directionResponses,'ro','MarkerFaceColor','r','MarkerSize',12);
        xlim([minPlotContrast maxPlotContrast]); ylim([minPlotResponse maxPlotResponse]);
        xlabel('Contrast'); ylabel('Response');
        title(sprintf('Direction: (%0.3f %0.3f)',uniqueDirections(1,directionIndex),uniqueDirections(2,directionIndex)));
        directionIndex = directionIndex+1;
    end
end

%% Set up plot contrasts
nPlotContrasts = 100;
plotContrasts = linspace(minPlotContrast,maxPlotContrast,nPlotContrasts);

%% Fit each CRF with its own Naka-Rushton but with common offset
%
% These go nicely through the data but are crazy ass otherwise.
indNRDirectionObj = tfeNakaRushtonDirection(uniqueDirections, ...
        'lockOffsetToZero',false,'commonAmp',false,'commonSemi',false,'commonExp',false,'commonOffset',true);
[indNRParams,~,indNRResponses] = indNRDirectionObj.fitResponse(theDirectionPacket);  
fprintf('*** NR independent fit params\n');
indNRDirectionObj.paramPrint(indNRParams);
directionIndex = 1;
figure(crfFig);
for ii = 1:subRows
    for jj = 1:subCols
        % Getting the object to compute for plots is hard, because of the
        % need to dummy up the right sort of packet with all the directions
        % in it.  This is easier.
        plotResponses = tfeNRForward(indNRParams(directionIndex),{plotContrasts});    
        subplot(subRows,subCols,directionIndex); hold on;
        plot(plotContrasts,plotResponses{1},'r','LineWidth',4);
        directionIndex = directionIndex+1;
    end
end    
   
%% Fit each CRF with common amplitude/offest Naka-Rushton but with common offset
commonAmpExpNRObj = tfeNakaRushtonDirection(uniqueDirections, ...
        'lockOffsetToZero',false,'commonAmp',true,'commonSemi',false,'commonExp',false,'commonOffset',true);
[commonAmpExpNRParams,~,commonAmpExpNRResponses] = commonAmpExpNRObj.fitResponse(theDirectionPacket);  
fprintf('*** NR common amp fit params\n');
commonAmpExpNRObj.paramPrint(commonAmpExpNRParams);
directionIndex = 1;
figure(crfFig);
for ii = 1:subRows
    for jj = 1:subCols
        plotResponses = tfeNRForward(commonAmpExpNRParams(directionIndex),{plotContrasts});
        subplot(subRows,subCols,directionIndex); hold on;
        plot(plotContrasts,plotResponses{1},'g','LineWidth',4);
        directionIndex = directionIndex+1;
    end
end 

%% Fit each CRF with common amplitude/offest/exp Naka-Rushton but with common offset
commonAmpExpNRObj = tfeNakaRushtonDirection(uniqueDirections, ...
        'lockOffsetToZero',false,'commonAmp',true,'commonSemi',false,'commonExp',true,'commonOffset',true);
[commonAmpExpNRParams,~,commonAmpExpNRResponses] = commonAmpExpNRObj.fitResponse(theDirectionPacket);  
fprintf('*** NR common amp/exp fit params\n');
commonAmpExpNRObj.paramPrint(commonAmpExpNRParams);
directionIndex = 1;
figure(crfFig);
for ii = 1:subRows
    for jj = 1:subCols
        plotResponses = tfeNRForward(commonAmpExpNRParams(directionIndex),{plotContrasts});
        subplot(subRows,subCols,directionIndex); hold on;
        plot(plotContrasts,plotResponses{1},'b','LineWidth',4);
        directionIndex = directionIndex+1;
    end
end 

%% Fit with the QCM direction model
%
% This starts with the amp/exp/offset NR parameters found from the common fit
% above, but that doesn't actually seem to be necessary.
QCMObj = tfeQCMDirection('verbosity','none','dimension',theDimension);
defaultParams = QCMObj.defaultParams;
defaultParams.crfAmp = commonAmpExpNRParams(1).crfAmp;
defaultParams.crfExponent = commonAmpExpNRParams(1).crfExponent;
defaultParams.crfOffset = commonAmpExpNRParams(1).crfOffset;
defaultParamsInfo = [];
[QCMParams,~,QCMResponses] = QCMObj.fitResponse(theDirectionPacket,'defaultParams',defaultParams,'defaultParamsInfo',defaultParamsInfo);
fprintf('\nQCM parameters from direction fit:\n');
QCMObj.paramPrint(QCMParams);
directionIndex = 1;
figure(crfFig);
for ii = 1:subRows
    for jj = 1:subCols
        plotDirection = uniqueDirections(:,directionIndex);
        plotStimuli = tfeQCMDirectionsContrastsToStimuli(plotDirection(:,ones(1,nPlotContrasts)),plotContrasts);
        plotResponses = tfeQCMForward(QCMParams,plotStimuli);
        subplot(subRows,subCols,directionIndex); hold on;
        plot(plotContrasts,plotResponses,'k','LineWidth',4);
        directionIndex = directionIndex+1;
    end
end 

%% Make a plot of QCM fit in contrast form
%
% Because of the offset in the data, a criterion response
% of 0 is reasonable
criterionResponse = 0;
ellPlot = figure; clf; hold on
nTheta = 500;
circleDirections = UnitCircleGenerate(nTheta);
[contrasts1,stimuli1] = tfeQCMInvertDirection(QCMParams,circleDirections,criterionResponse);
figure; hold on
plot(stimuli1(1,:),stimuli1(2,:),'k','LineWidth',3);
xlim([-0.3 0.3]); ylim([-0.3 0.3]);
xlabel('L contrast');
ylabel('M contrast');
axis('square');
for ii = 1:nUniqueDirections
    criterionContrast = tfeNRInvert(indNRParams(ii),{criterionResponse});
    % InvertNakaRushton([indNRParams(ii).crfAmp,indNRParams(ii).crfSemi,indNRParams(ii).crfExponent],{criterionResponse-indNRParams(ii).crfOffset);
    plotDirection = uniqueDirections(:,ii);
    plotStimulus = tfeQCMDirectionsContrastsToStimuli(plotDirection,criterionContrast{1});
    plot(plotStimulus(1),plotStimulus(2),'ro','MarkerFaceColor','r','MarkerSize',12);
    
    criterionContrast = tfeNRInvert(commonAmpExpNRParams(ii),{criterionResponse});
    %InvertNakaRushton([commonAmpExpNRParams(ii).crfAmp,commonAmpExpNRParams(ii).crfSemi,commonAmpExpNRParams(ii).crfExponent],criterionResponse-commonAmpExpNRParams(ii).crfOffset);
    plotDirection = uniqueDirections(:,ii);
    plotStimulus = tfeQCMDirectionsContrastsToStimuli(plotDirection,criterionContrast{1});
    plot(plotStimulus(1),plotStimulus(2),'bo','MarkerFaceColor','b','MarkerSize',12);
end
