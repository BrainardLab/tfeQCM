% Show how we can use functions available with the QCM to fit contrast-resp functions.
%
% Description:

% History:
%  08/28/18  dhb  Wrote it.

%% Clear and close
clear; close all;

%% Set parameters
theDimension = 2;
generatePlots = true;

%% Construct the model object
tfe = tfeQCM('verbosity','none','dimension',theDimension);

%% Specify some stimuli
% Set up stim order info to create LM contrasts in
% several color directions.
%
% The magic kron has the effect of producing the matrix
% with LM for all contrasts in all directions.
maxContrast = 0.6;
contrastCoding = [0.0625, 0.125, 0.25, 0.5, 1];
directionCoding = maxContrast*[ [1 0]',  [1 1]',  [0 1]', [1 -1]' ]; 
stimulusStruct.values = kron(directionCoding',contrastCoding')';
nValues = size(stimulusStruct.values,2);

%% Dummy up stimulus struct 
%
% We're just using the tfe for convenience, there is
% no time here.  Dummy up a timebase.
stimulusStruct.timebase = 1:nValues;
nTimeSamples = size(stimulusStruct.timebase,2);

%% Set parameters and simulate responses
params1 = tfe.defaultParams;
params1.Qvec = [5 45];
params1.crfAmp = 2;
params1.crfSemi = 1;
params1.crfExponent = 2;
params1.noiseSd = 0.01;
params1.offset = 0;
fprintf('Simulated model parameters:\n');
tfe.paramPrint(params1);
modelResponseStruct = tfe.computeResponse(params1,stimulusStruct,[],'AddNoise',true);

%% Each color direction is plotted in increasing contrast order.
if (generatePlots)
    tfe.plot(modelResponseStruct);
end

%% Construct a packet
thePacket.stimulus = stimulusStruct;
thePacket.response = modelResponseStruct;
thePacket.kernel = [];
thePacket.metaData = [];

%% Fit
[paramsFit,fVal,fitResponseStruct] = tfe.fitResponse(thePacket);
fprintf('Model parameter from fits:\n');
tfe.paramPrint(paramsFit);

%% Plot fit on top of data
if (generatePlots)
    tfe.plot(fitResponseStruct,'Color',[0 1 0],'NewWindow',false);
end


