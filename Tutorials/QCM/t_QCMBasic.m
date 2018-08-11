function validationData = t_QCMBasic(varargin)
% validationData = t_QCMBasic(varargin)
%
% Demonstrate some basic functionalty for the quadratic color model.
%
% Optional key/value pairs
%  'generatePlots' - true/fale (default true).  Make plots?

%% Parse vargin for options passed here
p = inputParser;
p.addParameter('generatePlots',true,@islogical);
p.parse(varargin{:});

%% Construct the model object
tfe = tfeQCM('verbosity','none');

%% Set parameters
%
% Six parameters define a quadratic form in three dimensions, but
% we normalize the first to 1 so we only need five numbers here.
params0 = tfe.defaultParams;
fprintf('Default model parameters:\n');
tfe.paramPrint(params0);

%% Set the timebase we want to compute on
deltaT = 1;
totalTime = 1000;

%% Specify the stimulus struct. 
%
% We'll specify this as a 3 by size(timebase,2) matrix,
% where each column is the signed L,M,S contrast of the stimulus
% at the specified time.  And then we'll blur it so that we have
% a smoothish signal to look at.
stimulusStruct.timebase = 0:deltaT:totalTime;
nTimeSamples = size(stimulusStruct.timebase,2);
filter = fspecial('gaussian',[1 nTimeSamples],6);
stimulusStruct.values = rand(3,nTimeSamples);
for i = 1:3
    stimulusStruct.values(i,:) = ifft(fft(stimulusStruct.values(i,:)) .* fft(filter)); 
end

%% Test that we can obtain a neural response
params1 = params0;
params1.crfAmp = 2;
params1.crfSemi = 0.5;
params1.crfExponent = 3;
params1.noiseSd = 0.02;
params1.offset = .5;
fprintf('Simulated model parameters:\n');
tfe.paramPrint(params1);
modelResponseStruct = tfe.computeResponse(params1,stimulusStruct,[],'AddNoise',true);
if (p.Results.generatePlots)
    tfe.plot(modelResponseStruct);
end

%% Construct a packet
thePacket.stimulus = stimulusStruct;
thePacket.response = modelResponseStruct;
thePacket.kernel = [];
thePacket.metaData = [];

%% Test the fitter
[paramsFit,fVal,fitResponseStruct] = tfe.fitResponse(thePacket);
fprintf('Model parameter from fits:\n');
tfe.paramPrint(paramsFit);
if (p.Results.generatePlots)
    tfe.plot(fitResponseStruct,'Color',[0 1 0],'NewWindow',false);
end

%% Set returned validationData structure
if (nargout > 0)
    validationData.params1 = params1;
    validationData.modelResponseStruct = modelResponseStruct;
    validationData.thePacket = thePacket;
end

end

