function modelResponseStruct = computeResponse(obj,params,stimulusStruct,kernelStruct,varargin)
% modelResponseStruct = computeResponse(obj,params,stimulusStruct,kernelStruct,varargin)
%
% Compute method for the quadratic color model with stimuli in
% direction/contrast form.
%
% In addition to timebase and values fields, the returned response structure contains a metaData field.
%   modelResponseStruct.metaData,quadtraticFactors.
% These factors may be used to adjust the Naka-Rushton parameters per color
% direction.  See tfeQCMForward.
%
% Optional key/value pairs
%   'addNoise'         - true/false (default false).  Add noise to computed
%                        response?  Useful for simulations.

%% Parse input
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true;
p.addRequired('params',@isstruct);
p.addRequired('stimulusStruct',@isstruct);
p.addRequired('kernelStruct',@(x)(isempty(x) || isstruct(x)));
p.addParameter('addNoise',false,@islogical);
p.parse(params,stimulusStruct,kernelStruct,varargin{:});
params = p.Results.params;

%% Convert stimulus values to useful format
%
% Handle stimulus dimension
switch obj.dimension
    case 3
        directions = stimulusStruct.values(1:3,:);
        contrasts = stimulusStruct.values(4,:);
    case 2
        directions = stimulusStruct.values(1:2,:);
        contrasts = stimulusStruct.values(3,:);
end

% The converstion is done just once over a series a calls if params.fitting
% is true, because when we're fitting this routine gets called over
% and over for the same stimuli.  The caller handles setting and
% clearing the flag.
if (obj.fitting & isempty(obj.stimuli))
    params.fitting = true;
    stimuli = tfeQCMDirectionsContrastsToStimuli(directions,contrasts);
    obj.stimuli = stimuli;
elseif (obj.fitting)
    stimuli = obj.stimuli;
else
    params.fitting = false;
    stimuli = tfeQCMDirectionsContrastsToStimuli(directions,contrasts);
end

%% What I want to do is replace the code below with a call to the tfeQCM method.
% But I get an error about no matching signature when I try that. So for
% now just doing the same calculation here.
%
% stimulusStruct.values = stimuli;
% neuralResponse = computeResponse@tfeQCM(params,stimulusStruct,kernelStruct,varargin{:});

%% Get neural response from QCM model
[neuralResponse,quadraticFactors] = tfeQCMForward(params,stimuli);

%% Make the neural response structure
modelResponseStruct.timebase = stimulusStruct.timebase;
modelResponseStruct.values = neuralResponse;
modelResponseStruct.metaData.quadraticFactors = quadraticFactors;

%% Optionally, convolve with a passed kernel
modelResponseStruct = obj.applyKernel(modelResponseStruct,kernelStruct,varargin{:});

%% Optional add noise
if (p.Results.addNoise)
    modelResponseStruct.values = modelResponseStruct.values + normrnd(0,params.noiseSd,size(modelResponseStruct.values));
end


end
