function modelResponseStruct = computeResponse(obj,params,stimulusStruct,kernelStruct,varargin)
% modelResponseStruct = computeResponse(obj,params,stimulusStruct,kernelStruct,varargin)
%
% Compute method for the quadratic color model.
%
% Optional key/value pairs
%   'addNoise' - true/false (default false).  Add noise to computed
%                response?  Useful for simulations.
%
% The field 'noiseInverseFrequencyPower' can be set in the parameters, and
%   will dictate the form of the stimulated noise. Noise is generated using
%   the dsp.ColoredNoise. This parameter determines alpha, the
%   exponent of the noise generation property, following the formula:
%   noise = 1/|f|^alpha. The default value of zero produces white noise.
%   A value of 1 produces 'pink' noise, and a value of 2 'brown'(ian)
%   noise. Functional MRI data have autocorrelated noise, with an
%   alpha value on the order of 1 or less.

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

%% Fitting?
if (obj.fitting)
    params.fitting = true;
else
    params.fitting = false;
end

%% Get neural response from QCM model
neuralResponse = tfeQCMForward(params,stimulusStruct.values);

%% Make the neural response structure
modelResponseStruct.timebase = stimulusStruct.timebase;
modelResponseStruct.values = neuralResponse;

%% Optionally, convolve with a passed kernel
modelResponseStruct = obj.applyKernel(modelResponseStruct,kernelStruct,varargin{:});

%% Optional add noise
if p.Results.addNoise
    if ~isfield(params, 'noiseInverseFrequencyPower')
        params.noiseInverseFrequencyPower=0;
    end
    cn = dsp.ColoredNoise(params.noiseInverseFrequencyPower,length(modelResponseStruct.timebase),1);
    noise = step(cn)';
    noise = noise/std(noise)*params.noiseSd;
    modelResponseStruct.values = modelResponseStruct.values + noise;
end

end
