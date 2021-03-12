function modelResponseStruct = computeResponse(obj,params,stimulusStruct,kernelStruct,varargin)
% modelResponseStruct = computeResponse(obj,params,stimulusStruct,kernelStruct,varargin)
%
% Compute method for the LCM (linear channel model), with stimuli in
% direction/contrast form.
%
% Optional key/value pairs
%   'addNoise'         - true/false (default false).  Add noise to computed
%                        response?  Useful for simulations.

%% History
%    03/01/21  dhb  Wrote from QCM version.

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

% The conversion is done just once over a series a calls if params.fitting
% is true, because when we're fitting this routine gets called over
% and over for the same stimuli.  The caller handles setting and
% clearing the flag.
if ((obj.fitting & isempty(obj.angles)) | ~obj.fitting)
    params.fitting = obj.fitting;
    
    % Convert stimulus values to useful format
    directions = stimulusStruct.values(1:2,:);
    contrasts = stimulusStruct.values(3,:);
    angles = ones(size(contrasts));
    for ii = 1:size(directions,2)
        % Convert x,y direction into integer angle in degrees; one
        % degree accuracy seems sufficient and simplifies computation
        % of responses since we can just index into the mechanism
        % sensitivities.
        angles(ii) = round(atan2d(directions(2,ii),directions(1,ii)));
        
        % Wrap into 0 to 360
        while (angles(ii) <= 0)
            angles(ii) = angles(ii) + 360;
        end
        while (angles(ii) > 360)
            angles(ii) = angles(ii) - 360;
        end
        
        % Make sure we understand how direction is specified
        if (abs(norm(directions(:,ii))-1) > 1e-4)
            error('Direction vector not normalized as we expect');
        end
        obj.angles = angles;
        obj.contrasts = contrasts;
    end
elseif (obj.fitting)
    angles = obj.angles;
    contrasts = obj.contrasts;
else
    error('Logic error in code.  Should not ever get to this else statement');
end

% Create a single linear channel as a weighted combination
% of the underlying channels. 
%
% Since we use symmetric modulations in our experiments, we'll keep
% this symmetric.  We specify the first three weights and then duplicat
% for the second three.
theChannelWeights = [[params.channelWeightsPos]' ; [params.channelWeightsPos]'];
theChannel = (obj.underlyingChannels'*theChannelWeights)';

%% Get linear response from LCM params
%
% We count on fact that the angle support is integers between 1 and 360,a
% and that we've rounded angles to nearest degree.
unitSensitivities = theChannel(angles);
linearResponse = contrasts.*unitSensitivities;

% Deal with small value problem
linearResponse(linearResponse <= 1e-6) = 1e-6;

%% Push the LCM linear response through a Naka-Rushton non-linearity
%
% Only do this if NR parameters are there.  Otherwise return equivalent
% contrasts.
if (isfield(params,'crfAmp'))
    neuralResponse = ComputeNakaRushton([params.crfAmp,params.crfSemi,params.crfExponent],linearResponse) + params.crfOffset;
else
    neuralResponse = linearResponse;
end

%% Make the neural response structure
modelResponseStruct.timebase = stimulusStruct.timebase;
modelResponseStruct.values = neuralResponse;

%% Optionally, convolve with a passed kernel
modelResponseStruct = obj.applyKernel(modelResponseStruct,kernelStruct,varargin{:});

%% Optional add noise
if (p.Results.addNoise)
    modelResponseStruct.values = modelResponseStruct.values + normrnd(0,params.noiseSd,size(modelResponseStruct.values));
end


end
