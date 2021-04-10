function modelResponseStruct = computeResponse(obj,params,stimulusStruct,kernelStruct,varargin)
% Compute LCM model responses
%
% Syntax:
%    modelResponseStruct = computeResponse(obj,params,stimulusStruct,kernelStruct,varargin)
%
% Description:
%   Compute method for the LCM (linear channel model), with stimuli in
%   direction/contrast form.
%
%   Note that we've implemented a more general form with a summation
%   exponent, which usually we set to 1 to match the LCM idea.
%
% Optional key/value pairs
%   'addNoise'         - true/false (default false).  Add noise to computed
%                        response?  Useful for simulations.
%   'noNakaRushton'    - true/false (default false). If true, don't apply
%                        Naka-Rushton, just return equivalent contrast.

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
p.addParameter('noNakaRushton',false,@islogical);
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
theChannel = ((obj.underlyingChannels)'*theChannelWeights)';

%% Get linear response from LCM params
%
% We count on fact that the angle support is integers between 1 and 360,a
% and that we've rounded angles to nearest degree.
theResponse1 = contrasts.*theChannel(angles);

%% This is another way of doing the same calculation when summation exponent is 1
%
% But it's what we want when that exponent is not 1.
theResponse2 = zeros(size(angles));
for cc = 1:size(obj.underlyingChannels,1)
    theResponse2 = theResponse2 + ((contrasts.*obj.underlyingChannels(cc,angles)).^(obj.summationExponent))*theChannelWeights(cc);
end

% Deal with small value and negative value problem
theResponse1(theResponse1 <= 1e-6) = 1e-6;
theResponse2(theResponse2 <= 1e-6) = 1e-6;

%% Check that both ways give same answer when summation exponent is 1
%
% For reasons I don't understand, choosing one or the other can matter
% even when they are the same to 1e-10. This has to do with fmincon
% and its evaluation of stopping point.  Using a different error function
% scale factor makes the two versions work the same. Keeping old way
% when summationExponent is 1, for backwards compatibility.
if (obj.summationExponent == 1)
    if (max(abs(theResponse1-theResponse2)./theResponse1) > 1e-8)
        error('Don''t compute same response two ways');
    end
    theResponse = theResponse1;
else
    theResponse = theResponse2;
end

%% Take the summed response and raise to 1/summationExponent
theResponse = theResponse.^(1/obj.summationExponent);

%% Push the LCM response through a Naka-Rushton non-linearity
%
% Only do this if NR parameters are there.  Otherwise return equivalent
% contrasts.
if (isfield(params,'crfAmp') & ~p.Results.noNakaRushton)
    neuralResponse = ComputeNakaRushton([params.crfAmp,params.crfSemi,params.crfExponent],theResponse) + params.crfOffset;
else
    neuralResponse = theResponse;
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
