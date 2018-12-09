function modelResponseStruct = computeResponse(obj,params,stimulusStruct,kernelStruct,varargin)
% modelResponseStruct = computeResponse(obj,params,stimulusStruct,kernelStruct,varargin)
%
% Compute method for the quadratic color model with stimuli in
% direction/contrast form.
%
% Optional key/value pairs
%   'addNoise' - true/false (default false).  Add noise to computed
%                response?  Useful for simulations.

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
switch obj.dimension
    case 3
        directions = stimulusStruct.values(1:3,:);
        contrasts = stimulusStruct.values(4,:);
    case 2
        directions = stimulusStruct.values(1:2,:);
        contrasts = stimulusStruct.values(3,:);
end

%% Figure out passed directions and make sure they match up with directions
% this was initialized with.

[indDirectionsTemp,~,whichColumnsOut] = unique(stimulusStruct.values','rows','stable');  % uniquetol(directions',unique_tolerance,'ByRows',true);
indDirectionsTemp = indDirectionsTemp';
nIndDirections = size(indDirectionsTemp,2);
if (nIndDirections > obj.nDirections)
    error('Passed stimuli with more directions than we know about');
end

%% Match up directions found in stimuli with those that we have from initializtion
objDirectionIndices = zeros(nIndDirections,1);
for ii = 1:nIndDirections
    for jj = 1:obj.nDirections
        if (max(abs(indDirectionsTemp(:,ii) - obj.directions(jj,:))) < 1e-10)
            indDirectionsIntoObjDirectionsIndices(ii) = jj;
        end
    end
end
if (any(objDirectionIndices == 0))
    error('Passed stimulus direciton not in set of directions we know about');
end

%% Working here.  Each direction in object must be in stimulus and vice
% versa. Check this a bit more here.

%% Parse stimuli in terms of which stimulus belongs to which direction, and
% match up order with parameters for initialized directions
for ii = 1:nIndDirections
    whichColumns = find(whichColumnsOut == ii);
    indDirectionIndices{ii} = whichColumns;
    indDirectionResponses{ii} = responses(whichColumns);
    indDirectionDirections{ii} = indDirectionsTemp(:,ii);
    indDirectionContrasts{ii} = contrasts(whichColumns);
end

%% What I want to do is replace the code below with a call to the tfeQCM method.
% But I get an error about no matching signature when I try that. So for
% now just doing the same calculation here.
%
% neuralResponse = computeResponse@tfeQCM(params,stimuli);

%% Get neural response from QCM model
neuralResponse = tfeQCMForward(params,stimuli);

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
