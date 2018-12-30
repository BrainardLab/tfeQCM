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
[indDirectionsTemp,directionIndices] = tfeQCMParseDirections(directions);
nIndDirections = size(indDirectionsTemp,2);
if (nIndDirections ~= obj.nDirections)
    error('Passed stimulus array does not have same number of directions as object');
end

%% Check on parameters
for ii = 1:nIndDirections
    if (abs(params(ii).noiseSd - params(1).noiseSd) > 1e-6)
        error('Noise sd parameter not matched across directions in parameters struct array');
    end
    if (abs(params(ii).expFalloff - params(1).expFalloff) > 1e-6)
        error('Exp falloff parameter not matched across directions in parameters struct array');
    end
end

%% Match up directions found in stimuli with those that we have from initializtion
objDirectionIndices = zeros(nIndDirections,1);
for ii = 1:nIndDirections
    for jj = 1:obj.nDirections
        if (max(abs(indDirectionsTemp(:,ii) - obj.directions(:,jj))) < 1e-10)
            objDirectionIndices(ii) = jj;
        end
    end
end
if (any(objDirectionIndices == 0))
    error('A passed stimulus direction is not in set of directions we know about');
end
checkIndices = unique(objDirectionIndices);
if (length(checkIndices) ~= obj.nDirections)
    error('One direction in object set not in passed stimulus array')
end

%% Parse stimuli in terms of which stimulus belongs to which direction
%
% Match up order with parameters for initialized directions
for ii = 1:nIndDirections
    theDirection = objDirectionIndices(ii);
    indDirectionIndices{theDirection} = directionIndices{ii};
    indDirectionDirections{theDirection} = indDirectionsTemp(:,ii);
    indDirectionContrasts{theDirection} = contrasts(directionIndices{ii});
end

%% Get neural response from NR forward model
neuralResponseCell = tfeNRForward(params,indDirectionContrasts);

% Convert cell array for each direction back into form that matches passed
% input.
for ii = 1:nIndDirections
    neuralResponse(indDirectionIndices{ii}) = neuralResponseCell{ii};
end

%% Make the neural response structure
modelResponseStruct.timebase = stimulusStruct.timebase;
modelResponseStruct.values = neuralResponse;

%% Optional, convolve with a passed kernel
modelResponseStruct = obj.applyKernel(modelResponseStruct,kernelStruct,varargin{:});

%% Optional, add noise
if (p.Results.addNoise)
    modelResponseStruct.values = modelResponseStruct.values + normrnd(0,params(1).noiseSd,size(modelResponseStruct.values));
end


end
