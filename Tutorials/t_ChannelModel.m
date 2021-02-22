% Work out various features of a channel model
%
% Work out isorespone contour ideas for the linear channel model
% of Brouwer and Heeger.
%
% 2/16/21  dhb Finish up and add comments
% 2/22/21  dhb Modularize the functions.

%% Clear
clear; close all;

%% Set number of channels
nChannels = 6;
if (rem(nChannels,2) ~= 0)
    error('nChannels must be even');
end

%% Set channel center points
centerSpacing = 360/nChannels;
centerLocations = 0:centerSpacing:360-centerSpacing;

%% Create and plot underlying channel tuning
angleSupport = 0:1:360;
figure; clf; hold on
colors = ['r' 'g' 'k' 'b' 'y' 'c'];
colorIndex = 1;
for ii = 1:nChannels
    underlyingChannels(ii,:) = cosd(angleSupport-centerLocations(ii));
    underlyingChannels(ii,sign(underlyingChannels(ii,:)) == -1) = 0;
    underlyingChannels(ii,:) = underlyingChannels(ii,:).^2;
    plot(angleSupport,underlyingChannels(ii,:),colors(colorIndex),'LineWidth',2);
    colorIndex = colorIndex+1;
    if (colorIndex > length(colors))
        colorIndex = 1;
    end
end
xlim([0 360]);
xlabel('Angle');
ylabel('Sensitivity');
title('Channel Sensitivities');

%% Create a single linear channel as a weighted combination of the underlying channels
%
% Since we use symmetric modulations in our experiments, we'll keep
% this symmetric.  We specify the first three weights and then duplicat
% for the second three.
rng(3);
theChannelWeightsPos = rand(nChannels/2,1);
theChannelWeights = [theChannelWeightsPos ; theChannelWeightsPos];
theChannel = (underlyingChannels'*theChannelWeights)';

% Plot sensitivity of the single channel as a function of angle
figure; clf;
subplot(1,2,1); hold on;
plot(angleSupport,theChannel,'k','LineWidth',2);
xlim([0 360]);
xlabel('Angle');
ylabel('Sensitivity');
title('Channel Sensitivity');

%% Get isoresponse contour from sensitivity
%
% If you wanted to plot this isoresponse contour, you could using the same
% code as we do to plot the elliptical and fit isoresponse contours below.
criterionResp = 2;
theChannelIsoContrast = ChannelWeightsPosToIsoContrast(theChannelWeightsPos,underlyingChannels,criterionResp);
subplot(1,2,2); hold on;
plot(theChannelIsoContrast.*cosd(angleSupport),theChannelIsoContrast.*sind(angleSupport),'k','LineWidth',2);
axis('square');
xlabel('Cone 1 Contrast');
ylabel('Cone 2 Contrast');
title('Channel IsoContrast');

%% Generate an elliptical isoresponse contour and plot
%
% Set ellipse parameters
angle = 45;         % Angle
aspectRatio = 0.15; % Minor axis aspect ratio
ellipticalIsoContrast = EllipticalIsoContrast(angle,aspectRatio,angleSupport,criterionResp);

% Plot
figure; clf; hold on;
plot(ellipticalIsoContrast.*cosd(angleSupport),ellipticalIsoContrast.*sind(angleSupport),'k','LineWidth',2);
axis('square');

%% If we have an isoresponse contour, we can solve for the channel weights
[isoContrastPred,weightsPosPred] = FitIsoContrast(ellipticalIsoContrast,underlyingChannels,criterionResp);

% Get the predicted isoresonse contrasts and add to plot.
plot(isoContrastPred.*cosd(angleSupport),isoContrastPred.*sind(angleSupport),'r:','LineWidth',2);
axis('square');
xlabel('Cone 1 Contrast');
ylabel('Cone 2 Contrast');
title('IsoContrast');

%% ChannelWeightsPosToIsoContrast
function isoContrast = ChannelWeightsPosToIsoContrast(theChannelWeightsPos,underlyingChannels,criterionResp)

isoContrast = ChannelWeightsToIsoContrast([theChannelWeightsPos ; theChannelWeightsPos],underlyingChannels,criterionResp);

end

%% ChannelWeightsToIsoContrast
%
% For each stimulus direction response is given
% by stimulus magnitude times sensitivity in that
% direction. That is:
%   resp = theChannel*stimMag
% So for stimuli in just one direction theta,%
%   theChannelIsoContrast = resp/theChannel(theta);
%
% We can get the stimulus magnitude for the channel
% above by inverting this equation.
function isoContrast = ChannelWeightsToIsoContrast(theChannelWeights,underlyingChannels,criterionResp)

% Form channel
theChannel = (underlyingChannels'*theChannelWeights)';

% Avoid divide by zero
theChannel(abs(theChannel) < 1e-6) = 1e-6;

% Get isocontrast
isoContrast = criterionResp./theChannel;

end

%% FitIsoContrast
%
% Use fmincon to find weights that produce isoresponse contour closest
% to the specified one.
function [isoContrastPred,weightsPosPred] = FitIsoContrast(isoContrast,underlyingChannels,criterionResp)

% You need to make a good guess about the initial weights for the search to
% converge. So, find initial weights using regression.  We find all x
% weights, but only specify half for the search (as above, since we are
% using bipolar modulations).  Average coresponding weights to get starting
% point.  Have to pay a bit of attention to what's a row vector and what's
% a column vector through here.
nChannels = size(underlyingChannels,1);
weightsRegress0 = (underlyingChannels'\(1./isoContrast'));
weightsPos0 = zeros(nChannels/2,1);
for ii = 1:nChannels/2
    weightsPos0(ii) = mean([weightsRegress0(ii) weightsRegress0(ii+nChannels/2)]);
end

% Set up search bounds
vlb = [-1e6 -1e6 -1e6]';
vub = [1e6 1e6 1e6]';

% Search
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','active-set');
weightsPosPred = fmincon(@(weightsPos)FitIsoContrastFun(weightsPos,underlyingChannels,isoContrast,criterionResp),weightsPos0,[],[],[],[],vlb,vub,[],options);

% Get prediction
[~,isoContrastPred] = FitIsoContrastFun(weightsPosPred,underlyingChannels,isoContrast,criterionResp);

end

%% FitIsoContrastFun
%
% Error function for search
function [f,isoContrastPred] = FitIsoContrastFun(theChannelWeightsPos,underlyingChannels,isoContrast,criterionResp)

% Get predicted isoresponse contour
isoContrastPred = ChannelWeightsPosToIsoContrast(theChannelWeightsPos,underlyingChannels,criterionResp);

% Take RMSE with what is being fit
f = sqrt(sum((isoContrastPred-isoContrast).^2)/length(isoContrast));

end

%% EllipticalIsoContrast
%
% Generate an elliptical isoresponse contour.
%
% This is pulled from earlier code we developed to check
% the appendix in the paper, t_EllipseCheck
function isoContrast = EllipticalIsoContrast(angle,aspectRatio,angleSupport,criterionResp)

% Set up matrices V and S
V = [cosd(angle) -sind(angle) ; sind(angle) cosd(angle)];
S = [1 0 ; 0 1/aspectRatio];

% Compute A and Q
A = S'*V';
Q = A'*A;

% Points on the ellipse satisfy c'*Q*c = resp^2;
% One way to find the ellipse is to go around the circle and adjust to desired length
resp2 = criterionResp^2;
circleVecs = zeros(2,length(angleSupport));
ellipseVecs = zeros(2,length(angleSupport));
for tt = 1:length(angleSupport)
    % Create a point on a circle
    circleVecs(1,tt) = cosd(angleSupport(tt));
    circleVecs(2,tt) = sind(angleSupport(tt));
    
    % Transform by scaling to produce a point on the ellipse
    rawResp2 = circleVecs(:,tt)'*Q*circleVecs(:,tt);
    ellipseVecs(:,tt) = sqrt(resp2)*circleVecs(:,tt)/sqrt(rawResp2);
    
    % Get vector length of point on ellipse (aka contrast).  This is
    % the contrast at each stimulus direction angle that produces
    % the constant criterion response.
    isoContrast(tt) = norm(ellipseVecs(:,tt));
end

end

