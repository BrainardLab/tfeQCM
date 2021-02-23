% Work out various features of a channel model
%
% Work out isorespone contour ideas for the linear channel model
% of Brouwer and Heeger.
%
% 2/16/21  dhb Finish up and add comments
% 2/22/21  dhb Modularize the functions.

%% Clear
clear; close all;

%% Set channel parameters
%
% channelType, choices are:
%   'BHSin2'           - Brouwer and Heeger 6 sin^2 channels
%   'OpponentSin2'     - Two opponent sin^2 channels
%   'OpponentLin2'     - Linear combination of cones, then squared
channelType = 'OpponentLin2';
switch (channelType)
    case 'BHSin2'
        nChannels = 6;
        startCenter = 0;
        theChannelWeightsPos = [0.2 0.4 0.6]';
        
    case 'OpponentSin2'
        nChannels = 4;
        startCenter = 45;
        theChannelWeightsPos = [0.2 0.8]';
        
    case 'OpponentLin2'
        nChannels = 4;
        coneContrastWeights = [ [1 1]' [1 -1]' [-1 -1]' [-1 1]' ];
        %coneContrastWeights = [ [1 0]' [0 1]' [-1 0]' [0 -1]' ];
   
    otherwise
        error('Unknown channelType specified');
end

% Check
if (rem(nChannels,2) ~= 0)
    error('nChannels must be even');
end

% Set criterion response
criterionResp = 1;

%% Create channel tuning
angleSupport = 0:1:360;
switch (channelType)
    case {'BHSin2', 'OpponentSin2'}
        % These have tuning described as half wave rectified sinusoids
        % squared. Compute response to unit contrast in each color
        % direction.
        %
        % Set channel center points
        centerSpacing = 360/nChannels;
        centerLocations = startCenter:centerSpacing:360-centerSpacing+startCenter;
        
        for ii = 1:nChannels
            underlyingChannels(ii,:) = cosd(angleSupport-centerLocations(ii));
            underlyingChannels(ii,sign(underlyingChannels(ii,:)) == -1) = 0;
            underlyingChannels(ii,:) = underlyingChannels(ii,:).^2;
        end
        
        % Just for fun, create a single linear channel as a weighted combination of the underlying channels
        %
        % Since we use symmetric modulations in our experiments, we'll keep
        % this symmetric.  We specify the first three weights and then duplicat
        % for the second three.
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
        
        % And then get isoresponse contour from the sensitivity we created
        % and also plot this.
        scalar = 1;
        theChannelIsoContrast = ChannelWeightsPosToIsoContrast([theChannelWeightsPos ; 1],underlyingChannels,criterionResp,channelType);
        subplot(1,2,2); hold on;
        plot(theChannelIsoContrast.*cosd(angleSupport),theChannelIsoContrast.*sind(angleSupport),'k','LineWidth',2);
        axis('square');
        xlabel('Cone 1 Contrast');
        ylabel('Cone 2 Contrast');
        title('Channel IsoContrast');
        
    case {'OpponentLin2'}
        % These are linear combinations of the cone contrasts.  Compute
        % response to unit contrast in each direction, half wave rectifid
        % and squared.
        for ii = 1:nChannels
            for jj = 1:length(angleSupport)
                cone1Contrast = cosd(angleSupport(jj));
                cone2Contrast = sind(angleSupport(jj));
                linResp = coneContrastWeights(1,ii)*cone1Contrast + coneContrastWeights(2,ii)*cone2Contrast;
                if (linResp < 0), linResp = 0; end
                underlyingChannels(ii,jj) = linResp;
            end
        end
    otherwise
        error('Unknown channelType specified');
end


% Plot channel sensitivities
figure; clf; hold on
colors = ['r' 'g' 'k' 'b' 'y' 'c'];
colorIndex = 1;
for ii = 1:nChannels
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
[isoContrastPred,weightsPosPred,isoContrastPred0] = FitIsoContrast(ellipticalIsoContrast,underlyingChannels,criterionResp,channelType);

% Get the starting isoresonse contrasts and add to plot.
%plot(isoContrastPred0.*cosd(angleSupport),isoContrastPred0.*sind(angleSupport),'g:','LineWidth',2);

plot(isoContrastPred.*cosd(angleSupport),isoContrastPred.*sind(angleSupport),'r:','LineWidth',2);
axis('square');
xlabel('Cone 1 Contrast');
ylabel('Cone 2 Contrast');
title('IsoContrast');

%% ChannelWeightsPosToIsoContrast
function isoContrast = ChannelWeightsPosToIsoContrast(theChannelWeightsPos,underlyingChannels,criterionResp,channelType)

isoContrast = ChannelWeightsToIsoContrast([theChannelWeightsPos(1:end-1) ; theChannelWeightsPos(1:end-1) ; theChannelWeightsPos(end)],underlyingChannels,criterionResp,channelType);

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
function isoContrast = ChannelWeightsToIsoContrast(theChannelWeights,underlyingChannels,criterionResp,channelType)

switch (channelType)
    case {'BHSin2', 'OpponentSin2'}
        
        % Form overall channel sensitivity as linear combination of
        % underlying channels.
        theChannelResp = theChannelWeights(end)*(underlyingChannels'*theChannelWeights(1:end-1))';
               
        % Avoid divide by zero
        theChannelResp(abs(theChannelResp) < 1e-6) = 1e-6;

        % Get isocontrast
        isoContrast = criterionResp./theChannelResp;

    case {'OpponentLin2'}
        nChannels = size(underlyingChannels,1);
        if (nChannels ~= 4)
            error('This model is specific for four single-ended channels in two pairs');
        end
        
        % Collapse to two channels
        twoChannels(1,:) = underlyingChannels(1,:)-underlyingChannels(3,:);
        twoChannels(2,:) = underlyingChannels(2,:)-underlyingChannels(4,:);
        
        % Compute the quadratic response.  This is the ellipse quadratic
        % form, with no constraint of positive definite.
        theChannelResp = (theChannelWeights(1)*twoChannels(1,:).^2) + ...
                         (theChannelWeights(2)*twoChannels(2,:).^2) + ...
                         (2*theChannelWeights(end)*twoChannels(1,:).*twoChannels(2,:));
                     
                           
        % Avoid divide by zero
        theChannelResp(abs(theChannelResp) < 1e-6) = 1e-6;

        % Get isocontrast
        isoContrast = sqrt(criterionResp./theChannelResp);
          
    otherwise
        error('Unknown channelType specified');
end


end

%% FitIsoContrast
%
% Use fmincon to find weights that produce isoresponse contour closest
% to the specified one.
function [isoContrastPred,weightsPosPred,isoContrastPred0] = FitIsoContrast(isoContrast,underlyingChannels,criterionResp,channelType)

% You need to make a good guess about the initial weights for the search to
% converge. So, find initial weights using regression.  We find all x
% weights, but only specify half for the search (as above, since we are
% using bipolar modulations).  Average coresponding weights to get starting
% point.  Have to pay a bit of attention to what's a row vector and what's
% a column vector through here.

% Set up base search bounds
nChannels = size(underlyingChannels,1);

% Slightly different heuristics used depending on model type.
switch (channelType)
    case {'BHSin2', 'OpponentSin2'}
        weightsRegress0 = ((underlyingChannels)'\(criterionResp./isoContrast'));
        weightsPos0 = zeros(nChannels/2,1);
        for ii = 1:nChannels/2
            weightsPos0(ii) = mean([weightsRegress0(ii) weightsRegress0(ii+nChannels/2)]);
        end
        vlb = -1e6*ones(nChannels/2,1);
        vub = 1e6*ones(nChannels/2,1);
        
        scalar0 = 1;
        weightsPos0 = [weightsPos0 ; scalar0];
        vlb = [vlb ; 1];
        vub = [vub ; 1];
    case {'OpponentLin2'}
        weightsRegress0 = ((underlyingChannels)'\(criterionResp./isoContrast'));
        weightsPos0 = zeros(nChannels/2,1);
        for ii = 1:nChannels/2
            weightsPos0(ii) = mean([weightsRegress0(ii) weightsRegress0(ii+nChannels/2)]);
        end
        %weightsPos0(2) = weightsPos0(2)*4;
        vlb = 0*ones(nChannels/2,1);
        vub = 1e6*ones(nChannels/2,1);
        
        scalar0 = 0;
        weightsPos0 = [weightsPos0 ; scalar0];
        vlb = [vlb ; -1e6];
        vub = [vub ; 1e6];
        % vlb = [vlb ; scalar0];
        % vub = [vub ; scalar0];
        
        Q0(1,1) = weightsPos0(1);
        Q0(2,2) = weightsPos0(2);
        Q0(1,2) = weightsPos0(3);
        Q0(2,1) = Q0(1,2);
        fprintf('Initial Q det: %g\n',det(Q0));
        

    otherwise
        error('Unknown channelType specified');
end

% Initial prediction (useful for debugging)
[~,isoContrastPred0] = FitIsoContrastFun(weightsPos0,underlyingChannels,isoContrast,criterionResp,channelType);

% Search
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','iter','LargeScale','off','Algorithm','active-set');
weightsPosPred = fmincon(@(weightsPos)FitIsoContrastFun(weightsPos,underlyingChannels,isoContrast,criterionResp,channelType),weightsPos0,[],[],[],[],vlb,vub,[],options);

Q(1,1) = weightsPosPred(1);
Q(2,2) = weightsPosPred(2);
Q(1,2) = weightsPosPred(3);
Q(2,1) = Q(1,2);
fprintf('Final Q det: %g\n',det(Q));

% Get prediction
[~,isoContrastPred] = FitIsoContrastFun(weightsPosPred,underlyingChannels,isoContrast,criterionResp,channelType);

end

%% FitIsoContrastFun
%
% Error function for search
function [f,isoContrastPred] = FitIsoContrastFun(theChannelWeightsPos,underlyingChannels,isoContrast,criterionResp,channelType)

% Get predicted isoresponse contour
isoContrastPred = ChannelWeightsPosToIsoContrast(theChannelWeightsPos,underlyingChannels,criterionResp,channelType);

% Take RMSE with what is being fit
f = 100*sqrt(sum((isoContrastPred-isoContrast).^2)/length(isoContrast));

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

