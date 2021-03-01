% Work out various features of a channel model
%
% Work out isorespone contour ideas for the linear channel model
% of Brouwer and Heeger.
%
% 2/16/21  dhb Finish up and add comments
% 2/22/21  dhb Modularize the functions.
% 2/23/21  dhb Fix some bugs, finalize

%% Clear
clear; close all;

%% Set channel parameters
%
% channelType, choices are:
%   'BH_6_Sin2'        - Brouwer and Heeger 6 sin^2 channels, which respond linearly.
%   'BH_6_Lin2'        - Linear then squared version of BH_6_Sin2 model
%   'BH_4_Sin2'        - Brouwer and Heeger 4 sin^2 channels, which respond linearly.
%   'BH_4_Lin2'        - Linear then squared version of BH_4_Sin2 model.
%   'BH_8_Sin2'        - Brouwer and Heeger 8 sin^2 channels, which respond linearly.
%   'BH_8_Lin2'        - Linear then squared version of BH_8_Sin2 model
%   'OpponentCone2'    - Linear combination of cones, outpus then squared and added
channelType = 'BH_8_Lin2';
switch (channelType)
    case {'BH_6_Sin2', 'BH_6_Lin2'}
        nChannels = 6;
        startCenter = 0;
        theChannelWeightsPos = [0.2 0.4 0.6]';
        
    case {'BH_4_Sin2', 'BH_4_Lin2'}
        nChannels = 4;
        startCenter = 45;
        theChannelWeightsPos = [0.2 0.8]';
        
    case {'BH_8_Sin2', 'BH_8_Lin2'}
        nChannels = 8;
        startCenter = 0;
        theChannelWeightsPos = [0.2 0.8 0.3 0.1]';
        
    case 'OpponentCone2'
        nChannels = 4;
        coneContrastWeights = [ [1 1]' [1 -1]' [-1 -1]' [-1 1]' ];
        %coneContrastWeights = [ [1 0]' [0 1]' [-1 0]' [0 -1]' ];
   
    otherwise
        error('Unknown channelType specified');
end

%% Ellipse parameters
angle = 45;         % Angle
aspectRatio = 0.15; % Minor axis aspect ratio

% Check
if (rem(nChannels,2) ~= 0)
    error('nChannels must be even');
end

% Set criterion response
criterionResp = 2;

%% Create channel tuning
angleSupport = 0:1:360;
switch (channelType)
    case {'BH_6_Sin2', 'BH_6_Lin2', 'BH_4_Sin2', 'BH_4_Lin2', 'BH_8_Sin2', 'BH_8_Lin2'}
        % These have tuning described as half wave rectified sinusoids
        % squared. Compute response to unit contrast in each color
        % direction by regarding these as linear channels tuned to hue
        % angle.
        %
        % Set channel center points
        centerSpacing = 360/nChannels;
        centerLocations = startCenter:centerSpacing:360-centerSpacing+startCenter;
        
        for ii = 1:nChannels
            underlyingChannels(ii,:) = cosd(angleSupport-centerLocations(ii));
            underlyingChannels(ii,sign(underlyingChannels(ii,:)) == -1) = 0;
            switch (channelType)
                case {'BH_6_Sin2', 'BH_4_Sin2', 'BH_8_Sin2'}
                    underlyingChannels(ii,:) = underlyingChannels(ii,:).^2;
            end
        end
        
        % Just for fun, create a single linear channel as a weighted combination
        % of the underlying channels.  This is mainly to look at a sample
        % isoresponse contour created from this model.
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
        [theChannelIsoContrast,theChannelIsoResp] = ChannelWeightsPosToIsoContrast(theChannelWeightsPos,angleSupport,underlyingChannels,criterionResp,channelType);
        if (max(abs(theChannelIsoResp-criterionResp)) > 1e-6)
            error('Fit iso ontrast does not produce criterion response');
        end
        subplot(1,2,2); hold on;
        plot(theChannelIsoContrast.*cosd(angleSupport),theChannelIsoContrast.*sind(angleSupport),'k','LineWidth',2);
        axis('square');
        xlabel('Cone 1 Contrast');
        ylabel('Cone 2 Contrast');
        title('Channel IsoContrast');
        
    case {'OpponentCone2'}
        % These are linear combinations of the cone contrasts.  Compute
        % response to unit contrast in each direction, half wave rectified
        % but not squared.  For this thread, we do the squaring as part of
        % the linear combination across channels.
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

%% Plot channel sensitivities
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
ellipticalIsoContrast = EllipticalIsoContrast(angle,aspectRatio,angleSupport,criterionResp);

% Plot
figure; clf; hold on;
plot(ellipticalIsoContrast.*cosd(angleSupport),ellipticalIsoContrast.*sind(angleSupport),'k','LineWidth',2);
axis('square');

%% If we have an isoresponse contour, we can solve for the channel weights
[ellipticalIsoContrastFit,ellipticalIsoReponseFit,weightsPosPred,isoContrastPred0] = FitIsoContrast(ellipticalIsoContrast,angleSupport,underlyingChannels,criterionResp,channelType);
if (max(abs(ellipticalIsoReponseFit-criterionResp)) > 1e-6)
    error('Fit iso ontrast does not produce criterion response');
end

% Get the starting isoresonse contrasts and add to plot.
%plot(isoContrastPred0.*cosd(angleSupport),isoContrastPred0.*sind(angleSupport),'g:','LineWidth',2);

% Add resulting isoresponse contour to the plot
plot(ellipticalIsoContrastFit.*cosd(angleSupport),ellipticalIsoContrastFit.*sind(angleSupport),'r:','LineWidth',2);
axis('square');
xlim([-2 2]); ylim([-2 2]);
xlabel('Cone 1 Contrast');
ylabel('Cone 2 Contrast');
title('IsoContrast');

%% ChannelWeightsPosToIsoContrast
function [isoContrast,isoResp] = ChannelWeightsPosToIsoContrast(theChannelWeightsPos,angleSupport,underlyingChannels,criterionResp,channelType)

[isoContrast,isoResp] = ChannelWeightsToIsoContrast([theChannelWeightsPos ; theChannelWeightsPos],angleSupport,underlyingChannels,criterionResp,channelType);

end

%% ChannelWeightsToIsoContrast
function [isoContrast,isoResp] = ChannelWeightsToIsoContrast(theChannelWeights,angleSupport,underlyingChannels,criterionResp,channelType)

% Get channel response to unit contrast in each angular direction
unitContrast = ones(size(angleSupport));
theChannelResp = ChannelWeightsToChannelResponse(theChannelWeights,unitContrast,underlyingChannels,channelType);

% Get isocontrast.  The sin^2 tuned channels respond linearly with
% contrast, so for each stimulus direction response is given
% by stimulus magnitude times sensitivity in that
% direction. That is:
%   resp = theConrast.*theChannel*stimMag
% So for stimuli in just one direction theta, the contrast to produce a
% given response is:
%   theContrast = resp/theChannel(theta);
%
% We can get the stimulus magnitude for the channel
% above by inverting this equation.
%
% The linear squared channels are a little different, but we do return the
% square root of the summed squared response for those, so that response
% is linearly proportional to contrast in each direction.  That means the
% same formula applies.
isoContrast = criterionResp./theChannelResp;
isoResp = ChannelWeightsToChannelResponse(theChannelWeights,isoContrast,underlyingChannels,channelType);

end

%% ChannelWeightsPosToChannelResponse
function [theChannelResp] = ChannelWeightsPosToChannelResponse(theChannelWeightsPos,theContrast,underlyingChannels,channelType)

% Expand weights full set of channels and pass to full routine.  Need to be
% a little careful about the extra parameter that is tacked onto the end.
theChannelResp = ChannelWeightsToChannelResponse([theChannelWeightsPos ; theChannelWeightsPos],theContrast,underlyingChannels,channelType);

end

%% ChannelWeightsToChannelResponse
function [theChannelResp] = ChannelWeightsToChannelResponse(theChannelWeights,theContrast,underlyingChannels,channelType)
    
switch (channelType)
    case {'BH_6_Sin2', 'BH_4_Sin2', 'BH_8_Sin2'}
        
        % Form overall channel sensitivity as linear combination of
        % underlying channels.
        theChannelResp = theContrast .* (underlyingChannels'*theChannelWeights)';

    case {'OpponentCone2', 'BH_6_Lin2', 'BH_4_Lin2', 'BH_8_Lin2'}
        % Collapse underlying channels.  This just turns the rectified linear
        % channels back into full (unrectified) linear channels.
        nChannels = size(underlyingChannels,1);
        for ii = 1:nChannels/2
            linChannels(ii,:) = underlyingChannels(ii,:)-underlyingChannels(ii+nChannels/2,:);
        end

        % Compute the quadratic response.  This is the ellipse quadratic
        % form, with no constraint of positive definite.
        theChannelResp2 = 0;
        for ii = 1:nChannels/2
            theChannelResp2 = theChannelResp2 + (theChannelWeights(ii)*(theContrast.*linChannels(ii,:)).^2);
        end
                     
        % Take square root at the end so that response is proportional to
        % contrast.
        theChannelResp = sqrt(theChannelResp2);
      
    otherwise
        error('Unknown channelType specified');
end
                       
% Avoid divide by zero in later calculations
theChannelResp(abs(theChannelResp) < 1e-6) = 1e-6;

end

%% FitIsoContrast
%
% Use fmincon to find weights that produce isoresponse contour closest
% to the specified one.
function [isoContrastPred,isoRespPred,weightsPosPred,isoContrastPred0] = FitIsoContrast(isoContrast,angleSupport,underlyingChannels,criterionResp,channelType)

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
    case {'BH_6_Sin2', 'BH_4_Sin2', 'BH_8_Sin2'}
        weightsRegress0 = ((underlyingChannels)'\(criterionResp./isoContrast'));
        weightsPos0 = zeros(nChannels/2,1);
        for ii = 1:nChannels/2
            weightsPos0(ii) = mean([weightsRegress0(ii) weightsRegress0(ii+nChannels/2)]);
        end
        vlb = -1e6*ones(nChannels/2,1);
        vub = 1e6*ones(nChannels/2,1);
     
    case {'OpponentCone2', 'BH_6_Lin2', 'BH_4_Lin2', 'BH_8_Lin2'}
        weightsRegress0 = ((underlyingChannels)'\(criterionResp./isoContrast'));
        weightsPos0 = zeros(nChannels/2,1);
        for ii = 1:nChannels/2
            weightsPos0(ii) = mean([weightsRegress0(ii) weightsRegress0(ii+nChannels/2)]);
        end
        %weightsPos0(2) = weightsPos0(2)*4;
        vlb = 0*ones(nChannels/2,1);
        vub = 1e6*ones(nChannels/2,1);
       
    otherwise
        error('Unknown channelType specified');
end

% Initial prediction (useful for debugging)
[~,isoContrastPred0] = FitIsoContrastFun(weightsPos0,angleSupport,underlyingChannels,isoContrast,criterionResp,channelType);

% Search
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','active-set');
weightsPosPred = fmincon(@(weightsPos)FitIsoContrastFun(weightsPos,angleSupport,underlyingChannels,isoContrast,criterionResp,channelType),weightsPos0,[],[],[],[],vlb,vub,[],options);

% Get prediction
[~,isoContrastPred,isoRespPred] = FitIsoContrastFun(weightsPosPred,angleSupport,underlyingChannels,isoContrast,criterionResp,channelType);

end

%% FitIsoContrastFun
%
% Error function for search
function [f,isoContrastPred,isoRespPred] = FitIsoContrastFun(theChannelWeightsPos,angleSupport,underlyingChannels,isoContrast,criterionResp,channelType)

% Get predicted isoresponse contour
[isoContrastPred,isoRespPred] = ChannelWeightsPosToIsoContrast(theChannelWeightsPos,angleSupport,underlyingChannels,criterionResp,channelType);

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

