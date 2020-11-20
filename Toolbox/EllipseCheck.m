%% Check that our ellipse code matches the write up and the analysis code

%% History
%    11/19/20  dhb  Wrote it.

%% Initialize
clear; close all;

%% Set ellipse parameters
p1 = 45;   % Angle of major axis, counterclockwise from x-axis
p2 = 0.25; % Minor axis aspect ratio, relative to major.

%% Set up matrix V.
%
% We use deg2rotm in the analysis code, here we write the matrix explicitly
% and verify that deg2rotm gives the same answer.  This matrix rotates
% vectors by the specified angle in degrees, counterclockwise.
V = [cosd(p1) -sind(p1) ; sind(p1) cosd(p1)];
V2 = deg2rotm(p1);
if (max(abs(V(:)-V2(:)) > 1e-6))
    error('Two ways of defining rotation matrix do not match');
end

% Verify counterclockwise rotation convention of V.  This is only checked
% when p1 == 45.
if (p1 == 45)
    check = V*[1 0]';
    if (abs(check(1)-sqrt(2)/2) > 1e-6 | abs(check(2)-sqrt(2)/2) > 1e-6)
        error('We do not understand rotation matrix convention');
    else
        fprintf('Rotation matrix has expected counterclockwise convention\n');
    end
else
    fprintf('Not checking rotation matrix convention. Set p1 = 45 to execute this check\n');
end

%% Set up S
%
% Using 1/p2 means that a p2 less than 1 shrinks the minor axis relative to
% the major axis length of 1.
S = [1 0 ; 0 1/p2];

%% Compute A and Q
%
% A is what we call M in the model appendix.
A = S'*V';
Ainv = inv(A);
Q = A'*A;

%% Points on the ellipse satisfy c'*Q*c = k^2.
% 
% The value of k is the equivalent contrast.  
k2 = 4;
nTheta = 1000;
circleVecs = UnitCircleGenerate(nTheta);

%% One way to find the ellipse is to go around the circle and adjust to desired length
%
% Do that here in each direction. We find the output of the quadratic form
% and scale the vector so that it will will go through the quadratic form
% and result in a value of k^2.  This means that in all directions, we get
% a vector that lies on the ellipse with equivalent contrast k.
ellipseVecs = zeros(size(circleVecs));
for tt = 1:nTheta
    rawK2 = circleVecs(:,tt)'*Q*circleVecs(:,tt);
    ellipseVecs(:,tt) = sqrt(k2)*circleVecs(:,tt)/sqrt(rawK2);
    equivalentContrasts(tt) = sqrt(ellipseVecs(:,tt)'*Q*ellipseVecs(:,tt));
end

%% Check lengths
if (max(abs(equivalentContrasts(:)-sqrt(k2))) > 1e-6)
    error('Do not get expected vector length for equivalent contrasts');
end

%% Another way is to apply the inverse of A to the circle
%
% You can't easily check the comparison between this way and the above
% numerically because the order of the points around the ellipse doesn't
% match up.  But you can see that they overlay in the plot below.  The
% factor of sqrt(k2) scales up the circle, to match the scale of the
% desired equivalent contrast.
ellipseVecs2 = sqrt(k2)*Ainv*circleVecs;

%% Now let's do it using the code in the tfeQCM
% 
% This takes minor axis ratio and angle as parameters
params.Qvec = [p2 p1];
[tfeEquivalentContrasts] = tfeQCMForward(params,ellipseVecs);
if (max(abs(tfeEquivalentContrasts-sqrt(k2) > 1e-6)))
    error('tfeQCM does not return expected equivalent contrasts for the isoresponse ellipse');
end

%% Plot ellipses both ways.
theLim = 4;
figure; clf; hold on
plot(ellipseVecs(1,:),ellipseVecs(2,:),'ro','MarkerSize',10,'MarkerFaceColor','r');
plot(ellipseVecs2(1,:),ellipseVecs2(2,:),'ko','MarkerSize',6,'MarkerFaceColor','k');
%plot(ellipseVecs3(1,:),ellipseVecs3(2,:),'go','MarkerSize',8,'MarkerFaceColor','g');

xlim([-theLim theLim]); ylim([-theLim theLim]);
axis('square');
    

