function isoContrast = tfeEllipticalIsoContrast(ellipseAngle,ellipseAspectRatio,angleSupport,criterionResp)
% Generate an elliptical isoresponse contour.
%
% Syntax:
%    isoContrast = tfeEllipticalIsoContrast(ellipseAngle,ellipseAspectRatio,angleSupport,criterionResp)
%
% Description:
%    Generate an elliptical isoresponse contour corresponding to the passed
%    major axis angle and major/minor axis ratio.  The contour is computed
%    on the angular support vector passed, and corresponds to response
%    of criterionResp when unit contrast vector is passed in, for the QCM
%    model.  See the appendix of Barnett et al.
%
%    This is a convenient way to generate an isoresponse contour from QCM,
%    either for plotting or for use in scaling LCM isoresponse contour.
%
%    It was pulled from earlier code we wrote to check
%    the appendix in the paper, t_EllipseCheck
%
% Inputs:
%    ellipseAngle       - Angle counterclockwise from x-axis in degrees of
%                         major axis of the ellipse.
%    ellipseAspectRatio - Ratio of length of major to minor axis of the
%                         ellipse. 
%    anglularSupport    - Angular support in degrees on which the contour
%                         is created. Often 1:1:360.
%    criterionResp      - Unit contrast in gives this value out when the
%                         square root of the ellipse quadratic form is computed.
%
% Outputs:
%    isoContrast       - Contrast in that produces square root of quadratic
%                        form value of criterionResp.
%
% Optional key/value pairs:
%     None.
%
% See also: t_EllipseCheck, t_LCMDirectionFit, t_ChannelModel


% Set up matrices V and S
V = [cosd(ellipseAngle) -sind(ellipseAngle) ; sind(ellipseAngle) cosd(ellipseAngle)];
S = [1 0 ; 0 1/ellipseAspectRatio];

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



