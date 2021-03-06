classdef tfeLCMDirection < tfeQCM
% tfeLCMDirection
%  Create at tfeLCMDirection (Linear Channels Model) object
%
% Syntax:
%    tfe = tfeLCMDirection(varargin);
% 
% Description:
%     Implements a model where the response is a linear sum of angularly tuned
%     mechanisms, uses this to compute equivalent contrast, and then applies
%     Naka-Rushton function.  This is a generalization of the forward model
%     introduced by Brouwer and Heeger (2009, J. Neuro.,
%     29(44):13992–14003).
%
%     Stimulus specified as unit vector direction and contrast, and are
%     assumed to be bipolor modulations, so that we symmetrerize channel
%     sensitivities and weights.
%
%     Inherits optional key/value pairs from parent class tfeQCM, plus
%     those specified below.
%
%     Using a summation exponent other than 1 turns the linear sum into a
%     non-linear sum.  We have a check that if summation exponent isn't
%     one, it needs to be two.  And if it is two, then the channel exponent
%     should be one. If you want to explore other choices, you can by
%     commenting out the check.
%
% Inputs:
%    None.
%
% Outputs:
%    obj         - The object.
%
% Optional key/value pairs:
%    'dimension' - The dimension of the stimuli. Can be 2 or 3. Default
%                  3, but this only works if it is set to 2. (Not entirely
%                  clear how to build the channesl in a 3D angular space.)
%    'nChannels' - Number of cos^n angular tuning linear channels.  Scalar.
%                  Default 6.
%    'channelExponent'  - Exponent n in cos^n above. Scalar. Default 2.
%    'summationExponent' - Channel responses raised to this number prior to
%                  summation. Scalar. Default 1.
%    'startCenter' - Angle of center of first channel in degrees. Scalar.
%                  Default 0.
%    'criterionResp' - Criterion response for isoresponse controu. Scalar.
%                  Default 1.
%
% See also:
%

% History:
%   03/01/21 dhb       Wrote it.

    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
        % Number of unipolar channels
        nChannels = [];
        
        % Channel center start
        startCenter = [];
        
        % Angle support
        angleSupport = 1:1:360;
        
        % Mechanisms
        underlyingChannels;
        channelExponent;
        summationExponent;
        
        % Cache stimuli in desired form for fitting here
        angles = [];
        contrasts = [];
        
    end
    
    % Public, read-writ properties.
    properties (SetAccess = public, GetAccess = public)
        
        % Criterion response to go from response to isocontour      
        criterionResp = [];
        
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods  
    end
    
    properties (Dependent)
    end
    
    % Methods.  Most public methods are implemented in a separate function,
    % but we put the constructor here.
    methods (Access=public)
        % Constructor
        function obj = tfeLCMDirection(varargin)
            
            % Parse input. Need to add any key/value pairs that need to go
            % to the tfe parent class, as well as any that are LCM
            % specific.
            p = inputParser; p.KeepUnmatched = true;
            p.addParameter('nChannels',6,@(x) (isnumeric(x) & isscalar(x)));
            p.addParameter('channelExponent',2,@(x) (isnumeric(x) & isscalar(x)));
            p.addParameter('summationExponent',1,@(x) (isnumeric(x) & isscalar(x)));
            p.addParameter('startCenter',0,@(x) (isnumeric(x) & isscalar(x)));
            p.addParameter('criterionResp',1,@(x) (isnumeric(x) & isscalar(x)));
            p.parse(varargin{:});
                        
            % Base class constructor
            obj = obj@tfeQCM(varargin{:});
            
            % Set basic parameters
            obj.nChannels = p.Results.nChannels;
            obj.channelExponent = p.Results.channelExponent;
            obj.summationExponent = p.Results.summationExponent;
            obj.startCenter = p.Results.startCenter;
            obj.criterionResp = p.Results.criterionResp;
            
            % Don't allow parameter choices to get too crazy
            if (obj.summationExponent ~= 1)
                if (obj.summationExponent ~= 2)
                    error('Using summation exponent other than 1 or 2 is dangerous');
                end
                if (obj.channelExponent ~= 1)
                    error('Not sure you want to use channelExponent ~= 1 with summationExponent == 2');
                end
            end
            
            % Checks
            if (obj.dimension ~= 2)
                error('The LCM model only works in two dimensions');  
            end
            if (rem(obj.nChannels,2) ~= 0)
                error('nChannels must be even');
            end
            
            % Create the channels
            %
            % These have tuning described as half wave rectified sinusoids
            % squared. Compute response to unit contrast in each color
            % direction by regarding these as linear channels tuned to hue
            % angle.
            %
            % Set channel center points
            centerSpacing = 360/obj.nChannels;
            centerLocations = obj.startCenter:centerSpacing:360-centerSpacing+obj.startCenter;
            for ii = 1:obj.nChannels
                obj.underlyingChannels(ii,:) = cosd(obj.angleSupport-centerLocations(ii));
                obj.underlyingChannels(ii,:) = obj.underlyingChannels(ii,:).^(obj.channelExponent);
            end
        end
    end
    
    % Get methods for dependent properties
    methods
    end
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
    end
    
end
