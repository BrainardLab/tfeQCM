classdef tfeNakaRushtonDirection < tfeQCM
% tfeNakaRushtonDirection
%  Create at tfeNakaRushtonDirection object
%
% Syntax:
%    tfe = tfeNakaRushtonDirection(varargin);
% 
% Description:
%     Implements a model that is quadratic in the color contrast of the
%     stimulus. Stimulus specified as unit vector direction and contrast.
%
%     Inherits optional key/value pairs from parent class tfeQCM, plus
%     those specified below.
%
% Inputs:
%    None.
%
% Outputs:
%    obj         - The object.
%
% Optional key/value pairs:
%      'nDirections'           - Scalar (default 1). Number of color direcitons that will be fit.
%      'lockOffsetToZero'      - Logical (default false). Force fits to go through 0 at 0 contrast
%      'commonAmplitude'       - Logical (default false). Force common amplitude across directions.
%      'commonSemi'            - Logical (default false). Force common semi-saturation across directions.
%      'commonExp'             - Logical (default false). Force common exponent across directions.
%      'commonOffset'          - Logical (default false). Force common offset across directions.
%
% See also:
%

% History:
%   12/09/18 dhb       Wrote it.

    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
        % Fitting constraints
        lockOffsetToZero = false;
        commonAmplitude = false;
        commonSemi = false; 
        commonExp = false;
        commonOffset = false;
        
        % Number of color directions.  We need to know this to set up
        % parameters.
        nDirections = 1;
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
        function obj = tfeNakaRushtonDirection(varargin)
           
            % Parse input. Need to add any key/value pairs that need to go
            % to the tfe parent class, as well as any that are QCM
            % specific.
            p = inputParser; p.KeepUnmatched = true;
            p.addParameter('nDirections',1,@isnumeric);
            p.addParameter('lockOffsetToZero',false,@islogical);
            p.addParameter('commonAmplitude',false,@islogical);
            p.addParameter('commonSemi',false,@islogical);
            p.addParameter('commonExp',false,@islogical);
            p.addParameter('commonOffset',false,@islogical);
            p.parse(varargin{:});
            
            % Base class constructor
            obj = obj@tfeQCM(varargin{:});
            
            % Set properties for this class
            nDirections = p.Results.nDirections;
            lockOffsetToZero = p.Results.lockOffsetToZero;
            commonAmplitude = p.Results.commonAmplitude;
            commonSemi = p.Results.commonSemi;
            commonExp = p.Results.commonExp;
            commonOffset = p.Results.commonOffset;
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
