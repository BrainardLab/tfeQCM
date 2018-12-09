classdef tfeQCMDirection < tfeQCM
% tfeQCMDirection
%  Create at tfeQCMDirection (Quadratic Color Model) object
%
% Syntax:
%    tfe = tfeQCMDirection(varargin);
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
%    'dimension' - The dimension of the ellipsoid. Can be 2 or 3. Default
%                  3.
% See also:
%

% History:
%   12/09/18 dhb       Wrote it.

    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
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
        function obj = tfeQCMDirection(varargin)
           
            % Parse input. Need to add any key/value pairs that need to go
            % to the tfe parent class, as well as any that are QCM
            % specific.
            p = inputParser; p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            % Base class constructor
            obj = obj@tfeQCM(varargin{:});
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
