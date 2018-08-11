classdef tfeQCM < tfe
% tfeQCM
%
%   tfe = tfeQCM(varargin);
% 
% Implements a model that is quadratic in the color contrast of the
% stimulus.
%
% Inherits optional key/value pairs from parent class tfe.
%
% 6/26/16  dhb  Started in on this.

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
        function obj = tfeQCM(varargin)
            
            % Base class constructor
            obj = obj@tfe(varargin{:});
            
            % You could do model specific stuff here.

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
