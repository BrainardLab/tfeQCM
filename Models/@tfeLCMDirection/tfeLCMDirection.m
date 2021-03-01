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
%
% See also:
%

% History:
%   03/01/21 dhb       Wrote it.

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
            % to the tfe parent class, as well as any that are LCM
            % specific.
            p = inputParser; p.KeepUnmatched = true;
            p.addParameter('dimension',3,@(x) (isnumeric(x) & isscalar(x)));
            p.parse(varargin{:});
            
            % Check dimension
            if (p.Results.dimension ~= 2)
                error('The LCM model only works in two dimensions');
            end
            
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
