classdef tfeQCM < tfe
% tfeQCM
%  Create at tfeQCM (Quadratic Color Model) object
%
% Syntax:
%    tfe = tfeQCM(varargin);
% 
% Description:
%     Implements a model that is quadratic in the color contrast of the
%     stimulus.
%
%     Inherits optional key/value pairs from parent class tfe.
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
%   06/26/16 dhb       Started in on this.
%   08/24/18 dhb, mab  Get key/value dimension pair working.

    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
        % Specify dimension of ellipsoid.  Can be 2 or 3.
        dimension = 3;
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
           
            % Parse input. Need to add any key/value pairs that need to go
            % to the tfe parent class, as well as any that are QCM
            % specific.
            p = inputParser;
            p.addParameter('verbosity','none',@ischar);
            p.addParameter('dimension',3,@(x) (isnumeric(x) & isscalar(x)));
            p.parse(varargin{:});
            
            % Base class constructor
            obj = obj@tfe('verbosity',p.Results.verbosity);
            
            % Set dimension
            obj.dimension = p.Results.dimension;
            if (obj.dimension ~= 2 & obj.dimension ~= 3)
                error('Can only handle dimension 2 or 3');
            end

        end
        
        % Override tfe fitError so we can scale appropriately for this
        % problem.  This makes a difference to fmincon's behavior, but
        % I did not want to change the tfe function itself. 
        function [fVal,modelResponseStruct] = fitError(obj,paramsVec,thePacket,varargin)

            [fVal,modelResponseStruct] = fitError@tfe(obj,paramsVec,thePacket,varargin{:});
            fVal= 1000*fVal;
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
