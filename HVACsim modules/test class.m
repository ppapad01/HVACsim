classdef hvac
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        m; %mikos;
        p; %platos;
        u; %upsos
    end
    
    methods
        function obj=hvac(varargin)
      
        end
        function v = getvolume(obj,m,p,u)
            
            v=m*p*u;
        end
    end
    
end

