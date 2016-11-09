classdef Building
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        DATA
        areaUnit
        volumeUnit
        lengthUnit
        temperatureUnit
        %Type
        %Area
        %Volume
        %Location
        %Zones 
        %Name
    end
    
    methods
        function obj=Building(varargin)
        end
         
        function areaUnit=getareaunit(obj,DATA)
           areaUnit=getfield(DATA.Attributes(1).Value,{})
        end
        
    end
    
end

