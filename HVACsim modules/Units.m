classdef Units
    %Units < Building
    %   Units is a subclass of Building.
    %   Units contains the units that we use, like area Unit,
    %   length Unit, temperature Unit and volume Unit. 
    
    properties
        getareaUnit;
        getlengthUnit;
        gettemperatureUnit;
        getvolumeUnit;
    end
    
    methods
        function obj=Units(varargin,Building)
        end
        
        function getareaUnit=areaUnit(obj,str)
            areaUnit=str.Attributes(1).Value
        end
        
        function getlengthUnit=lengthUnit(obj,str)
            lengthUnit=str.Attributes(2).Value
        end
       
        function gettemperatureUnit=temperatureUnit(obj,str)
            temperatureUnit=str.Attributes(3).Value
        end
        
        function getvolumeUnit=volumeUnit(obj,str)
            volumeUnit=str.Attributes(6).Value
        end
        
    end
    
end

