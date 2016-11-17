classdef Building
    %   Building is the main class.
    %   Building contains the type, the address and the area of the
    %   building, the units that we use and the space that it
    %   contains details about the zones of the building. 
    %   Units and space are other classes that they may contain
    %   other classes, fields, structs and functions.
    
    properties
        getType;
        getAddress;
        getArea;
        getUnits=struct();
        getSpace;
    end
    
    methods
        function obj=Building(varargin)
        end
        
        function getType=Type(obj,str)
            Type=str.Children(2).Children(4).Attributes.Value
        end
        
        function getAddress=Address(obj,str)
            Address=str.Children(2).Children(4).Children(2).Children.Data
        end
        
        function getArea=Area(obj,str)
            Area=str.Children(2).Children(4).Children(4).Children.Data
        end
        
        function getUnits=Units(obj,str)
            
        end
        
    end
    
end

