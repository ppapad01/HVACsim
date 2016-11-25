classdef Spaces
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        get_id
        get_numberSpace
        get_areaSpace
        get_volumeSpace
        get_walls
    end
     
    methods
        function obj=Spaces(varargin,Building)
        end
     
     function get_numberSpace=numberSpace(obj,str)
         count=0;
         x=length(str.Children(2).Children(4).Children);
         s=char('Space');
        for i=1:x
               
                h=char(str.Children(2).Children(4).Children(i).Name);
                 if size(h)==size(s)
                    if h == s
                    count=count+1;
                    end
                 end
                
            end
            numberSpace=count
     end
     
    end
    
end

