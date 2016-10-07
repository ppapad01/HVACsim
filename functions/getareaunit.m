function obj=getareaunit(DATA)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    for i=1:length(DATA.Attributes)
        if DATA.Attributes.Name(i)=='areaUnit'
            obj=DATA.Attributes.Value(i)
        else
            error('Invalid input')
        end
    end
end

