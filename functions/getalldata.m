classdef getalldata
    % getalladata: Me auti tin klasi ginetai diavasma olon ton stoixeion
    % ton zonon pou iparxoun sta XML files pou einai ston fakelo
    % Benchmarks/House_8x23
    
    % ENTOLES
    % data = getalldata;
    % zonei = data.zoneidata;
    
    % paradeigma gia zone 1
    % zone1 = data.zone1data;
    
    properties
       zone1_data;
       zone2_data;
       zone3_data;
       zone4_data;
       zone5_data;
       zone6_data;
       zone7_data;
       zone8_data;
       zone9_data;
       zone10_data;
    end
    
    methods
        function obj = getalldata(varargin);
        end
        
        function groundkitchen = zone1data(zone1_data)
            zone1_data = getdata('Ground_kitchen.xml');
            groundkitchen = zone1_data;
        end
        
        function groundstorage = zone2data(zone2_data)
            zone2_data = getdata('Ground_storage.xml');
            groundstorage = zone2_data;
        end
    
        function floorbedroom2 = zone3data(zone3_data)
            zone3_data = getdata('Floor_bedroom2.xml');
            floorbedroom2 = zone3_data;
        end
        
        function floortoilet2 = zone4data(zone4_data)
            zone4_data = getdata('Floor_toilet2.xml');
            floortoilet2 = zone4_data;
        end
        
        function floorstairs = zone5data(zone5_data)
            zone5_data = getdata('Floor_stairs.xml');
            floorstairs = zone5_data;
        end
        
        function floorbedroom1 = zone6data(zone6_data)
            zone6_data = getdata('Floor_bedroom1.xml');
            floorbedroom1 = zone6_data;
        end
        
        function floortoilet1 = zone7data(zone7_data)
            zone7_data = getdata('Floor_toilet1.xml');
            floortoilet1 = zone7_data;
        end
        
        function sofitaroom1 = zone8data(zone8_data)
            zone8_data = getdata('Sofita_room1.xml');
            sofitaroom1 = zone8_data;
        end        
        
        function sofitaroom2 = zone9data(zone9_data)
            zone9_data = getdata('Sofita_room2.xml');
            sofitaroom2 = zone9_data;
        end
        
        function sofitastairs = zone10data(zone10_data)
            zone10_data = getdata('Sofita_stairs.xml');
            sofitastairs = zone10_data;
        end    
        
    end
    
    
    
end

