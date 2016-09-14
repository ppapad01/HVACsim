function indices = indices_new2(indices);

%indices = 1:80;
indices = [indices indices];

[t_zones ~] = size(indices);

for i=1:t_zones;

    if (mod(indices(i,2),10)==6)
        indices(i,2) = indices(i,2)+4;      % XX6 --> XY0
    elseif(mod(indices(i,2),10)==7)
        indices(i,2) = indices(i,2)+2;      % XX7 --> XX9
    elseif (mod(indices(i,2),10)==9)
        indices(i,2) = indices(i,2)-2;      % XX9 --> XX7
    elseif(~mod(indices(i,2),10))
        indices(i,2) = indices(i,2)-4;      % XY0 --> XX6
    end

end