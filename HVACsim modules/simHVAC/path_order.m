function paths_n = path_order(paths_o,indices)      % o --> old     n --> new

paths_n = zeros(size(paths_o));

[t_paths ~] = size(paths_o);            % total # of paths/doors

paths_n(:,1) = 1:t_paths;

paths_o(:,2:3) = paths_o(:,2:3)-83;


paths_o;
ff1 = paths_o(:,2:3);
max(ff1);
s1 = size(ff1);

% Getting rid of paths involving zones 167-170
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
tmp_indices = find(paths_o(:,2)>=84);     % getting rid of the "exterior doors" --> correspond to negative zones
tmp1 = paths_o(tmp_indices,1);
paths_o(tmp_indices,:) = [];

tmp_indices = find(paths_o(:,3)>=84);
tmp2 = paths_o(tmp_indices,1);
paths_o(tmp_indices,:) = [];

clear tmp_indices;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

[t_paths_inside ~] = size(paths_o);         % Getting rid of the paths which are considered to be outside
paths_n = paths_n(1:t_paths_inside,:);

paths_o;
ff2 = paths_o(:,2:3);
max(ff2);
s2 = size(ff2);

p2 = paths_o(:,2);

p3 = paths_o(:,3);


%% At this point, I have all the necessary paths, w/ the "old" indices

for i=1:length(indices)
     % Replace new path-zones in the 1st row
% - - - - - - - - - - - - - - - - - - - - - - -
    n_ind = indices(i,1);       % "new" index
    o_ind = indices(i,2);       % "old" index
    
    tmp_indices = find(p2==o_ind);
    paths_n(tmp_indices,2) = n_ind;
    
    clear tmp_indices;
% - - - - - - - - - - - - - - - - - - - - - - -

     % Replace new path-zones in the 2nd row
% - - - - - - - - - - - - - - - - - - - - - - -
    tmp_indices = find(p3==o_ind);
    paths_n(tmp_indices,3) = n_ind;
    
    clear tmp_indices;
% - - - - - - - - - - - - - - - - - - - - - - -
end

%% Rearrange 'path_n", s.t. the smaller-indexed zone is in the left column 

for i=1:t_paths_inside
    entr_l = paths_n(i,2);      % left entry
    entr_r = paths_n(i,3);      % right entry
    if (entr_r<entr_l)          % swaps the entries
        paths_n(i,2) = entr_r;
        paths_n(i,3) = entr_l;
    end
end

paths_n = sortrows(paths_n,2);

paths_n(:,1) = 1:t_paths_inside;





