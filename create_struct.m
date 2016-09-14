function [str_paths, str_zones] = create_struct(Nzones, paths, vol, height, p_air, eps1, Ta)

[t_paths ~] = size(paths);
[t_zones ~] = size(vol);

noise_perc = 0.03;

Tz = 20; %%%%%%%%%%%%%%%%%%% initial zone temperature

umax=3700;
az=740;



%% - - - - Creates the Path Structure - - - -

for i=1:t_paths
    str_paths(i).idx = paths(i,1);
    str_paths(i).zone1 = paths(i,2);        % 1st zone of the path
    str_paths(i).zone2 = paths(i,3);        % 2nd zone of the path
    str_paths(i).Ad = 1.95096;              % Path/Door area
end

%% - - - - Creates the Zone Structure - - - -

for i=1:80
    str_zones(i).idx = vol(i,1);
    str_zones(i).vol = vol(i,2);        % zone volume
    str_zones(i).height = 2.438;        % zone "height"
end

h = height;
area = vol(:,2)/h;      % array w/ zone areas
Cv=.717;

l1a = 2.7432;
l1b = 3.6576;

l2a = l1b;
l2b = 5.30352;
l2c = l2b-l1a;

l3a = 6.2179;
l3b = 1.2192;

l4a = 2.4384;
l4b = l1b;

l5a = l1b;
l5b = l4a;
l5c = l3b;
l2d = l2a-l5b;

l_stairs_1 = 2.4384;
l_stairs_2 = 4.5052;
l_stairs_3 = 2.4384;
    
for i=1:16               % 8 apartments to the right (top)
    
    f = 5*(i-1);            % f --> factor

    if ( (i==1)||(i==9) )       % "Corner" Apartment
        
      % - - - - - - Bedroom2 - - - - - -
        z = f+1;                         % "Overall" zone index
        str_zones(z).ConnWalls = 2;
        str_zones(z).connected = f+[2 3];
        str_zones(z).wAreas = h*[(l1a) (l1b)];
	  % - - - - - - Bedroom1 - - - - - -
        z = f+2;
        str_zones(z).ConnWalls = 4;
        str_zones(z).connected = f+[1 3 5 7];
        str_zones(z).wAreas = h*[(l1a) (l2c+l2d) (l5b) (l2b)];
	  % - - - - - - Living Room - - - - - -
        z = f+3;
        str_zones(z).ConnWalls = 5;
        str_zones(z).connected = f+[1 2 4 5 (83-f)];
        str_zones(z).wAreas = h*[(l1b) (l2c+l2d) (l4a+l4b) (l5c) (l3b)];
   	  % - - - - - - Kitchen - - - - - -
        z = f+4;
        str_zones(z).ConnWalls = 3;
        str_zones(z).connected = f+[3 5 (83-f)];
        str_zones(z).wAreas = h*[(l4a+l4b) (l4a) (l4b)];
   	  % - - - - - - Bathroom - - - - - -
        z = f+5;
        str_zones(z).ConnWalls = 5;
        str_zones(z).connected = f+[2 3 4 10 (83-f)];
        str_zones(z).wAreas = h*[(l5b) (l5c) (l4a) (l5a) (l5b)];
        
    elseif ( (i==8)||(i==16) )       % "Corner" Apartment
        
      % - - - - - - Bedroom2 - - - - - -
        z = f+1;                         % "Overall" zone index
        str_zones(z).ConnWalls = 2;
        str_zones(z).connected = f+[2 3];
        str_zones(z).wAreas = h*[(l1a) (l1b)];
	  % - - - - - - Bedroom1 - - - - - -
        z = f+2;
        str_zones(z).ConnWalls = 4;
        str_zones(z).connected = f+[(-3) 1 3 5];
        str_zones(z).wAreas = h*[(l2b) (l1a) (l2c+l2d) (l5b)];
	  % - - - - - - Living Room - - - - - -
        z = f+3;
        str_zones(z).ConnWalls = 5;
        str_zones(z).connected = f+[1 2 4 5 (83-f)];
        str_zones(z).wAreas = h*[(l1b) (l2c+l2d) (l4a+l4b) (l5c) (l3b)];
   	  % - - - - - - Kitchen - - - - - -
        z = f+4;
        str_zones(z).ConnWalls = 3;
        str_zones(z).connected = f+[3 5 (83-f)];
        str_zones(z).wAreas = h*[(l4a+l4b) (l4a) (l4b)];
   	  % - - - - - - Bathroom - - - - - -
        z = f+5;
        str_zones(z).ConnWalls = 5;
        str_zones(z).connected = f+[0 2 3 4 (83-f)];
        str_zones(z).wAreas = h*[(l5a) (l5b) (l5c) (l4a) (l5b)];
        
        
    elseif ( (i>1)&&(i<8) || (i>9)&&(i<16) )       % "Central" Apartment

     if (~mod(i,2))                  % "Even" numbered apartments
         
      % - - - - - - Bedroom2 - - - - - -
        z = f+1;                         % "Overall" zone index
        str_zones(z).ConnWalls = 3;
        str_zones(z).connected = f+[2 3 6];
        str_zones(z).wAreas = h*[(l1a) (l1b) (l1a)];
	  % - - - - - - Bedroom1 - - - - - -
        z = f+2;
        str_zones(z).ConnWalls = 4;
        str_zones(z).connected = f+[(-3) 1 3 5];
        str_zones(z).wAreas = h*[(l2b) (l1a) (l2c+l2d) (l5b)];
	  % - - - - - - Living Room - - - - - -
        z = f+3;
        str_zones(z).ConnWalls = 6;
        str_zones(z).connected = f+[1 2 4 5 8 (83-f)];
        str_zones(z).wAreas = h*[(l1b) (l2c+l2d) (l4a+l4b) (l5c) (l3a) (l3b)];
      % - - - - - - Kitchen - - - - - -
        z = f+4;
        str_zones(z).ConnWalls = 3;
        str_zones(z).connected = f+[3 5 (83-f)];
        str_zones(z).wAreas = h*[(l4a+l4b) (l4a) (l4b)];
   	  % - - - - - - Bathroom - - - - - -
        z = f+5;
        str_zones(z).ConnWalls = 5;
        str_zones(z).connected = f+[0 2 3 4 (83-f)];
        str_zones(z).wAreas = h*[(l5a) (l5b) (l5c) (l4a) (l5b)];
        
     elseif (mod(i,2))               % "Odd" numbered apartments
         
      % - - - - - - Bedroom2 - - - - - -
        z = f+1;                         % "Overall" zone index
        str_zones(z).ConnWalls = 3;
        str_zones(z).connected = f+[(-4) 2 3];
        str_zones(z).wAreas = h*[(l1a) (l1a) (l1b)];
	  % - - - - - - Bedroom1 - - - - - -
        z = f+2;
        str_zones(z).ConnWalls = 4;
        str_zones(z).connected = f+[1 3 5 7];
        str_zones(z).wAreas = h*[(l1a) (l2c+l2d) (l5b) (l2b)];
      % - - - - - - Living Room - - - - - -
        z = f+3;
        str_zones(z).ConnWalls = 6;
        str_zones(z).connected = f+[(-2) 1 2 4 5 (83-f)];
        str_zones(z).wAreas = h*[(l3a) (l1b) (l2c+l2d) (l4a+l4b) (l5c) (l3b)];
   	  % - - - - - - Kitchen - - - - - -
        z = f+4;
        str_zones(z).ConnWalls = 3;
        str_zones(z).connected = f+[3 5 (83-f)];
        str_zones(z).wAreas = h*[(l4a+l4b) (l4a) (l4b)];
   	  % - - - - - - Bathroom - - - - - -
        z = f+5;
        str_zones(z).ConnWalls = 5;
        str_zones(z).connected = f+[2 3 4 10 (83-f)];
        str_zones(z).wAreas = h*[(l5b) (l5c) (l4a) (l5a) (l5b)];

     end
     
   end     

end

clear z f;

    % - - - - - Entering the "surrounding-wall" area - - - - -
    % --> Perimeter
  m=100;
    
for i=1:Nzones-3

    if (mod(i,5)==3)
        str_zones(i).surf_perimeter = (l1b + l2c + l2d + l5c + l4b + l4a + l3b + l3a);
    elseif (mod(i,5)==0)
        str_zones(i).surf_perimeter = 2*(l5a+l5b);
    elseif (mod(i,5)==1)
        str_zones(i).surf_perimeter = 2*(l1a+l1b);
    elseif (mod(i,5)==2)
        str_zones(i).surf_perimeter = 2*(l2a+l2b);
    elseif (mod(i,5)==4)
        str_zones(i).surf_perimeter = 2*(l4a+l4b);
    else
        fprintf('\n- - - - - ERROR In "surf_perimeter" - - - - -\n')
    end
    
    str_zones(i).cz = str_zones(i).vol*p_air*Cv*m;
  % MIGHT HAVE TO CHANGE!!!
    str_zones(i).Ad = h*str_zones(i).surf_perimeter; %total area of the walls in zone i
    %str_zones(i).Umax = (str_zones(i).vol/min_vol)*umax; %3700
    str_zones(i).az = az;
    str_zones(i).Afloor = str_zones(i).vol/str_zones(i).height;
    str_zones(i).az_ij = 50*(ones(size(str_zones(i).connected)));
    %str_zones(i).Ad_ij = 1.95096*(ones(size(str_zones(i).connected)));
end


    % - - - - - - Paths - SubStructure - - - - - -

for i=1:Nzones-3
    
   % fprintf('\n- - - - - - - - - - -\n')
   tmp_ind_1 = find(paths(:,2)==i); p1 = paths(tmp_ind_1,3);
   tmp_ind_2 = find(paths(:,3)==i); p2 = paths(tmp_ind_2,2);
   str_zones(i).paths.ConnDoors = sort([p1' p2']);
   [tmp num_paths] = size(str_zones(i).paths.ConnDoors);
   str_zones(i).paths.TotalPaths = num_paths;
   str_zones(i).paths.Ad_ij = 1.95096*ones(1,num_paths);
   
   str_zones(i).Tz = Tz;
   str_zones(i).T1 = 10;        % nodes

end

clear tmp_ind_1 tmp_ind_2 tmp num_paths i

% - - - - - Add "stair-connections" to zones 3 and 38 - - - - -

str_zones(3).connected = [1 2 4 5 81 83];
str_zones(3).ConnWalls = length(str_zones(3).connected);
str_zones(3).az_ij = 50*(ones(size(str_zones(3).connected)));
str_zones(3).wAreas = h*[(l1b) (l2c+l2d) (l4a+l4b) (l5c) (l_stairs_2) (l3b)];

str_zones(38).connected = [36 37 39 40 82 83];
str_zones(38).ConnWalls = length(str_zones(38).connected);
str_zones(38).az_ij = 50*(ones(size(str_zones(38).connected)));
str_zones(38).wAreas = h*[(l1b) (l2c+l2d) (l4a+l4b) (l5c) (l_stairs_2) (l3b)];

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    % - - - - - Create struct for zones 81, 82, 83 - - - - - 

for z=Nzones-2:Nzones-1

    str_zones(z).idx = z;
    str_zones(z).vol = vol(z,2);        % zone volume
    str_zones(z).height = 2.438;        % zone "height"
    str_zones(z).ConnWalls = 2;

    if (z==81)
        str_zones(z).connected = [3 83];
    elseif (z==82)
        str_zones(z).connected = [38 83]; 
    end

    str_zones(z).wAreas = h*[l_stairs_2 l_stairs_1];
    str_zones(z).surf_perimeter = 2*(l_stairs_1 + l_stairs_2 + l_stairs_3);
    str_zones(z).cz = str_zones(z).vol*p_air*Cv*m;
    str_zones(z).Ad = h*str_zones(z).surf_perimeter; %total area of the walls in zone i
    %str_zones(z).Umax = umax;
    str_zones(z).az = az;
    str_zones(z).Afloor = str_zones(z).vol/str_zones(z).height;
    str_zones(z).az_ij = 50*(ones(size(str_zones(z).connected)));
    str_zones(z).paths.ConnDoors = 83;
    str_zones(z).paths.TotalPaths = 1;
    str_zones(z).paths.Ad_ij = 1.95096;
    str_zones(z).Tz = Tz;
    str_zones(z).T1 = 10;        % nodes

end




for z=Nzones
    str_zones(z).idx = z;
    str_zones(z).vol = vol(z,2);        % zone volume
    str_zones(z).height = 2.438;        % zone "height"

  % - - - - - - Find connected zones - - - - - -
  tmp_conn = [];
  for i = 1:Nzones-3
      if ( mod(i,5)==3 || mod(i,5)==4 || mod(i,5)==0 )
          tmp_conn = [tmp_conn i];
      end
  end
  tmp_conn = [tmp_conn 81 82];
  str_zones(z).connected = tmp_conn;
  str_zones(z).ConnWalls = length(tmp_conn);
  clear tmp_conn i
  % - - - - - - - - - - - - - - - - - - - - - - -

  % - - - - - - - - Find wall Areas - - - - - - - - 
  tmp_wPer = [l3b l4a l4b];                                   % 1st apartment
  tmp_wPer = [tmp_wPer tmp_wPer tmp_wPer tmp_wPer];   % first 4 apartments
  tmp_wPer = [tmp_wPer tmp_wPer];                         % first 8 apartments
  tmp_wPer = [tmp_wPer tmp_wPer];                         % all 16 apartments
  tmp_wPer = [tmp_wPer l_stairs_1 l_stairs_1];              % stair connected walls

  str_zones(z).wAreas = h*tmp_wPer;
  str_zones(z).surf_perimeter = sum(tmp_wPer);
  clear tmp_wPer
  % - - - - - - - - - - - - - - - - - - - - - - -

    str_zones(z).cz = str_zones(z).vol*p_air*Cv*m;
    str_zones(z).Ad = h*str_zones(z).surf_perimeter; %total area of the walls in zone i
    %str_zones(z).Umax = umax;
    str_zones(z).az = az;
    str_zones(z).Afloor = str_zones(z).vol/str_zones(z).height;
    str_zones(z).az_ij = 50*(ones(size(str_zones(z).connected)));

  % - - - - - - - - Find connected Doors - - - - - - - -
  tmp_cDoors = 0:15;
  tmp_cDoors = 3 + tmp_cDoors*5;
  tmp_cDoors = [tmp_cDoors 81 82];
  
  str_zones(z).paths.ConnDoors = tmp_cDoors;
  str_zones(z).paths.TotalPaths = length(tmp_cDoors);
  str_zones(z).paths.Ad_ij = 1.95096*ones(1,length(tmp_cDoors));
  clear tmp_cDoors;  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
    str_zones(z).Tz = Tz;
    str_zones(z).T1 = 10;        % nodes

end
clear z

for i=1:t_zones
    
    str_zones(i).Tref = 24;     % Referenced Temperature
    str_zones(i).k = -25;        % Controller gain
    
        % - - - Bounds - - -
    str_zones(i).min = -noise_perc*str_zones(i).Tref;
    str_zones(i).max = +noise_perc*str_zones(i).Tref;
    
    str_zones(i).seed = round(rand*1000);   % Generating a random "seed"
    str_zones(i).st = 0.05;
    
        % - - - Additional Parameters - - -
    str_zones(i).T_hat = 0;
    str_zones(i).L = 20;
    str_zones(i).F_time = 10;
    str_zones(i).F_value = 0;
    str_zones(i).x_bar = str_zones(i).Tz + eps1;
    str_zones(i).r_bar = 0.50 * Ta + eps;
    str_zones(i).n_bar = (str_zones(i).max-str_zones(i).min)/2;
    str_zones(i).y_bar = str_zones(i).Tref + str_zones(i).max;
    str_zones(i).he_bar = eps;
    str_zones(i).p = 3/2;
    str_zones(i).xi = 6;
    str_zones(i).Omega = 0;
    str_zones(i).xq = 0;
    str_zones(i).theta = 0;
    str_zones(i).gamma = 5;
    str_zones(i).delta = (0.15*str_zones(i).Tref) + eps1;
    str_zones(i).x_barq = str_zones(i).Tz +eps1;
    str_zones(i).lambda = 10; %4;
    str_zones(i).kappa = 1;
    
    
        str_zones(i).xis = 12;
    str_zones(i).ps = 1;
    str_zones(i).Ls = 25;
    str_zones(i).x_bars = str_zones(i).Tz +eps1;
    str_zones(i).gammas = 20;
    str_zones(i).dfs_bar=0.001;
    str_zones(i).fs_bar=0.001;
    str_zones(1).f_bar=0.001;
    str_zones(1).thetas=0;
    
    str_zones(i).xia = 12;
    str_zones(i).pa = 1;
    str_zones(i).La = 25;
    str_zones(i).x_bara = str_zones(i).Tz +eps1;
    str_zones(i).gammaa = 500;
    str_zones(i).dfa_bar=eps;
    str_zones(i).fa_bar=eps;
    str_zones(1).f_bar=0.001;
    str_zones(1).thetaa=0;
        
end

min_vol=1e15;
for j=1:Nzones
    if str_zones(j).vol < min_vol
        min_vol=str_zones(j).vol;
    end
    str_zones(j).Umax = (str_zones(j).vol/min_vol)*umax; %3700
end


clear i

%% Correct zone-54.. As we get one additional door/path

str_zones(54).paths = str_zones(64).paths;
str_zones(54).paths.ConnDoors = str_zones(64).paths.ConnDoors - 10;

%%

% - - - - - - - - - What was used for the wall-areas approximations - - - - - - - - -
%{
for i=1:16               % 8 apartments to the right (top)
    
    f = 5*(i-1);            % f --> factor
   
  % - - Estimating wall lengths - -
	  % Theses are only rough estimations!!!

    l1 = sqrt(area(f+1));
    l2 = sqrt(area(f+2)/2);
    l4 = sqrt(area(f+4));
	l5 = sqrt(area(f+5)/3);
    l3 = (2*l2+3*l5)-l1;clc
    
    l6 = l2-l5;
    l7 = 2*l2-l1;
    l8 = 3*l5-l4;
    l9 = (l1+l2)-(l4+l5);
    
    
   if ( (i==1)||(i==9) )       % "Corner" Apartment
        
      % - - - - - - Bedroom2 - - - - - -
        z = f+1;                         % "Overall" zone index
        str_zones(z).ConnWalls = 2;
        str_zones(z).connected = f+[2 3];
        str_zones(z).wAreas = h*[(l1) (l1)];
	  % - - - - - - Bedroom1 - - - - - -
        z = f+2;
        str_zones(z).ConnWalls = 4;
        str_zones(z).connected = f+[1 3 5 7];
        str_zones(z).wAreas = h*[(l1) (l6+l7) (l5) (2*l2)];
	  % - - - - - - Living Room - - - - - -
        z = f+3;
        str_zones(z).ConnWalls = 5;
        str_zones(z).connected = f+[1 2 4 5 (83-f)];
        str_zones(z).wAreas = h*[(l1) (l6+l7) (2*l4) (l8) (l9)];
   	  % - - - - - - Kitchen - - - - - -
        z = f+4;
        str_zones(z).ConnWalls = 3;
        str_zones(z).connected = f+[3 5 (83-f)];
        str_zones(z).wAreas = h*[(2*l4) (l4) (l4)];
   	  % - - - - - - Bathroom - - - - - -
        z = f+5;
        str_zones(z).ConnWalls = 5;
        str_zones(z).connected = f+[2 3 4 10 (83-f)];
        str_zones(z).wAreas = h*[(l5) (l8) (l4) (3*l5) (l5)];
        
    elseif ( (i==8)||(i==16) )       % "Corner" Apartment
        
      % - - - - - - Bedroom2 - - - - - -
        z = f+1;                         % "Overall" zone index
        str_zones(z).ConnWalls = 2;
        str_zones(z).connected = f+[2 3];
        str_zones(z).wAreas = h*[(l1) (l1)];
	  % - - - - - - Bedroom1 - - - - - -
        z = f+2;
        str_zones(z).ConnWalls = 4;
        str_zones(z).connected = f+[(-3) 1 3 5];
        str_zones(z).wAreas = h*[(2*l2) (l1) (l6+l7) (l5)];
	  % - - - - - - Living Room - - - - - -
        z = f+3;
        str_zones(z).ConnWalls = 5;
        str_zones(z).connected = f+[1 2 4 5 (83-f)];
        str_zones(z).wAreas = h*[(l1) (l6+l7) (2*l4) (l8) (l9)];
   	  % - - - - - - Kitchen - - - - - -
        z = f+4;
        str_zones(z).ConnWalls = 3;
        str_zones(z).connected = f+[3 5 (83-f)];
        str_zones(z).wAreas = h*[(2*l4) (l4) (l4)];
   	  % - - - - - - Bathroom - - - - - -
        z = f+5;
        str_zones(z).ConnWalls = 5;
        str_zones(z).connected = f+[2 3 4 (-5) (83-f)];
        str_zones(z).wAreas = h*[(l5) (l8) (l4) (3*l5) (l5)];
        
        

        % I THINK THERE WERE A FEW ISSUES FROM HERE ON! (mixed up even/odd cases)
    elseif ( (i>1)&&(i<8) )       % "Central" Apartment

     if (~mod(i,2))                  % "Even" numbered apartments
         
      % - - - - - - Bedroom2 - - - - - -
        z = f+1;                         % "Overall" zone index
        str_zones(z).ConnWalls = 3;
        str_zones(z).connected = f+[(-4) 2 3];
        str_zones(z).wAreas = h*[(l1) (l1) (l1)];
	  % - - - - - - Bedroom1 - - - - - -
        z = f+2;
        str_zones(z).ConnWalls = 4;
        str_zones(z).connected = f+[1 3 5 7];
        str_zones(z).wAreas = h*[(l1) (l6+l7) (l5) (2*l2)];
	  % - - - - - - Living Room - - - - - -
        z = f+3;
        str_zones(z).ConnWalls = 5;
        str_zones(z).connected = f+[1 2 4 5 8 (83-f)];
        str_zones(z).wAreas = h*[(l1) (l6+l7) (2*l4) (l8) (l3) (l9)];
      % - - - - - - Kitchen - - - - - -
        z = f+4;
        str_zones(z).ConnWalls = 3;
        str_zones(z).connected = f+[3 5 (83-f)];
        str_zones(z).wAreas = h*[(2*l4) (l4) (l4)];
   	  % - - - - - - Bathroom - - - - - -
        z = f+5;
        str_zones(z).ConnWalls = 5;
        str_zones(z).connected = f+[0 2 3 4 (83-f)];
        str_zones(z).wAreas = h*[(3*l5) (l5) (l8) (l4) (l5)];
        
     elseif (mod(i,2))               % "Odd" numbered apartments
         
      % - - - - - - Bedroom2 - - - - - -
        z = f+1;                         % "Overall" zone index
        str_zones(z).ConnWalls = 3;
        str_zones(z).connected = f+[2 3 6];
        str_zones(z).wAreas = h*[(l1) (l1) (l1)];
	  % - - - - - - Bedroom1 - - - - - -
        z = f+2;
        str_zones(z).ConnWalls = 4;
        str_zones(z).connected = f+[(-3) 1 3 5];
        str_zones(z).wAreas = h*[(2*l2) (l1) (l6+l7) (l5)];
      % - - - - - - Living Room - - - - - -
        z = f+3;
        str_zones(z).ConnWalls = 5;
        str_zones(z).connected = f+[(-2) 1 2 4 5 (83-f)];
        str_zones(z).wAreas = h*[(l3) (l1) (l6+l7) (2*l4) (l8) (l9)];
   	  % - - - - - - Kitchen - - - - - -
        z = f+4;
        str_zones(z).ConnWalls = 3;
        str_zones(z).connected = f+[3 5 (83-f)];
        str_zones(z).wAreas = h*[(2*l4) (l4) (l4)];
   	  % - - - - - - Bathroom - - - - - -
        z = f+5;
        str_zones(z).ConnWalls = 5;
        str_zones(z).connected = f+[2 3 4 10 (83-f)];
        str_zones(z).wAreas = h*[(l5) (l8) (l4) (3*l5) (l5)];

     end
     
   end     

end
%}
