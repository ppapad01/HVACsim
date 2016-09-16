% - - - - - - - - - Collect the data from the csv file - - - - - - - - -
% function [n_ind]=collect_data_building_v1()
format compact;
clc; format short;

tic

M = csvread('paths.csv');

level = 1;                           % enter the floor/level you want to work with
levels = M(:,11);                    % array w/ corresponding levels
l_ind = find(levels==level);         % finds all the indices of the "paths" which correspond to the given level

M_new = M(l_ind,:);
mult = M_new(:,15);

length(l_ind);
length(mult);

p_ind = find(mult==1);               % finds all the indices of the "paths" which correspond to a door,
                                     % WITHIN l_ind
path_indices = l_ind(p_ind);         % path indices
t_paths = length(p_ind);             % total # of doors - door_paths
clear t_paths;

paths = M_new(p_ind,3:4);

path_results = [ path_indices paths ];

tmp_indices = find(path_results(:,2)<=0);     % getting rid of the "exterior doors" --> correspond to negative zones
path_results(tmp_indices,:) = [];
tmp_indices = find(path_results(:,3)<=0);
path_results(tmp_indices,:) = [];
clear tmp_indices;

[t_doors ~] = size(path_results);         % total # of "interior" doors on given level

path_results;

Ad = 1.95096;       % Door Area  Adij
Height = 2.438;     % Height of ceiling
H = Height;
p_air = 1.225;     % Air Pressure
Cp=1.004;

% level_1 --> 393    level_2 --> 468    level_3 --> 79  PATHS
% level_1 --> 118    level_2 --> 118                    DOORS

%%  Reading zones

% z_data, z_type --> type string
% Mz_new --> type double

Mz = xlsread('zones.xls');                      % z for zones
[~, ~, zones_str] = xlsread('zones.xls',1);   % zone type/name

z_type = zones_str(:,11);
z_data = [zones_str(:,1) zones_str(:,6) zones_str(:,8) zones_str(:,11)];
%z_data = [cell2mat(zones_str(:,1)) cell2mat(zones_str(:,6)) cell2mat(zones_str(:,8)) zones_str(:,11)];


% ONLY FOR THE 1ST FLOOR!!!
Mz = xlsread('zones_1st_floor.xls');                      % z for zones
[~, ~, zones_str_1st] = xlsread('zones_1st_floor.xls',1);   % zone type/name
z_type_1st = zones_str_1st(:,11);
z_data_1st = [zones_str_1st(:,1) zones_str_1st(:,6) zones_str_1st(:,8) zones_str_1st(:,11)];



v1 = Mz(:,1);                   % zone #
v2 = Mz(:,6);                   % level #
l_ind_2 = find(v2==level);      % finds all the indices of the "paths" which correspond to the given level
v3 = Mz(:,8);                   % volume in cubic metres

Mz_new = [v1 v2 v3];            % [zone level volume]
Mz_new = Mz_new(l_ind_2,:);

z_type = z_type(l_ind_2,:);     % for only the slected level
z_data = z_data(l_ind_2,:);

[t_zones ~] = size(Mz_new);     % total # of zones in the given level
Mz_new(:,1) =  1:t_zones;       % "re-indexes" the zones of the given level, starting from 1

% old_indexes = 1:t_zones; old_indexes = old_indexes';
% new_indexes = -99*ones(t_zones,1);     % use "-99" for sanity check

%%  Rearranging the zones, and "assigning" the paths
clc

indices = func_ind;

paths_new = path_order(path_results,indices);

%% Assigning the volumes to the rearranged zones

[t_paths ~] = size(paths_new);  % total # of paths/doors
volumes = zeros(t_paths,2);
volumes(:,1) = 1:t_paths;

for i=1:t_paths
    o_ind = indices(i,2);       % "old" index
    n_ind = indices(i,1);       % "new" index.. basically n_ind:=i
    volumes(n_ind,2) = Mz_new(o_ind,3);
end

%% Defining parameters, taken from initial simulation

COPmax=3.5;
aw=12;
awz=0.6;
h=8.29;   % heat tranfer coefficient W/m^-2 C
Cw=8370;
Tpl=20;
Ta=5; %ambient
Cp=1.004;
Cv=.717;
R=(Cp-Cv);
Ustmax=27.36e5;
COPmax=3.5;
DTmax=45;
To=5;
Twi=30;
K=10000;
% - - - - - - - ??
% ust=1000;
% u1=0.6;
% u2=0.5;
% u3=.4;
% u4=.7;
% u5=.5;
% u6=0.6;
% u7=0.5;
% u8=.4;
% u9=.7;
% u10=.5;
time=24; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end time
kst=-12; 
Trefw=55;
% - - - - - - -
minst=-0.03*Trefw;
maxst=0.03*Trefw;
seedst=2;
stst=0.05;
% - - - - - - - ??
% Uncertainty
amp_s=0.10*Tpl;
f_s=2;
amp_out=0.10*Ta;
f_out=2;
% - - - - - - -

  % - - - Additional Parameters - - -
T_hat_s = 0;
L_s = 50;

F_time_s=20; %%%%%%%%%%%%%%%%%%%%%%% Time of fault occurence
F_value_s=0; %%%%%%%%%%%%%%%%%%%%%%% Magnitude of the fault

eps1 = 0.08;
x_bar_s = Twi+eps1;
r_bar_s=0.50*Tpl+eps;
n_bar_s=(maxst-minst)/2;
p_s=1.3;
x_s=30;
Omega_s=0;
xq_s=0;
theta_s=0;
gamma_s=8;
delta_s=(0.15*Trefw)+eps1;
x_barq_s=Twi+eps1;
lambda_s=25; %8;
kappa_s=1;


toc


%% Creating the structures

% struct-paths... struct-zones...
[str_p, str_z] = create_struct(t_zones, paths_new, volumes, Height, p_air, eps1, Ta)

str_z(1).F_value = 0.1*str_z(1).Tref ;
str_z(1).F_time = 0.5;
% sim Aprtmnts_floor1_v11
