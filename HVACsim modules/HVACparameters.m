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
Ad = 1.95096;       % Door Area  Adij
Height = 2.438;     % Height of ceiling
H = Height;
p_air = 1.225;     % Air Pressure
Cp=1.004;
noise_perc = 0.03;

Tz = 20; %%%%%%%%%%%%%%%%%%% initial zone temperature

umax=3700;
az=740;
time=24; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end time
kst=-12; 
Trefw=55;
% - - - - - - -
minst=-0.03*Trefw;
maxst=0.03*Trefw;
seedst=2;
stst=0.05;


% Uncertainty
amp_s=0.10*Tpl;
f_s=2;
amp_out=0.10*Ta;
f_out=2;



  % - - - Additional Parameters - - -
T_hat_s = 0;
L_s = 50;


K=10000;
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



for i=1:Nzones
    str_z(i)=id;
    str_z(i).width;
    str_z(i).height;
    str_z(i).str_z(i).height*str_z(i).width;         % volume of each zone m^3
    str_z(i).cz = str_z(i).vol*p_air;
    str_z(i).h= sum_of_hi;
    str_z(i).Ad=sum                               % sum_of_area_of_all_walls_inzone_i
    str_z(i).az = az;
    str_z(i).Afloor = str_z(i).vol/str_zones(i).height;
    
    str_z(i).connected = [1 2 4 5 81 83];  % ? i.e., zone i is wall connected with zones 1 2 4 5 81 83 

    str_z(i).az_ij = 50*(ones(size(str_zones(i).connected)));
    str_z(i).paths.Ad_ij =                 % ? area of door connecting zone i with j
    str_z(i).Tz = Tz;
    str_z(i).T1 = 10;
   
end