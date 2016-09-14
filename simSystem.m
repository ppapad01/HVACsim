clc;
clear all;
close all;
format long g
collect_data_building_v1;

Leng=length(str_z);

 for a=1:(Leng)
    xi(a)= str_z(a).Tz;
    %x(a)= str_z(a).Tz;   
 end


E(Leng+1)=0;
xi(Leng+1)=Twi;
xi(Leng+2:2*(Leng)+1)=0;
xi(2*Leng+3:3*(Leng)+3)=0;

% Xref=zeros(1,Leng+1);
% gain=zeros(1,Leng+1);
% itta=zeros(1,Leng+1);
% g=zeros(Leng+1,Leng+1);
% A=zeros(Leng+1,Leng+1);
% hh=zeros(1,Leng+1);

% for k=1:(Leng)
%     
%     Xref(k)= str_z(k).Tref;
%     gain(k)= 2; %str_z(k).k; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     g(k,k)=(((str_z(k).Umax * awz)) /(str_z(k).cz)) * (xi(Leng+1) - xi(k));
%     itta(k)=(str_z(k).az/str_z(k).cz)*Ta - ((h*str_z(k).Ad)/str_z(k).cz)*str_z(k).T1;
%     
%     sum4=0;
%     sum2=0;
%     sum3=0;
%     for c=1:(length(str_z(k).connected))     
%           sum4=sum4+str_z(k).az_ij(c)/str_z(k).cz; %sum for A
%           sum2=sum2+(str_z(k).az_ij(c)/str_z(k).cz)*xi(str_z(k).connected(c));
%           sum3 = sum3 + sign( xi(str_z(k).connected(c)) - xi(k) ) * str_z(k).paths.Ad_ij(1)...
%               * xi(k) * sqrt(2 * (Cp - Cv) * abs(xi(str_z(k).connected(c)) - xi(k)));
%     end
%      % str_z(k).Ad is the same this A_w_i
%     A(k,k)=(h*str_z(k).Ad-(str_z(k).az))/str_z(k).cz - sum4; 
%     hh(k) = sum2 + ((p_air * Cp)/str_z(k).cz) * sum3;
%     
% end
% 
% 
% Uz=(-(xi(1:Leng))*(A(1:Leng,1:Leng))-hh(1:Leng)-itta(1:Leng)+gain(1:Leng).*(Xref(1:Leng)-xi(1:Leng)))/(g(1:Leng,1:Leng));
% 
% sum1=0;
% for k=1:Leng
% sum1=sum1 + str_z(k).Umax * (xi(Leng+1)-xi(k))*Uz(k);
% end
% 
% A(Leng+1, Leng+1)= -aw/Cw; %aw
% g(Leng+1, Leng+1)= (Ustmax/Cw)* (1+(COPmax-1)*(1-((xi(Leng+1)-To))/DTmax));
% hh(Leng+1) = (awz/Cw)* sum1;%aw
% itta(Leng+1) = (aw/Cw) * Tpl;%aw
% kst=10;
% gain(Leng+1)=kst;
% Xref(Leng+1)=Trefw;
% 
% 
% Us=(1/g(Leng+1, Leng+1)) * (-A(Leng+1, Leng+1)*xi(Leng+1) - hh(Leng+1) - itta(Leng+1) + gain(Leng+1)*(Xref(Leng+1)-xi(Leng+1)));
% 
% U=[Uz Us];

%U=zeros(1,Leng+1);

% detection decision
dd(1:Leng+1)=0;

 for a=3*Leng+4:4*Leng+3
    xi(a)= str_z(a-(3*Leng+3)).Tz;
    %x(a)= str_z(a-(3*Leng+3)).Tz;   
 end
xi(4*Leng+4)=Twi;
%x(4*Leng+4)=Twi;
 

global mm U tt D E Ebar
mm=1; % counter
U={};
tt={};
D={};
E={};
Ebar={};
time=1;
options=odeset('RelTol',1e-3,'AbsTol',1e-2);%,'OutputFcn',@Controller)
%status = Controller(t,x, str_z, awz, h, p_air, Cp, kst, Trefw, Ta, Cv, aw, Cw, Ustmax, COPmax, DTmax, To, Tpl);
[t1,x1] = ode45(@(t,x) systemmv4(t, x, str_z, h, awz, Cw, Ustmax, COPmax, To,...
    DTmax, Cp, Cv, p_air, Ta, Tpl, aw, maxst, minst, ...
    kst, Trefw, L_s, n_bar_s, p_s, lambda_s, x_bar_s, r_bar_s, K, Leng, amp_out, amp_s, f_out, f_s), [0 time], xi, options);

%U1=Controller(t1,x1, str_z, awz, h, p_air, Cp, kst, Trefw, Ta, Cv, aw, Cw, Ustmax, COPmax, DTmax, To, Tpl);


figure
subplot(3,2,1)
plot(t1, x1(:,Leng+1));
subplot(3,2,2)
plot(t1, x1(:,1:Leng));
subplot(3,2,3)
plot(t1, x1(:,2*(Leng+1)));
subplot(3,2,4)
plot(t1, x1(:,Leng+2:2*(Leng)));
subplot(3,2,5)
plot(t1, x1(:,3*(Leng)+3));
subplot(3,2,6)
plot(t1, x1(:,2*Leng+3:3*(Leng)+2));


figure;hold all;
for k=1:Leng+1
    for j=1:mm-1
        Utmp(j) = U{j}(k);
        ttnew(j)=tt{j};
        Dtmp(j) = D{j}(k);
        Etmp(j) = E{j}(k);
        Ebartmp(j) = Ebar{j}(k);
    end
    Unew{k} = Utmp;
    Dnew{k} = Dtmp;
    Enew{k} = Etmp;
    Ebarnew{k} = Ebartmp;
    plot(ttnew,Enew{k},ttnew, Ebarnew{k})
end

figure
plot(ttnew,abs(Enew{1}),ttnew,Ebarnew{1})


figure
plot(ttnew,abs(Enew{84}),ttnew,Ebarnew{84})


