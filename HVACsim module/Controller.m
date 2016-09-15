function  U  = Controller(t, y, str_z, awz, h, p_air, Cp,...
    kst, Trefw, Ta, Cv, aw, Cw, Ustmax, COPmax, DTmax, To, Tpl, Leng)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% collect_data_building_v1
% 
%  for a=1:(Leng)
%     xi(a)= str_z(a).Tz;
%     x(a)= str_z(a).Tz;
%     
%  end
% 
%  for b=1:(Leng+1)   %temporary U vector
%     U(b)=0;
%     A(b,b)=0;
%     g(b,b)=0;
%     hh(b)=0;
%     itta(b)=0;
% 
%  end

% yi(Leng+1)=Twi;
% y(Leng+1)=Twi;


sum1=0;
for k=1:(Leng)
    
    Xref(k)= str_z(k).Tref;
    gain(k)= 5; %str_z(k).k; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g(k,k)=(((str_z(k).Umax * awz)) /(str_z(k).cz)) * (y(Leng+1) - y(k));
    itta(k)=(str_z(k).az/str_z(k).cz)*Ta - ((h*str_z(k).Ad)/str_z(k).cz)*str_z(k).T1;
    
    sum4=0;
    sum2=0;
    sum3=0;
    for c=1:(length(str_z(k).connected))     
          sum4=sum4+str_z(k).az_ij(c)/str_z(k).cz; %sum for A
          sum2=sum2+(str_z(k).az_ij(c)/str_z(k).cz)*y(str_z(k).connected(c));
          sum3 = sum3 + sign( y(str_z(k).connected(c)) - y(k) ) * str_z(k).paths.Ad_ij(1)...
              * y(k) * sqrt(2 * (Cp - Cv) * abs(y(str_z(k).connected(c)) - y(k)));
    end
     % str_z(k).Ad is the same this A_w_i
    A(k,k)=(h*str_z(k).Ad-(str_z(k).az))/str_z(k).cz - sum4; 
    hh(k) = sum2 + ((p_air * Cp)/str_z(k).cz) * sum3;
   
    %sum1=sum1 + str_z(k).Umax * (y(Leng+1)-y(k))*U(k);
    
end


% Uz=((ones(1,83)/(g(1:Leng,1:Leng)))) .* (-(y(1:Leng))*A(1:Leng,1:Leng)-hh(1:Leng)-itta(1:Leng)+gain.*(Xref-y(1:Leng)));
% Uz=(-(A(1:Leng,1:Leng)*((y(1:Leng))'))'-hh(1:Leng)-itta(1:Leng)+gain.*(Xref-y(1:Leng)))/(g(1:Leng,1:Leng));
Uz=(-(y(1:Leng))*(A(1:Leng,1:Leng))-hh(1:Leng)-itta(1:Leng)+gain.*(Xref-y(1:Leng)))/(g(1:Leng,1:Leng));


 sum1=0;
 for k=1:Leng
 sum1=sum1 + str_z(k).Umax * (y(Leng+1)-y(k))*Uz(k);
 end

A(Leng+1, Leng+1)= -aw/Cw; %aw
g(Leng+1, Leng+1)= (Ustmax/Cw)* ...
    (1+(COPmax-1)*(1-((y(Leng+1)-To))/DTmax));
hh(Leng+1) = (awz/Cw)* sum1;%aw
itta(Leng+1) = (aw/Cw) * Tpl;%aw
kst=10;
gain(Leng+1)=kst;
Xref(Leng+1)=Trefw;


%Us=(1/g(Leng+1, Leng+1)) .*...
%(-A(Leng+1, Leng+1)*y(Leng+1) - hh(Leng+1) - itta(Leng+1) + kst*(Trefw-y(Leng+1)));

%U=[Uz Us];
U=(-y(1:Leng+1)*A-hh-itta+gain.*(Xref-y(1:Leng+1)))/g;

% % U(Leng+1)=U1(Leng+1);
% for j=1:Leng+1
%     if U1(j)>1
%     U(j)=1;
%     elseif U1(j)<0 
%     U(j)=0.001;
%     else 
%     U(j)=U1(j);
%     end
% end



% Ac=A;
% gc=g;
% hhc=hh;
% ittac=itta;
% state_c=y(1:Leng+1);
% 
% save controller_Data.mat Ac gc hhc ittac state_c
end

