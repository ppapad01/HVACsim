
function dx= systemmv4(t, x, str_z, h, awz, Cw, Ustmax,...
                        COPmax, To, DTmax, Cp,...
                        Cv, p_air, Ta, Tpl, aw, maxst,... 
                        minst, kst, Trefw, L_s, n_bar_s,...
                        p_s, lambda_s, x_bar_s, r_bar_s,...
                        K, Leng, amp_out, amp_s, f_out, f_s)

     n=zeros(1,Leng+1);
     y=zeros(1,Leng+1);
     f=zeros(1,Leng+1);
if t==0
 y(1:Leng+1)=x(3*Leng+4:4*Leng+4);
else
    for k=1:Leng
        if t<str_z(k).F_time
            f(k)=0;
        else
            f(k)=str_z(k).F_value*(1-exp(-K*t)); 
        end
    n(k)=str_z(k).max+2*str_z(k).min*rand(1); %%%%%% noise at sensor of zones
    y(k) =x(k)+n(k)+f(k);
    end
     
end
     %noise for the storage tank
     n(Leng+1)=maxst+2*minst*rand(1);
     y(Leng+1) = x(Leng+1) + n(Leng+1);             
     

     
global mm U tt D E Ebar
tt{mm}=t;

%     if tt{mm}==tt{mm-1}
%         mm=mm;
%         
        
U{mm}  = Controller(t, y, str_z, awz, h, p_air, Cp, kst, Trefw, Ta, Cv, aw, Cw, Ustmax, COPmax, DTmax, To, Tpl, Leng);




% g=zeros(Leng+1,Leng+1);
% L=zeros(Leng+1,Leng+1);
% n=zeros(1,Leng+1);
% f=zeros(1,Leng+1);
% %y=zeros(1,Leng+1);
% itta=zeros(1,Leng+1);
% A=zeros(Leng+1,Leng+1);
% hh=zeros(1,Leng+1);

sum1=0;
for k=1:(Leng)
     
    
    g(k,k)=(((str_z(k).Umax * awz)) /(str_z(k).cz)) * (x(Leng+1) - x(k));
    
    L(k,k)=str_z(k).L;
    
    sum1=sum1 + str_z(k).Umax * (x(Leng+1)-x(k))*U{mm}(k);
        itta(k)=(str_z(k).az/str_z(k).cz)*Ta - ((h*str_z(k).Ad)/str_z(k).cz)*str_z(k).T1;
    sum=0;
    sum2=0;
    sum3=0;
    for c=1:(length(str_z(k).connected))     %sum for A
          sum=sum+str_z(k).az_ij(c)/str_z(k).cz;
          sum2=sum2+(str_z(k).az_ij(c)/str_z(k).cz)*x(str_z(k).connected(c));
          sum3 = sum3 + sign( x(str_z(k).connected(c)) - x(k) ) * str_z(k).paths.Ad_ij(1) * x(k) * sqrt(2 * (Cp - Cv) * abs(x(str_z(k).connected(c)) - x(k)));
    end
     
    A(k,k)=(h*str_z(k).Ad-(str_z(k).az))/str_z(k).cz - sum;   
    hh(k) = sum2 + ((p_air * Cp)/str_z(k).cz) * sum3;    
    rho(k) = (str_z(k).az/str_z(k).cz)*amp_out*sin(f_out*t);
end


A(Leng+1, Leng+1)= -aw/Cw;%aw


g(Leng+1, Leng+1)= (Ustmax/Cw)* (1+(COPmax-1)*(1-((x(Leng+1)-To))/DTmax));

hh(Leng+1) = (awz/Cw)* sum1; %aw

itta(Leng+1) = (aw/Cw) * Tpl;%aw

rho(Leng+1) = (awz/Cw)*amp_s*sin(f_s*t);

%system equations (zones & storage tank)
 dx1=x(1:Leng+1)'*A(1:Leng+1, 1:Leng+1)...
     + U{mm}*g(1:Leng+1, 1:Leng+1) + hh(1:Leng+1) + itta(1:Leng+1) + rho(1:Leng+1);


sum1=0;
for k=2*Leng+3 :3*Leng+2
    m=k-(2*Leng+2);
    g(k,k)=(((str_z(m).Umax * awz)) /(str_z(m).cz)) * (x(Leng+1) - x(k));

    sum1=sum1 + str_z(m).Umax * (x(Leng+1)-x(k))*U{mm}(m);
        itta(k)=(str_z(m).az/str_z(m).cz)*Ta - ((h*str_z(m).Ad)/str_z(m).cz)*str_z(m).T1;
    sum=0;
    sum2=0;
    sum3=0;
    for c=1:(length(str_z(m).connected))     %sum for A
          sum=sum+str_z(m).az_ij(c)/str_z(m).cz;
          sum2=sum2+(str_z(m).az_ij(c)/str_z(m).cz)*x(str_z(m).connected(c));
          sum3 = sum3 + sign( x(str_z(m).connected(c)) - x(k) ) * str_z(m).paths.Ad_ij(1) * x(k) * sqrt(2 * (Cp - Cv) * abs(x(str_z(m).connected(c)) - x(k)));
    end
     
    A(k,k)=(h*str_z(m).Ad-(str_z(m).az))/str_z(m).cz - sum;   
    hh(k) = sum2 + ((p_air * Cp)/str_z(m).cz) * sum3;    

end


A(3*Leng+3, 3*Leng+3)= -aw/Cw;%aw
g(3*Leng+3, 3*Leng+3)= (Ustmax/Cw)* (1+(COPmax-1)*(1-((x(Leng+1)-To))/DTmax));
hh(3*Leng+3) = (awz/Cw)* sum1; %aw
itta(3*Leng+3) = (aw/Cw) * Tpl;%aw


 
% storage tank observer gain
L(Leng+1,Leng+1)=L_s;


% obsevers (zones & storage tank)
dxe= x(2*Leng+3:3*Leng+3)'*A(2*Leng+3:3*Leng+3, 2*Leng+3:3*Leng+3)...
    + U{mm}*g(2*Leng+3:3*Leng+3, 2*Leng+3:3*Leng+3)...
    + hh(2*Leng+3:3*Leng+3) ...
    + itta(2*Leng+3:3*Leng+3)...
    + (y(1:Leng+1)-x(2*Leng+3:3*Leng+3)')*L;



E{mm}(1:Leng+1)=y(1:Leng+1)-(x(2*Leng+3:3*Leng+3))';

sum5=0;
for l=1:Leng
    sum5=sum5 + str_z(l).Umax * (n_bar_s + str_z(l).n_bar) * abs(U{mm}(l));
    sum6=0;
    sum7=0;
    for d=1:(length(str_z(l).connected))     %sum for A
 %       if y(l)==y(d)
            h_bar(l) = ( abs(y(l)) + str_z(l).n_bar) * sqrt( str_z(d).n_bar + str_z(l).n_bar);    
%         else
%             h_bar(l) =max( abs(y(l)) * sqrt( abs( y(d) - y(l) ))  ,  str_z(l).n_bar* (( abs(2*(y(d)) - 3*(y(l)))) / (2 * abs( y(d) - y(l)))) + str_z(d).n_bar* ( abs(y(l)) / (2 * abs( y(d) - y(l)))) + str_z(l).he_bar  );     
 %       end
        
         sum6 = sum6 + str_z(l).az_ij(d) * (str_z(l).n_bar  + str_z(d).n_bar);
         sum7 = sum7 + str_z(l).paths.Ad_ij(1) * h_bar(l);
    end
    
    
   % TF1=tf([0 str_z(l).p] , [1 str_z(l).lambda]);
    
    e1= abs( str_z(l).L ) * str_z(l).n_bar + str_z(l).r_bar...
        + ((str_z(l).Umax * awz) / str_z(l).cz) * (str_z(l).n_bar + n_bar_s) * abs(U{mm}(l))...
        + (1/str_z(l).cz) * sum6...
        + ((p_air * Cp)/str_z(l).cz) * sqrt( 2 * abs(Cp - Cv)) * sum7;
    
    funi=@(t) str_z(l).p*exp(-str_z(l).lambda * t)*e1;

% threshold of zone i
    Ebar{mm}(l) = str_z(l).p * exp(-str_z(l).lambda * t)*str_z(1).x_bar...
             +str_z(l).n_bar...
             +integral(funi,0,tt{mm});
%              +lsim(TF1,e1,t);
    
end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% replace length(t)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% with mm

%TF =tf([0 p_s] , [1  lambda_s]);

es= abs((Ustmax * (COPmax-1))/(Cw * DTmax)) * n_bar_s * abs( U{mm}(Leng+1) ) + abs(L_s) * n_bar_s;

funs=@(t) p_s*exp(-lambda_s * t)*es;

% threshold of storage tank

Ebar{mm}(Leng+1) =p_s * exp(-lambda_s * t) * x_bar_s ...
                    + integral(funs,0,tt{mm})...
                    + (aw/(Cw*1000)) *sum5 + r_bar_s...
                    + n_bar_s;
   
                
% Decision  

%dd=zeros(1,Leng+1);

for k=1:Leng+1        
    if abs(E{mm}(k))>Ebar{mm}(k)
        D{mm}(k)=1;
        %dd{mm}(k)=1;
    else
        D{mm}(k)=0;
    end
end
%dx1=x(1:Leng+1)'*A((1:Leng+1),1:Leng+1) + U*g((1:Leng+1),(1:Leng+1)) + hh((1:Leng+1)) + itta((1:Leng+1)); %state space equation of the system
% dx1 -- states of the system dx1+n
% U   -- control inputs 
% dxe -- states estimations of the observers
% U
% if t==time
%    save('var.mat','U')
% end

dx=[dx1 U{mm} dxe y(1:Leng+1)]';

%    else
        mm=mm+1;
 %   end
%dx=[dx1+n U]';
%dx(Leng+1)
%U(1)
% As=A;
% gs=g;
% hhs=hh;
% ittas=itta;
% state_s=x(1:Leng+1);
% 
% save system_Data.mat As gs hhs ittas state_s


% dx1=A*x(1:Leng+1) + g*U' + hh' + itta';
% 
% 
% dx2=[dx1+n' U'];
% 
% dx=dx2';
end
