classdef hvac
    
   properties
      x1,Leng,t1;
   end
   
    methods
%       function obj = hvac(varargin)
%       end
%       function obj=getnum(obj,num)
%           obj.num=num;
%       end
        
        %Constructor
        function obj=hvac(x1,Leng,t1)
            
            obj.x1=x1;
            obj.Leng=Leng;
            obj.t1=t1;
            
        
        end

      function plotresults(obj,num)
         if (num > 1 && num <=(obj.Leng+1))
             if (num == 1)
                figure
                subplot(3,2,1)
                plot(obj.t1, obj.x1(:,num),'LineWidth',2);
                title('Water temperature in Storage tank','Interpreter','latex')
                ylabel('$T_{st}$','Interpreter','latex')
                xlabel('time (hours)','Interpreter','latex')
             else
                figure
                subplot(3,2,1)
                plot(obj.t1, obj.x1(:,num),'LineWidth',2);
                title('Air temperature in Zone','Interpreter','latex')
                ylabel('$T_{z}$','Interpreter','latex')
                xlabel('time (hours)','Interpreter','latex')
             end
         else
            error('Invalid Subsystem')
         end
         
      end
      
   end
   
end