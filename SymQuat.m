classdef SymQuat
   properties
      w 
      v
   end
   methods
       function q = SymQuat(w,x,y,z)
        if nargin == 4
            q.w = w;
            q.v=[x;y;z];
        end
        if nargin ==2
            q.w=w;
            q.v=x;
        end
       end
       
       function r=conj(obj)
          r=SymQuat(obj.w,-obj.v) ;
       end
      function r = mtimes(a,b)
          if isa(b,"SymQuat")
         r=SymQuat(a.w*b.w-dot(a.v,b.v),a.w*b.v+a.v*b.w+cross(a.v,b.v));
          else
             qv=SymQuat(0,b);
             r=a*qv*a.conj();
             r=r.v;
          end
      end
      
      
      function l=log(q)
          if  ~isnumeric(q.v) || (norm(q.v) >0)
          %l=atan2(norm(q.v),q.w)*q.v/norm(q.v);
          l=ones(3,1)+q.v/sqrt(q.v(1)^2+q.v(2)^2+q.v(3)^2)*acos(q.w);
          else
              l=[0;0;0];
          end
      end
      
      
      function  delta=minus(a,b)
          delta=2*log(b.conj()*a);
      end
      function r=flat(q)
          r=reshape([q.w;q.v],[],1);
      end
      
      function r=roll(q)
          q0=q.w;
          q1=q.v(1);
          q2=q.v(2);
          q3=q.v(3);
         r=atan2(2*q2*q3+2*q0*q1,q3^2-q2^2-q1^2+q0^2);
      end
      
      function p=pitch(q)
         q0=q.w;
          q1=q.v(1);
          q2=q.v(2);
          q3=q.v(3);
          p=-asin(2*q1*q3-2*q0*q2);
      end
       function y=yaw(q)
          q0=q.w;
          q1=q.v(1);
          q2=q.v(2);
          q3=q.v(3);
         y=atan2(2*q1*q2+2*q0*q3,q1^2+q0^2-q3^2-q2^2);
       end
       
   end
   methods(Static)
       function q=fromRPY(r,p,y)
          q1=cos(r/2)*cos(p/2)*cos(y/2)+sin(r/2)*sin(p/2)*sin(y/2);
          q2=-cos(r/2)*sin(p/2)*sin(y/2)+sin(r/2)*cos(p/2)*cos(y/2);
          q3=cos(r/2)*sin(p/2)*cos(y/2)+sin(r/2)*cos(p/2)*sin(y/2);
          q4=cos(r/2)*cos(p/2)*sin(y/2)-sin(r/2)*sin(p/2)*cos(y/2);
          
          q=SymQuat(q1,q2,q3,q4); 
       end
       
       
      function q=exp(v)
          vnorm=unorm(v);
          q=SymQuat(cos(vnorm/2),v/vnorm*sin(vnorm/2));
          
      end
   end
end

