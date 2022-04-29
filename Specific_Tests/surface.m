function yaw_dot=surface(meas_eq, start_Omega)
createArgs();
qz=sym("qz","real"); %For any reason evalin does not overwrite qz factorization function
py=sym("py","real"); %same goes for py


    
syms v real
syms yaw(s)
assume(yaw(s),"real");
xdash=v*cos(yaw)
ydash=v*sin(yaw)
syms s real
syms x(s) y(s) z(x,y) 
r=[x(s); y(s); z(x(s),y(s))]
rv=diff(r,s);
rv=simplify(subs(rv,[diff(x,s),diff(y,s)],[xdash, ydash]), "IgnoreAnalyticConstraints",true)
syms roll(s)
assume(roll(s),"real");
theta=acos(rv(3)/unorm(rv))
syms pitch real
q=q.fromRPY(roll(s),theta,yaw(s))
%yaw_formula=simplify(q.yaw,"IgnoreAnalyticConstraints",true)
quat=SymQuat(qw,qx,qy,qz);
syms wx wy wz real;
w=[wx;wy;wz];
qdot=q*SymQuat(0,w);
yaw_jacob=jacobian(quat.yaw,[qw,qx,qy,qz]);
yaw_jacob=subs(yaw_jacob,[qw;qx;qy;qz],q.flat());
yaw_dot=simplify(yaw_jacob*qdot.flat(),"IgnoreAnalyticConstraints",true)
return
assume(x(s),"real");
assume(y(s),"real");
assume(z(x(s),y(s)),"real");
assume(r,"real");
T=diff(r,s)
T=simplify(subs(T,[diff(x,s),diff(y,s)],[xdash, ydash]), "IgnoreAnalyticConstraints",true)
assume(T,"real")
N=diff(T,s);
assume(N,"real");
N=simplify(N/unorm(N), "IgnoreAnalyticConstraints",true)
B=simplify(cross(T,N), "IgnoreAnalyticConstraints",true)
assume(B,"real");
torsion=diff(B,s);
assume(torsion,"real");
torsion=simplify(-unorm(torsion), "IgnoreAnalyticConstraints",true)
curvature=diff(T,s);
assume(curvature,"real");
curvature=simplify(unorm(curvature), "IgnoreAnalyticConstraints",true)
pretty(curvature)


return

syms yaw real
state=reshape([qw;qx;qy;qz; lambda],1,[]);
%cprintf('*black','State order: \n');
%disp(state);

%q=SymQuat(cos(lambda/2),[0;0;sin(lambda/2)]);
dotq=simplify(0.5*flat(q*SymQuat(0,w)))

%Part of the dynamic model that corresponds to the state
f0=[zeros(5,1);
    ];
%Part that corresponds to the acceleration (local frame)
syms ax ay az real;
a=[ax;ay;az];
wa=q*a;
%Part that corresponds to the angular rate (local frame)
f1=[simplify(diff(dotq,wx));zeros(1,1)];
f2=[simplify(diff(dotq,wy));zeros(1,1)];
f3=[simplify(diff(dotq,wz));1];
f4=[zeros(4,1); 1/sqrt(nx(lambda)^2+ny(lambda)^2+nz(lambda)^2)];
%f4=[zeros(4,1); 1];
input_dynamics=[f1 f2 f3 f4]
input_names=["wx","wy","wz", "|w|"];
  
if nargin < 2 %Use start Omega if given
%Add Unit norm constraint
%z=[(qw^2+qx^2+qy^2+qz^2)/2;n/sqrt(nx(lambda)^2+ny(lambda)^2+nz(lambda)^2)];
syms fqx(lambda) fqy(lambda) fqz(lambda) fqw(lambda) real

syms a b c alph(lambda) real
alph(lambda)=lambda;
%z=[(qw^2+qx^2+qy^2+qz^2)/2;qw-cos(alph(lambda)/2);qx;qy; qz-sin(alph(lambda)/2)];
z=[(qw^2+qx^2+qy^2+qz^2)/2;qw-fqw(lambda);qx-fqx(lambda);qy-fqy(lambda)];
a=0
b=0
c=1
q_lambda=SymQuat(cos(alph(lambda)/2),[a*sin(alph(lambda)/2);b*sin(alph(lambda)/2);c*sin(alph(lambda)/2)]);
%z=[(qw^2+qx^2+qy^2+qz^2)/2;SymQuat(qw,qx,qy,qz)-q_lambda]
%Omega 0
Omega=simplify(jacobian(z,state));
start_line=2; %Never need to check unit norm constraint as it can only give 1dof
%Check if omega0 is span
if rank(Omega) ~= size(Omega,1)
    fprintf("Omega 0 is rank deficient. %i measurement equations are linear dependent from others.", size(Omega,1)-rank(Omega));
end
else
    Omega=[start_Omega;simplify(jacobian(meas_eq,state))];
    start_line=size(start_Omega,1);
end
Omega
% Ganz verwerfen
%Für alle möglichen Kombinationen durchrechen (2^6) und speichern
%dann Baumartig durchgehen und testen, ob ein Vorgänger bereits für
%Unbeobachtbarkeit sorgt
check_names=[ "Roll" "Pitch" "Yaw" "lambda"];
num_checks=size(check_names,2); 
global check_results;
check_results=-1*ones(num_checks,2^6);

[check_results(:,2^6),new_Omega]=checkObservability(2^6-1,f0,input_dynamics,Omega,num_checks,state,start_line);
%checkPartialObservability(2^6-1,f0,input_dynamics,Omega,num_checks, check_names,input_names);
%max_unobs_dof=check_results(24,2^6);
check_result=check_results(:,2^6);
%Print results
%cprintf('*black','Check of each state: \n');
% for index=1:num_checks
%    if(check_results(index,2^6))
%    cprintf('green','%s is observable.\n',check_names(index));
%    else
%    cprintf('red','%s is not observable.\n',check_names(index));
%     end
%    if index==16    
%     cprintf('*black','Check of combined states: \n');
%    end
% 
% end
end
   function [check_result, Omega]=checkObservability(f_index,f0,input_dynamics,Omega,num_checks,state, start_line)
    s=state;
    createArgs();
    state=s;
    py=sym("py","real"); %same goes for py
    state_size=size(Omega,2);
    dynamics=[f0 input_dynamics];
    check_result=zeros(num_checks,1);
    syms lambda
new_Omega=Omega;
old_rank=rank(Omega);
current_rank=old_rank;
while 1
for o_index=start_line:size(Omega,1)
    %o_index
    for f_index=1:size(dynamics,2)
        %f_index
        lie=Omega(o_index,:)*dynamics(:,f_index);
        gradient_lie=jacobian(lie, state);
        new_rank=rank([new_Omega;gradient_lie]);
        if new_rank > current_rank
            new_Omega=[new_Omega;simplify(gradient_lie)];
            current_rank=new_rank;
        end
        if new_rank ==state_size
            break
        end
    end
    if new_rank ==state_size
            break
    end
end
new_rank=current_rank;
if (old_rank==new_rank) 
    break
end
start_line=size(Omega,1);
Omega=new_Omega;
if new_rank==state_size
    break
end
old_rank=new_rank;
end
%Give basic information
%cprintf('*black','Observable Codistribution: \n');
%disp(Omega)
%cprintf('*black','Null space of Omega (transposed):  \n');
%disp(null(Omega).')
%Check how many DOFs are unobservable

rank_omega=old_rank;
unobs_dof=state_size-rank_omega;
if f_index==2^6-1
if unobs_dof==0
    fprintf("The system is fully observable.\n");
else
    fprintf("%i DOF are unobservable.\n",unobs_dof);
end
end
%Check each state
body_velocity=q.conj()*v;
check_result(1:num_checks)=[SingleCheck( q.roll,Omega, state,"Roll");
    SingleCheck( q.pitch,Omega, state,"Pitch");
    SingleCheck( q.yaw,Omega, state,"Yaw");
    SingleCheck(lambda,Omega, state);
 ];
%unobs_dof]; %unobs_dof at the end
%Omega

   end

   function check_result=checkPartialObservability(f_index,f0,input_dynamics,Omega,num_checks, check_names,input_names)
   global check_results;
   if (f_index ~=2^6-1) && (check_results(1,f_index+1) ~=-1)
       check_result=check_results(:,f_index+1);
       return
   end
   mask=NumberToMask(f_index);
   pred_result=zeros(num_checks,1);
   for index=1:6
       if mask(index) ==1
           mask(index)=0;
           pred_result=max(pred_result,checkPartialObservability(MaskToNumber(mask),f0,input_dynamics,Omega,num_checks, check_names,input_names));
           mask(index)=1;
       end
   end
   if pred_result==check_results(2^6) %Check if maximum observability is already reached
      check_result=pred_result;
   else
       if f_index==2^6-1
           check_result=check_results(:,f_index+1);
       else
           check_result=checkObservability(f_index,f0,input_dynamics,Omega,num_checks);
       end
       newly_observable=1:num_checks;
       newly_observable=newly_observable(check_result > pred_result);
       for index=1:size(newly_observable,2)
           if sum(mask) > 0
           cprintf('black',"%s is observable from f0,%s.\n",check_names(newly_observable(index)),join(input_names(mask==1),','));
           else
           cprintf('black',"%s is observable from f0.\n",check_names(newly_observable(index)));
           end
       end
   end
   check_results(:,f_index+1)=check_result;
   
   end
