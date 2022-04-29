

function [check_result, new_Omega, new_Momega]=obs_proof_Manifold(meas_eq, start_Omega)
createArgs();
qz=sym("qz","real"); %For any reason evalin does not overwrite qz factorization function
py=sym("py","real"); %same goes for py
%cprintf('*black','State order: \n');
%disp(state);
syms wx wy wz real;
w=[wx;wy;wz];
dotq=simplify(0.5*flat(q*SymQuat(0,w-bg)));

%Part of the dynamic model that corresponds to the state
f0=[simplify(subs(dotq,w,zeros(3,1)));
    v;
    q*(-ba)+g;
    zeros(6,1)
    ];
%Part that corresponds to the acceleration (local frame)
syms ax ay az real;
a=[ax;ay;az];
wa=q*a;
f1=[zeros(7,1);simplify(diff(wa,ax));zeros(6,1)];
f2=[zeros(7,1);simplify(diff(wa,ay));zeros(6,1)];
f3=[zeros(7,1);simplify(diff(wa,az));zeros(6,1)];
%Part that corresponds to the angular rate (local frame)
f4=[simplify(diff(dotq,wx));zeros(12,1)];
f5=[simplify(diff(dotq,wy));zeros(12,1)];
f6=[simplify(diff(dotq,wz));zeros(12,1)];
input_dynamics=[f1 f2 f3 f4 f5 f6 ];
input_names=["ax","ay","az","wx","wy","wz"];
  
if nargin < 2 %Use start Omega if given
%Add Unit norm constraint
z=[meas_eq];
%Omega 0
Omega=simplify(jacobian(z,state));
start_line=1; 
%Check if omega0 is span
if rank(Omega) ~= size(Omega,1)
    fprintf("Omega 0 is rank deficient. %i measurement equations are linear dependent from others.", size(Omega,1)-rank(Omega));
end
else
    Omega=[start_Omega;simplify(jacobian(meas_eq,state))];
    start_line=size(start_Omega,1);
end
%Omega
% Ganz verwerfen
%Für alle möglichen Kombinationen durchrechen (2^6) und speichern
%dann Baumartig durchgehen und testen, ob ein Vorgänger bereits für
%Unbeobachtbarkeit sorgt
check_names=["px" "py" "pz"  "Body velocity" "Non-Holonomic" "Forward v" "|v|" "vz" "vx" "vy" "Ground velocity" "Gravity direction" "Roll" "Pitch" "Yaw" "bax" "bay" "baz" "bgx" "bgy" "bgz"];
num_checks=size(check_names,2); 
global check_results;
check_results=-1*ones(num_checks,2^6);

Momega=jacobian(subs(z,flat(q),flat(q*q.exp(delta))), delta_state);
[check_results(:,2^6),new_Omega, new_Momega]=checkObservability(2^6-1,f0,input_dynamics,Omega,Momega,num_checks,start_line);
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
   function [check_result, Omega,Momega]=checkObservability(f_index,f0,input_dynamics,Omega, Momega,num_checks, start_line)
    createArgs();
    py=sym("py","real"); %same goes for py
    mask=NumberToMask(f_index);
    state_size=size(Omega,2);
    dynamics=[f0 input_dynamics(:,mask==1)];
    check_result=zeros(num_checks,1);
new_Omega=Omega;
new_Momega=Momega;
old_rank=rank(Momega);
current_rank=old_rank;
while 1
for o_index=start_line:size(Omega,1)
    for f_index=1:size(dynamics,2)
        lie=Omega(o_index,:)*dynamics(:,f_index);
        gradient_lie=jacobian(subs(lie,flat(q),flat(q*q.exp(delta))), delta_state);
        new_rank=rank([new_Momega;gradient_lie]);
        if new_rank > current_rank
            new_Omega=[new_Omega;simplify(jacobian(lie,state))];
            new_Momega=[new_Momega;gradient_lie]
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
Momega=new_Momega;
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
state_size=size( state,2);
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
check_result(1:num_checks)=[MSingleCheck( px,Momega);
    MSingleCheck( py,Momega);
    MSingleCheck( pz,Momega);
    MSingleCheck(body_velocity,Momega,"Body velocity");
    MSingleCheck(body_velocity(2:3),Momega,"Non Holonomic");
    MSingleCheck(body_velocity(1),Momega,"Forward velocity");
    MSingleCheck( vx^2+ vy^2+ vz^2,Momega,"|v|");
    MSingleCheck(vz,Momega); 
    MSingleCheck(vx,Momega); 
    MSingleCheck(vy,Momega);
    MSingleCheck( vx^2+ vy^2,Momega,"Ground Velocity");
    MSingleCheck( q.conj()* g,Momega,"Gravity direction");
    MSingleCheck( q.roll,Momega,"Roll");
    MSingleCheck( q.pitch,Momega,"Pitch");
    MSingleCheck( q.yaw,Momega,"Yaw");
    MSingleCheck(bax,Momega);
    MSingleCheck(bay,Momega);
    MSingleCheck(baz,Momega);
    MSingleCheck(bgx,Momega);
    MSingleCheck(bgy,Momega);
    MSingleCheck(bgz,Momega);
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
