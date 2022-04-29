function [check_result, new_Omega]=wheel_chair(meas_eq, start_Omega)
createArgs();
qz=sym("qz","real"); %For any reason evalin does not overwrite qz factorization function
py=sym("py","real"); %same goes for py
%cprintf('*black','State order: \n');
%disp(state);
syms wx wy wz real;
syms w_px w_py w_pz w_qx w_qy w_qz w_qw w_vx w_vy w_vz w_bax w_bay w_baz w_bgx w_bgy w_bgz real
w_q=SymQuat(w_qw,w_qx,w_qy,w_qz);
w_v=[w_vx; w_vy; w_vz];
w_ba=[w_bax; w_bay; w_baz];
w_bg=[w_bgy; w_bgy; w_bgz];
w_p=[w_px; w_py; w_pz];
wheel_state=reshape([qw;qx;qy;qz;v;ba;bg; w_qw; w_qx;w_qy;w_qz;w_v;w_ba;w_bg],1,[]);

w=[wx;wy;wz];
dotq=simplify(0.5*flat(q*SymQuat(0,w-bg)));
w_dotq=simplify(0.5*flat(w_q*SymQuat(0,w-w_bg)));

%Part of the dynamic model that corresponds to the state
f0=[simplify(subs(dotq,w,zeros(3,1)));
    q*(-ba)+g;
    zeros(6,1);
    simplify(subs(w_dotq,w,zeros(3,1)));
    w_q*(-w_ba)+g;
    zeros(6,1);
    ];
%Part that corresponds to the acceleration (local frame)
syms ax ay az real;
a=[ax;ay;az];
wa=q*a;
w_wa=w_q*a;
f1=[zeros(4,1);simplify(diff(wa,ax));zeros(6,1);zeros(13,1)];
f2=[zeros(4,1);simplify(diff(wa,ay));zeros(6,1);zeros(13,1)];
f3=[zeros(4,1);simplify(diff(wa,az));zeros(6,1);zeros(13,1)];
%Part that corresponds to the angular rate (local frame)
f4=[simplify(diff(dotq,wx));zeros(9,1);zeros(13,1)];
f5=[simplify(diff(dotq,wy));zeros(9,1);zeros(13,1)];
f6=[simplify(diff(dotq,wz));zeros(9,1);zeros(13,1)];

f7=[zeros(13,1);zeros(4,1);simplify(diff(w_wa,ax));zeros(6,1)];
f8=[zeros(13,1);zeros(4,1);simplify(diff(w_wa,ay));zeros(6,1)];
f9=[zeros(13,1);zeros(4,1);simplify(diff(w_wa,az));zeros(6,1)];
%Part that corresponds to the angular rate (local frame)
f10=[zeros(13,1);simplify(diff(w_dotq,wx));zeros(9,1)];
f11=[zeros(13,1);simplify(diff(w_dotq,wy));zeros(9,1)];
f12=[zeros(13,1);simplify(diff(w_dotq,wz));zeros(9,1)];
input_dynamics=[f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12];
input_names=["ax","ay","az","wx","wy","wz","w_ax","w_ay","w_az","w_wx","w_wy","w_wz"];
  
syms lever_x lever_y lever_z real;
lever=[lever_x; lever_y; lever_z];
if nargin < 2 %Use start Omega if given
%Add Unit norm constraint
meas_functions;
z=[(qw^2+qx^2+qy^2+qz^2)/2;(w_qw^2+w_qx^2+w_qy^2+w_qz^2)/2;v   w_v; body_velocity(2:3)];
%Omega 0
Omega=simplify(jacobian(z,wheel_state));
start_line=3; %never need to check unit norm constraint as it can only give 1dof
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
check_names=[  "Body velocity" "Non-Holonomic" "Forward v" "|v|" "vz" "vx" "vy" "Ground velocity" "Gravity direction" "Roll" "Pitch" "Yaw" "bax" "bay" "baz" "bgx" "bgy" "bgz"];
num_checks=size(check_names,2); 
global check_results;
num_inputs=size(input_dynamics,2)
check_results=-1*ones(num_checks,2^num_inputs);

[check_results(:,2^num_inputs),new_Omega]=checkObservability(2^6-1,f0,input_dynamics,Omega,num_checks,start_line,num_inputs);
checkPartialObservability(2^num_inputs-1,f0,input_dynamics,Omega,num_checks, check_names,input_names, start_line,num_inputs);
%max_unobs_dof=check_results(24,2^6);
check_result=check_results(:,2^num_inputs);
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
   function [check_result, Omega]=checkObservability(f_index,f0,input_dynamics,Omega,num_checks, start_line,num_inputs)
    createArgs();
    qz=sym("qz","real"); %For any reason evalin does not overwrite qz factorization function
    py=sym("py","real"); %same goes for py
    syms w_px w_py w_pz w_qx w_qy w_qz w_qw w_vx w_vy w_vz w_bax w_bay w_baz w_bgx w_bgy w_bgz real
    w_q=SymQuat(w_qw,w_qx,w_qy,w_qz);
    w_v=[w_vx; w_vy; w_vz];
    w_ba=[w_bax; w_bay; w_baz];
    w_bg=[w_bgy; w_bgy; w_bgz];
    w_p=[w_px; w_py; w_pz];
    wheel_state=reshape([qw;qx;qy;qz;v;ba;bg; w_qw; w_qx;w_qy;w_qz;w_v;w_ba;w_bg],1,[]);
    mask=NumberToMask(f_index,num_inputs);
    state_size=size(Omega,2);
    dynamics=[f0 input_dynamics(:,mask==1)];
    check_result=zeros(num_checks,1);
new_Omega=Omega;
old_rank=rank(Omega);
current_rank=old_rank;
while 1
for o_index=start_line:size(Omega,1)
    for f_index=1:size(dynamics,2)
        lie=Omega(o_index,:)*dynamics(:,f_index);
        gradient_lie=jacobian(lie, wheel_state);
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
state_size=size( wheel_state,2);
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
check_result(1:num_checks)=[
    SingleCheck(body_velocity,Omega, wheel_state,"Body velocity");
    SingleCheck(body_velocity(2:3),Omega, wheel_state,"Non Holonomic");
    SingleCheck(body_velocity(1),Omega, wheel_state,"Forward velocity");
    SingleCheck( vx^2+ vy^2+ vz^2,Omega, wheel_state,"|v|");
    SingleCheck(vz,Omega, wheel_state); 
    SingleCheck(vx,Omega, wheel_state); 
    SingleCheck(vy,Omega, wheel_state);
    SingleCheck( vx^2+ vy^2,Omega, wheel_state,"Ground Velocity");
    SingleCheck( q.conj()* g,Omega, wheel_state,"Gravity direction");
    SingleCheck( q.roll,Omega, wheel_state,"Roll");
    SingleCheck( q.pitch,Omega, wheel_state,"Pitch");
    SingleCheck( q.yaw,Omega, wheel_state,"Yaw");
    SingleCheck(bax,Omega, wheel_state);
    SingleCheck(bay,Omega, wheel_state);
    SingleCheck(baz,Omega, wheel_state);
    SingleCheck(bgx,Omega, wheel_state);
    SingleCheck(bgy,Omega, wheel_state);
    SingleCheck(bgz,Omega, wheel_state);
 ];
%unobs_dof]; %unobs_dof at the end
%Omega

   end

   function check_result=checkPartialObservability(f_index,f0,input_dynamics,Omega,num_checks, check_names,input_names,start_line,num_inputs)
   global check_results;
   if (f_index ~=2^num_inputs-1) && (check_results(1,f_index+1) ~=-1)
       check_result=check_results(:,f_index+1);
       return
   end
   mask=NumberToMask(f_index,num_inputs);
   pred_result=zeros(num_checks,1);
   for index=1:num_inputs
       if mask(index) ==1
           mask(index)=0;
           pred_result=max(pred_result,checkPartialObservability(MaskToNumber(mask),f0,input_dynamics,Omega,num_checks, check_names,input_names,start_line,num_inputs));
           mask(index)=1;
       end
   end
   if pred_result==check_results(2^num_inputs) %Check if maximum observability is already reached
      check_result=pred_result;
   else
       if f_index==2^num_inputs-1
           check_result=check_results(:,f_index+1);
       else
           check_result=checkObservability(f_index,f0,input_dynamics,Omega,num_checks,start_line,num_inputs);
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
