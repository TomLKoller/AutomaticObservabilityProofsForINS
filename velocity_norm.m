syms px py pz real;
syms vx vy vz real;
syms qw qx qy qz real;
syms bax bay baz real;
syms bgx bgy bgz real;
p=[px;py;pz];
v=[vx;vy;vz];
q=SymQuat(qw,qx,qy,qz);
ba=[bax;bay;baz];
bg=[bgx;bgy;bgz];
state=reshape([p;v;qw;qx;qy;qz;ba;bg],1,[]);
cprintf('*black','State order: \n');
disp(state);
g=[0;0;-9.81];

syms wx wy wz real;
w=[wx;wy;wz];
dotq=simplify(0.5*flat(q*SymQuat(0,w-bg)));

%Part of the dynamic model that corresponds to the state
f0=[v;
    q*(-ba)+g;
    subs(dotq,w,zeros(3,1));
    zeros(6,1)
    ];
%Part that corresponds to the acceleration (local frame)
syms ax ay az real;
a=[ax;ay;az];
wa=q*a;
f1=[zeros(3,1);diff(wa,ax);zeros(10,1)];
f2=[zeros(3,1);diff(wa,ay);zeros(10,1)];
f3=[zeros(3,1);diff(wa,az);zeros(10,1)];
%Part that corresponds to the angular rate (local frame)
f4=[zeros(6,1);diff(dotq,wx);zeros(6,1)];
f5=[zeros(6,1);diff(dotq,wy);zeros(6,1)];
f6=[zeros(6,1);diff(dotq,wz);zeros(6,1)];
dynamics=[f0 f1 f2 f3 f4 f5 f6 ];

%z=[atan2(y,x)-theta;atan2(y,x-1)-theta;]
%-----------Type your measurement equations here, concatenated with ;--------
%e.g. meas_eq=[p;q];
meas_eq=[vx^2+vy^2+vz^2];
%--------------------------
%Add Unit norm constraint
z=[qw^2+qx^2+qy^2+qz^2;meas_eq];
%Omega 0
Omega=simplify(jacobian(z,state));
%Check if omega0 is span
if rank(Omega) ~= size(Omega,1)
    fprintf("Omega 0 is rank deficient. %i measurement equations are linear dependent from others.", size(Omega,1)-rank(Omega));
end

start_line=1;
new_Omega=Omega;
while 1
for o_index=start_line:size(Omega,1)
    for f_index=1:7
        lie=simplify(Omega(o_index,:)*dynamics(:,f_index));
        gradient_lie=simplify(jacobian(lie,state));
        if rank([new_Omega;gradient_lie]) > rank(new_Omega)
            new_Omega=[new_Omega;gradient_lie];
        end
    end
end
if rank(Omega)==rank(new_Omega)
    break
end
start_line=size(Omega,1);
Omega=new_Omega;
end
%Give basic information
cprintf('*black','Observable Codistribution: \n');
disp(Omega)
cprintf('*black','Null space of Omega (transposed):  \n');
disp(null(Omega).')
%Check how many DOFs are unobservable
state_size=size(state,2);
rank_omega=rank(Omega);
unobs_dof=state_size-rank_omega;
if unobs_dof==0
    fprintf("The system is fully observable.\n");
    return
else
    fprintf("%i DOF are unobservable.\n",unobs_dof);
end
%Check each state
cprintf('*black','Check of each state: \n');
for index=1:state_size
   SingleCheck(state(index),Omega,state);
end
%Check Special states 
cprintf('*black','Check of combined states: \n');
SingleCheck(q.roll,Omega,state,"Roll");
SingleCheck(q.pitch,Omega,state,"Pitch");
SingleCheck(q.yaw,Omega,state,"Yaw");
SingleCheck(q.conj()*g,Omega,state,"Gravity direction");
SingleCheck(vx^2+vy^2+vz^2,Omega,state,"|v|");
SingleCheck(vx^2+vy^2,Omega,state,"Ground Velocity");
SingleCheck(q.conj()*v,Omega,state,"Body velocity");




