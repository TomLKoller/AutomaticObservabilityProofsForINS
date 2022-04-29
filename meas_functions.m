createArgs();
check_names=["px" "py" "pz"  "Body velocity" "Non-Holonomic" "Forward v" "|v|" "vz" "vx" "vy" "Ground velocity" "Gravity direction" "Roll" "Pitch" "Yaw" "bax" "bay" "baz" "bgx" "bgy" "bgz"];
%TODO non holonomic constraint
%Todo Knowledge about acceleration and angular rate
body_velocity=q.conj()*v;
meas_functions_={px; 
    py ; 
    pz; 
    [body_velocity];
    [body_velocity(2:3)];
    body_velocity(1); 
    vx^2+ vy^2+ vz^2;
    vz;
    vx; 
    vy; 
    vx^2+ vy^2;
    [q.conj()* g];
    q.roll; 
    q.pitch; 
    q.yaw; 
    bax; 
    bay; 
    baz; 
    bgx; 
    bgy; 
    bgz; 
    };