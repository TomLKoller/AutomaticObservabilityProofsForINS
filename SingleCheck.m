function is_observable=SingleCheck(toCheck,Omega,state,text)
    if nargin < 4 
        text=toCheck;
    end
    augmented_Omega=[Omega;jacobian(toCheck,state)];
    is_observable=  rank(augmented_Omega) ==rank(Omega);
end

