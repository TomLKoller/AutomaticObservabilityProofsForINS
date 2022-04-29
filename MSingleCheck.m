function is_observable=MSingleCheck(toCheck,Momega,~)
    createArgs();
    augmented_Momega=[Momega;jacobian(subs(toCheck,flat(q),flat(q*q.exp(delta))), delta_state);];
    is_observable=  rank(augmented_Momega) ==rank(Momega);
end

