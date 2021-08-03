function SS = steadyState_sc(v, dv, w, dw, dt, steps)

    SS = zeros(2,steps);
    SS(1,1) = v(1,1);
    SS(2,1) = w(1,1);
    for i=2:steps
        SS(1,i) = v(1,1) + dt*dv(v(1,1), w(1,1), 0);
        SS(2,i) = w(1,1) + dt*dw(v(1,1), w(1,1));
        v(1,1) = SS(1,i);
        w(1,1) = SS(2,i);
    end

    
end