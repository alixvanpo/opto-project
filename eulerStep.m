function S_next = eulerStep(v, dv, w, dw, iapp, dx, dt)

    lap = del2_noflux(v, dx);
    S_next(:,:,1) = v + dt*dv(v, w, iapp, lap);
    S_next(:,:,2) = w + dt*dw(v, w);
    
    
end