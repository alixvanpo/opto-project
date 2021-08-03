function lap = del2_noflux(v, dx)
    %%%%%%%%%%%%%%%%
    % Compute the Laplacian of matrix v using  no-flux boundary conditions
    %%%%%%%%%%%%%%%%
    
    % Append extra row on top and bottom, and extra column on left and
    % right of v, using same values as were on edge.
    
    % Append col to left and right
    v_temp = cat(2, v(:,1), v, v(:,end));
    
    % Append row to top and bottom
    v_temp2 = cat(1, v_temp(1,:), v_temp, v_temp(end,:));
    
    
    % Use del2 to get the Laplacian (that uses unwanted linear
    % extrapolation on boundary)
    lap_LE = 4*del2(v_temp2, dx);
   
    % Remove boundary values
    lap = lap_LE(2:end-1, 2:end-1);

end