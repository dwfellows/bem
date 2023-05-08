function [psi, dpsi_dx] = beam_deformation(col_pts, La, Lb, Nx, N)

  % Define 2D data
 c1 = 4.73;
 psi2d = @(x,Lp) cos(c1.*x./Lp) - cosh(c1.*x./Lp) - ((cos(c1) - cosh(c1))/(sin(c1) - sinh(c1))).*(sin(c1.*x./Lp) - sinh(c1.*x./Lp));
 psi_max_2d = psi2d(L_panel/2, L_panel);
 psi2d = @(x,Lp) (1/psi_max_2d) .* ( cos(c1.*x./Lp) - cosh(c1.*x./Lp) - ((cos(c1) - cosh(c1))/(sin(c1) - sinh(c1))).*(sin(c1.*x./Lp) - sinh(c1.*x.    /Lp)) );
 dpsi_dx2d = @(x,Lp) (1/psi_max_2d) .* (-c1/Lp) .* ( ((cos(c1) - cosh(c1))/(sin(c1) - sinh(c1))).*(cos(c1.*x./Lp) - cosh(c1.*x./Lp)) + sin(c1.*x./    Lp) + sinh(c1.*x./Lp) );

  psi = zeros(N,1); dpsi_dx = zeros(N,1);
  for j = 1:(Ny-1)

    m = (j-1)*Nx + 1;
    loc_amp = psi2d(col_pts(m,2), Lb);

    strip_indices = (j-1)*Nx+1:(j-1)*Nx+Nx
    psi(strip_indices) = loc_amp.*psi2d(col_pts(strip_indices,1), La);
    dpsi_dx(strip_indices) = loc_amp.*dpsi_dx(col_pts(strip_indices,1), La); 

  end

end
