function [phip_r, phip_i, dphi_dx_r, dphi_dx_i] = bem_quantity(col_pts, bndry_pts, N, Nx, psi, dpsi_dx, amp, omega, flow)

  % Unpackage needed flow quantities
  pref = flow.pref;
  rhoref = flow.rhoref;
  gamma = flow.gamma;
  Tref = flow.Tref;
  M = flow.M;
  La = flow.La;

  R = pref / (rhoref*Tref);
  aref = sqrt(Gamma*R*Tref);
  Uinf = M*aref;

  % Compute BEM solution on geometry
  [phibar_prime, pbar_prime] = bem3d(col_pts, bndry_pts, Nx, psi, dpsi_dx, amp, omega, flow)

  % Compute derivatives in streamwise direction
  Ny = N / Nx; 
  dx = bndry_pts(2,1) - bndry_pts(1,1);
  dphi_dx_r = zeros(N,1);
  dphi_dx_i = zeros(N,1);

  for l = 1:Ny
    for m = 1:Nx

      n = (l-1)*Nx + m;
      if ((m>1) && (m<Nx))
        dphi_dx_r(n) = (phi_r(n+1) - phi_r(n-1)) / (2*dx);
        dphi_dx_i(n) = (phi_i(n+1) - phi_i(n-1)) / (2*dx);
      elseif (m==1)
        dphi_dx_r(n) = (phi_r(n+1) - phi_r(n)) / dx;
        dphi_dx_i(n) = (phi_i(n+1) - phi_i(n)) / dx;
      else
        dphi_dx_r(n) = (phi_r(n) - phi_r(n-1)) / dx;
        dphi_dx_i(n) = (phi_i(n) - phi_i(n-1)) / dx;
      end

    end
  end

end
