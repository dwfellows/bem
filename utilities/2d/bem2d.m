function [phibar_prime, pbar_prime] =  bem2d(col_pts, bndry_pts, psi, dpsi_dx, amp, omega, flow)

  % Utilize boundary element method to obtain spatial solution to reduced velocity perturbation potential
  % and reduced pressure fluctuation for two-dimensional flow.
  % This method assumes an oscillating beam exposed to an external freestream.

  % Inputs:
  %   col_pts: collocation points at which to solve for spatial quantities
  %   bndry_pts: boundary element demarcating points (col_pts assumed to be at midponts
  %   psi: mode shape of harmonic deformation of beam
  %   dpsi_dx: first derivative of mode shape psi
  %   amp: amplitude of harmonic deformation
  %   omega: frequency (in radians, omega = 2*pi*f) of harmonic deformation
  %   flow: flow quantities of freestream {p_in, rho_in, gamma, T_in, M}, generally matching a PC2 simulation

  % Outputs: spatial solutions to reduced velocity perturbation potential (phibar_prime) and reduced pressure fluctuation (pbar_prime)

  % Unpackage flow quantities
  pref = flow.pref;
  rhoref = flow.rhoref;
  gamma = flow.gamma;
  Tref = flow.Tref;
  M = flow.M;
  L = flow.L;

  % Compute dependent flow quantities
  R = pref / (rhoref*Tref);
  aref = sqrt(gamma*R*Tref);
  Uref = sqrt(pref / rhoref);
  Uinf = M*aref;

  % Compute quantities of interest to convected wave equation
  k = omega / aref;      % Wavenumber
  beta = sqrt(1 - M^2);  % Prandtl-Glauert factor

  % Compute quantities of interest to boundary element method
  N = length(col_pts);   % Number of collocation points / boundary elements
  dx = bndry_pts(2) - bndry_pts(1);
  
  % Compute panel boundary condition from flow and deformation quantities
  phibar_prime = zeros(1,N);
  %dphi_dn = (Uinf*amp*dpsi_dx + amp*i*omega*psi);  % original method that seemed to work
  dphi_dn = (-(i/2)*amp*Uinf*dpsi_dx + (1/2)*omega*amp*psi); % computed from complex-conjugate method

  % Compute reduced velocity perturbation potential
  if (M < 1)  % If subsonic -- follow Lacerda et al.
    for j = 1:N
      x = col_pts(j); y = 0;
      for l = 1:N
        phibar_prime(j) = phibar_prime(j) - dphi_dn(l)*gq2d(x,y,bndry_pts(l:l+1),k,M);  % Compute Green's function via Gaussian quadrature
      end
    end
  else % If supersonic -- follow Jones

    % Compute required quantities
    mu = asin(1/M);
    lambda = omega * L / Uinf;

    % Perform transformation
    % dphi_dn = L .* exp(i*omega*(1/cos(mu)^2).*col_pts./Uinf) .* dphi_dn;
    bndry_pts = bndry_pts .* tan(mu) ./ L;
    col_pts = col_pts .* tan(mu) ./ L;
    dphi_dn = L .* exp(i*omega*(L/(cos(mu)*sin(mu))).*col_pts./Uinf) .* dphi_dn;

    % for j = 1:N
    %   f = @(x) besselj(0,(lambda/cos(mu)) * (col_pts(j) - x));
    %   for l = 1:j
    %     phibar_prime(j) = phibar_prime(j) - dphi_dn(l)*integral(f,bndry_pts(l),bndry_pts(l+1));
    %   end
    %   phibar_prime(j) = phibar_prime(j) * exp(-i*lambda*(1/(cos(mu)*sin(mu)))*col_pts(j));
    % end

    % phibar_prime = zeros(1,N);

    for j = 1:N
      x = col_pts(j); y = 0;
      for l = 1:j  % Integration is only performed up to current control point (supersonic flow)
        phibar_prime(j) = phibar_prime(j) - dphi_dn(l)*gq2d(x,y,bndry_pts(l:l+1),lambda/cos(mu),M); % Compute Green's function via Gaussian quadrature
      end
    end

    % Transform velocity perturbation potential back to physical space
    phibar_prime = phibar_prime .* exp(-i * lambda * (1/(cos(mu)*sin(mu))) .* col_pts);

  end

  % Compute pressure fluctuation via finite differences
  pbar_prime = zeros(1,N);
  for j = 1:N
    if ((j>1) && (j<N)) % Interior point
      pbar_prime(j) = -rhoref*( Uinf*((phibar_prime(j+1) - phibar_prime(j-1))/(2*dx)) + i*omega*phibar_prime(j) );
    elseif (j==1) % Leading edge
      pbar_prime(j) = -rhoref*( Uinf*((phibar_prime(j+1) - phibar_prime(j))/dx) + i*omega*phibar_prime(j) );
    else % Trailing edge
      pbar_prime(j) = -rhoref*( Uinf*((phibar_prime(j) - phibar_prime(j-1))/dx) + i*omega*phibar_prime(j) );
    end
  end

end
