function [phibar_prime, pbar_prime] = bem3d(col_pts, bndry_pts, Nx, psi, dpsi_dx, amp, omega, flow)

  % Utilize boundary element method to obtain spatial solution to reduced velocity perturbation potential
  % and reduced pressure fluctuation for three-dimensional flow.
  % This method assumes an oscillating panel exposed to an external freestream.

  % Inputs:
  %   col_pts: collocation points at which to solve for spatial quantities, size Nx2
  %   bndry_pts: boundary element demarcating points (col_pts assumed to be at midponts, size 4*Nx2
  %   Nx: number of points in x-directional discretization
  %   psi: mode shape of harmonic deformation of beam, size Nx1
  %   dpsi_dx: first derivative of mode shape psi, size Nx2
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
  La = flow.La;
  Lb = flow.Lb;
  psi_func = flow.psi;
  dpsi_func = flow.dpsi;

  % Compute dependent flow quantities
  R = pref / (rhoref*Tref);
  aref = sqrt(gamma*R*Tref);
  Uref = sqrt(pref / rhoref);
  Uinf = M*aref;

  % Compute quantities of interest to convected wave equation
  k = omega / aref;      % Wavenumber
  beta = sqrt(1 - M^2);  % Prandtl-Glauert factor

  % Compute quantities of interest to boundary element method
  N = max(size(col_pts)); Ny = N / Nx;
  dx = bndry_pts(2,1) - bndry_pts(1,1);
  dy = bndry_pts(3,2) - bndry_pts(1,2);

  % Compute panel boundary condition from flow and deformation quantities -- no-penetration
  phibar_prime = zeros(1,N);
  %dphi_dn = (Uinf*amp*dpsi_dx + amp*i*omega*psi); % original method
  dphi_dn = (-(i/2)*amp*Uinf*dpsi_dx + (1/2)*omega*amp*psi); % computed from complex-conjugate method

  % Compute reduced velocity perturbation potential
  if (flow.M < 1)  % Utilize subsonic Green's function -- Wu and Lee
    pts = [-1 / sqrt(3), 1 / sqrt(3)];
    for ind = 1:N
      x = col_pts(ind,1); y = col_pts(ind,2);
      for j = 1:N
        loc_bndry_pts = bndry_pts(4*(j-1)+1:4*(j-1)+4,:);
        a = loc_bndry_pts(1,1); b = loc_bndry_pts(2,1); c = loc_bndry_pts(1,2); d = loc_bndry_pts(3,2);
        for l = 1:2
          for m = 1:2
            int_point_loc = [((b-a)/2)*pts(l) + ((a+b)/2), ((d-c)/2)*pts(m) + ((d+c)/2)];
            [G] = greens3d_subsonic(x,y,0, int_point_loc(1), int_point_loc(2), 0, k, M);
            psi_loc = psi_func(int_point_loc(1), int_point_loc(2), La, Lb);
            dpsi_dx_loc = dpsi_func(int_point_loc(1), int_point_loc(2), La, Lb);
            dphi_dn = -(i/2)*amp*Uinf*dpsi_dx_loc + (1/2)*amp*omega*psi_loc;
            phibar_prime(ind) = phibar_prime(ind) - dphi_dn*G*((b-a)/2)*((d-c)/2);
          end
        end
      end
    end
    phibar_prime = (1/(4*pi)) .* phibar_prime;

  else  % Utilize supersonic Green's function -- Jones

    % Scale dphi_dn
    mu = asin(1 / M);
    lambda = omega * La / Uinf;

    % Perform transformation
    % dphi_dn = La .* exp(i*omega(1/cos(mu)^2) .* col_pts(:,1)) .* dphi_dn;
    bndry_pts = [bndry_pts(:,1) .* tan(mu) ./ La, bndry_pts(:,2) ./ La]; % convert to nondim
    col_pts = [col_pts(:,1) .* tan(mu) ./ La, col_pts(:,2) ./ La]; % convert to nondim
    %dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* col_pts(:,1) ./ Uinf) .* dphi_dn; % convert to nondim
    %dphi_dn = exp(i*omega*(La/(cos(mu)*sin(mu))) .* col_pts(:,1) ./ Uinf) .* dphi_dn;

    Ny = max(size(col_pts)) / Nx;
    phibar_prime = zeros(N,1);
    for ind = 1:N
      x = col_pts(ind,1); y = col_pts(ind,2);

      for j = 1:N
        loc_bndry_pts = bndry_pts(4*(j-1)+1:4*(j-1)+4,:);
        alphas = sqrt((loc_bndry_pts(:,1) - x).^2);
        test1 = y.*ones(4,1) + alphas; test2 = y.*ones(4,1) - alphas;

        %if ((ind == (floor(Ny/2)*Nx + Nx - 3)) && (j==1))
        %  disp(ind); disp(floor(Ny/2)*Nx);
        %  figure;
        %  scatter(bndry_pts(:,1), bndry_pts(:,2), 24, 'b'); hold on;
        %  scatter(col_pts(:,1), col_pts(:,2), 24, 'r');
        %  scatter(col_pts(ind,1), col_pts(ind,2), 24, 'r', 'filled');
        %  xt = linspace(0,col_pts(ind,1), 101); alphas = sqrt((xt-x).^2);
        %  ytu = y + alphas; ytl = y - alphas;
        %  plot(xt,ytu,'r'); plot(xt,ytl,'r');
        %end

        % Check for scenarios
        if (j == ind) % Contribution from own element
          if (loc_bndry_pts(1,2) < test2(1))
            contribution = kernel6(x,y,loc_bndry_pts,lambda,mu,flow,amp,omega,Uinf);
            phibar_prime(ind) = phibar_prime(ind) - contribution;
          else
            contribution = kernel7(x,y,loc_bndry_pts,lambda,mu,flow,amp,omega,Uinf);
            phibar_prime(ind) = phibar_prime(ind) - contribution;
          end
        elseif ( (x == col_pts(j,1)) && ( ((loc_bndry_pts(1,2) < test1(1)) && (col_pts(j,2) > y)) || ((loc_bndry_pts(3,2) > test2(3)) && (col_pts(j,2) < y)) || (col_pts(j,2) == y)  ))
          if ( ((loc_bndry_pts(1,2) < test1(1)) && (loc_bndry_pts(3,2) > test1(3))) || ((loc_bndry_pts(3,2) > test2(3)) && (loc_bndry_pts(1,2) < test2(1))) ) 
            contribution = kernel9(x,y,loc_bndry_pts,lambda,mu,flow,amp,omega,Uinf);
            phibar_prime(ind) = phibar_prime(ind) - contribution;
          else
            contribution = kernel8(x,y,loc_bndry_pts,lambda,mu,flow,amp,omega,Uinf);
            phibar_prime(ind) = phibar_prime(ind) - contribution;
          end
        elseif ( (x >= col_pts(j,1))  && ( ((loc_bndry_pts(1,2) < test1(1)) && (col_pts(j,2) > y)) || ((loc_bndry_pts(3,2) > test2(3)) && (col_pts(j,2) < y)) || (col_pts(j,2) == y)  )) % Element only contibutes if it is behind collocation point
          if ((loc_bndry_pts(1,2) < test2(1)) && (loc_bndry_pts(3,2) > test1(3))) % Element contains fore-cone
            contribution = kernel4(x,y,loc_bndry_pts,lambda,mu,flow,amp,omega,Uinf);
            phibar_prime(ind) = phibar_prime(ind) - contribution;
          elseif ((loc_bndry_pts(2,2) < test2(2)) && loc_bndry_pts(4,2) > test1(4)) % Scenario 5
            contribution = kernel5(x,y,loc_bndry_pts,lambda,mu,flow,amp,omega,Uinf);
            phibar_prime(ind) = phibar_prime(ind) - contribution;
          elseif ( ((loc_bndry_pts(1,2) < test1(1)) && (col_pts(j,2) > y)) || ((loc_bndry_pts(3,2) > test2(3)) && (col_pts(j,2) < y)) )
            if ((loc_bndry_pts(2,2) > test1(2)) || (loc_bndry_pts(4,2) < test2(4)))
              % Scenario 1: One boundary point contained in Mach cone
              if ((loc_bndry_pts(3,2) > test1(3)) || (loc_bndry_pts(1,2) < test2(1)))
                contribution = kernel1(x,y,loc_bndry_pts,lambda,mu,flow,amp,omega,Uinf);
                phibar_prime(ind) = phibar_prime(ind) - contribution;
              else
                contribution = kernel12(x,y,loc_bndry_pts,lambda,mu,flow,amp,omega,Uinf);
                phibar_prime(ind) = phibar_prime(ind) - contribution;
              end
            elseif ((loc_bndry_pts(3,2) > test1(3)) || (loc_bndry_pts(1,2) < test2(1)))
              % Scenario 2: Two boundary points contained in Mach cone
              contribution = kernel2(x,y,loc_bndry_pts,lambda,mu,flow,amp,omega,Uinf);
              phibar_prime(ind) = phibar_prime(ind) - contribution;
            elseif ((loc_bndry_pts(4,2) > test1(4)) || (loc_bndry_pts(2,2) < test2(2)))
              % Scenario 3: Three boundary points contained in Mach cone
              contribution = kernel3(x,y,loc_bndry_pts,lambda,mu,flow,amp,omega,Uinf);
              phibar_prime(ind) = phibar_prime(ind) - contribution;
            else
              % Otherwise: Entire element in Mach cone or out of it
              %contribution = kernel(x,y,loc_bndry_pts,lambda,mu,flow,amp,omega,Uinf);
              contribution = kernel_hires(x,y,loc_bndry_pts,lambda,mu,flow,amp,omega,Uinf);
              phibar_prime(ind) = phibar_prime(ind) - contribution;
            end
          end
        end

      end
    end
    %phibar_prime = (1/pi) .* phibar_prime .* exp(-i * lambda * (1/(cos(mu)*sin(mu))) .* col_pts(:,1)); % physical coords

    phibar_prime = (1/pi) .* phibar_prime .* exp(-i*omega*(La/(cos(mu)*sin(mu))) .* col_pts(:,1) ./ Uinf);

  end

  % Compute pressure fluctuation via finite differences -- physical pbar

  if (M< 1)

    for l = 1:Nx
      for m = 1:Ny
        j = (m-1)*Nx + l;
        if ((l>1) && (l<Nx)) % interior point in x-directional sense
          pbar_prime(j) = -rhoref*( Uinf*((phibar_prime(j+1) - phibar_prime(j-1))/(2*dx)) + i*omega*phibar_prime(j) );
        elseif (l==1) % Leading edge
          pbar_prime(j) = -rhoref*( Uinf*((phibar_prime(j+1) - phibar_prime(j))/dx) + i*omega*phibar_prime(j) );
        else % Training edge
          pbar_prime(j) = -rhoref*( Uinf*((phibar_prime(j) - phibar_prime(j-1))/dx) + i*omega*phibar_prime(j) );
        end
      end
    end

  else

    %for l = 1:Nx
    %  for m = 1:Ny
    %    j = (m-1)*Nx + l;
    %    if ((l>1) && (l<Nx)) % interior point in x-directional sense
    %      pbar_prime(j) = -rhoref*( Uinf*((phibar_prime(j+1) - phibar_prime(j-1))/(2*dx)) + i*omega*phibar_prime(j) );
    %      %pbar_prime(j) = -rhoref*( Uinf*((phibar_prime(j) - phibar_prime(j-1))/(dx)) + i*omega*phibar_prime(j) );
    %    elseif (l==1) % Leading edge
    %      pbar_prime(j) = -rhoref*( Uinf*((phibar_prime(j+1) - phibar_prime(j))/dx) + i*omega*phibar_prime(j) );
    %    else % Training edge
    %      pbar_prime(j) = -rhoref*( Uinf*((phibar_prime(j) - phibar_prime(j-1))/dx) + i*omega*phibar_prime(j) );
    %    end
    %  end
    %end

    for l = 1:Nx
      for m = 1:Ny
        j = (m-1)*Nx + l;
        if ((l>2) && (l<(Nx-1)))
          pbar_prime(j) = -rhoref * (Uinf*((-phibar_prime(j+2) + 8*phibar_prime(j+1) - 8*phibar_prime(j-1) + phibar_prime(j-2))/(12*dx)) + i*omega*phibar_prime(j));
        elseif ((l==2) || (l==(Nx-1)))
          pbar_prime(j) = -rhoref * (Uinf*((phibar_prime(j+1) - phibar_prime(j-1))/(2*dx)) + i*omega*phibar_prime(j) );
        elseif (l==1)
          pbar_prime(j) = -rhoref*( Uinf*((phibar_prime(j+1) - phibar_prime(j))/dx) + i*omega*phibar_prime(j) );
        else % Training edge
          pbar_prime(j) = -rhoref*( Uinf*((phibar_prime(j) - phibar_prime(j-1))/dx) + i*omega*phibar_prime(j) );
        end
      end
    end

  end

end
