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
    for j = 1:N
      loc = col_pts(j,:); x = loc(1); y = loc(2); z = 0;
      for l = 1:N
        phibar_prime(j) = phibar_prime(j) - dphi_dn(l)*gq3d_subsonic(x,y,z,bndry_pts(4*(l-1)+1:4*(l-1)+4,:),k,M); % Compute Green's function via Gaussian quadrature
      end
      phibar_prime(j) = (1/(4*pi)) * phibar_prime(j);
    end
  else  % Utilize supersonic Green's function -- Jones

    % Scale dphi_dn
    mu = asin(1 / M);
    lambda = omega * La / Uinf;

    % Perform transformation
    % dphi_dn = La .* exp(i*omega(1/cos(mu)^2) .* col_pts(:,1)) .* dphi_dn;
    bndry_pts = [bndry_pts(:,1) .* tan(mu) ./ La, bndry_pts(:,2) ./ La]; % convert to nondim
    col_pts = [col_pts(:,1) .* tan(mu) ./ La, col_pts(:,2) ./ La]; % convert to nondim
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* col_pts(:,1) ./ Uinf) .* dphi_dn; % convert to nondim
    %dphi_dn = exp(i*omega*(La/(cos(mu)*sin(mu))) .* col_pts(:,1) ./ Uinf) .* dphi_dn;

    test_pt = 11*40+38;
    phibar_prime = zeros(N,1);
    for j = 1:N
      loc = col_pts(j,:); x = loc(1); y = loc(2); z = 0;
      for l = 1:N

        % Identify element type -- either in Mach-fore cone or segmented by it
        loc_bndry_pts = bndry_pts(4*(l-1)+1:4*(l-1)+4, :);
        alphas = sqrt((loc_bndry_pts(:,1) - x).^2);
        test1 = y.*ones(4,1) + alphas; test2 = y.*ones(4,1) - alphas;

        if ( (j == (Nx*floor(Ny/2) + (Nx-3))) && (l==N) )
          figure;
          scatter(bndry_pts(:,1), bndry_pts(:,2), 12, 'b'); hold on;
          scatter(x,y,12,'r','filled');
          xt = linspace(0,x,101); ytu = y + (x-xt); ytl = y - (x-xt);
          plot(xt, ytu, 'r'); plot(xt,ytl,'r');
          pause;
        end

        if (l==j) % Contribution from own element
          if (loc_bndry_pts(1,2) < test2(1))
            phibar_prime(j) = phibar_prime(j) - dphi_dn(l)*kernel6(x,y,loc_bndry_pts,lambda,mu);
          else
            phibar_prime(j) = phibar_prime(j) - dphi_dn(l)*kernel7(x,y,loc_bndry_pts,lambda,mu);
          end

        elseif ( (x == col_pts(l,1)) && ( ((loc_bndry_pts(1,2) < test1(1)) && (col_pts(l,2) > y)) || ((loc_bndry_pts(3,2) > test2(3)) && (col_pts(l,2) < y)) || (col_pts(l,2) == y)  ))
          if ( ((loc_bndry_pts(1,2) < test1(1)) && (loc_bndry_pts(3,2) > test1(3))) || ((loc_bndry_pts(3,2) > test2(3)) && (loc_bndry_pts(1,2) < test2(1))) ) 
            disp('1');
            phibar_prime(j) = phibar_prime(j) - dphi_dn(l)*kernel9(x,y,loc_bndry_pts,lambda,mu);
          else
            disp('2');
            phibar_prime(j) = phibar_prime(j) - dphi_dn(l)*kernel8(x,y,loc_bndry_pts,lambda,mu);
          end

        elseif ( (x >= col_pts(l,1))  && ( ((loc_bndry_pts(1,2) < test1(1)) && (col_pts(l,2) > y)) || ((loc_bndry_pts(3,2) > test2(3)) && (col_pts(l,2) < y)) || (col_pts(l,2) == y)  )) % Element only contibutes if it is behind collocation point
          if ((loc_bndry_pts(1,2) < test2(1)) && (loc_bndry_pts(3,2) > test1(3))) % Element contains fore-cone
            disp('3');
            phibar_prime(j) = phibar_prime(j) - dphi_dn(l)*kernel4(x,y,loc_bndry_pts,lambda,mu);
          elseif ((loc_bndry_pts(2,2) < test2(2)) && loc_bndry_pts(4,2) > test1(4)) % Scenario 5
            phibar_prime(j) = phibar_prime(j) - dphi_dn(l)*kernel5(x,y,loc_bndry_pts,lambda,mu);
          elseif ( ((loc_bndry_pts(1,2) < test1(1)) && (col_pts(l,2) > y)) || ((loc_bndry_pts(3,2) > test2(3)) && (col_pts(l,2) < y)) )
            if ((loc_bndry_pts(2,2) > test1(2)) || (loc_bndry_pts(4,2) < test2(4)))
              % Scenario 1: One boundary point contained in Mach cone
              if ((loc_bndry_pts(3,2) > test1(3)) || (loc_bndry_pts(1,2) < test2(1)))
                phibar_prime(j) = phibar_prime(j) - dphi_dn(l)*kernel1(x,y,loc_bndry_pts,lambda,mu);
              else
                disp('6');
                phibar_prime(j) = phibar_prime(j) - dphi_dn(l)*kernel12(x,y,loc_bndry_pts,lambda,mu);
              end
            elseif ((loc_bndry_pts(3,2) > test1(3)) || (loc_bndry_pts(1,2) < test2(1)))
              % Scenario 2: Two boundary points contained in Mach cone
              phibar_prime(j) = phibar_prime(j) - dphi_dn(l)*kernel2(x,y,loc_bndry_pts,lambda,mu);
            elseif ((loc_bndry_pts(4,2) > test1(4)) || (loc_bndry_pts(2,2) < test2(2)))
              % Scenario 3: Three boundary points contained in Mach cone
              phibar_prime(j) = phibar_prime(j) - dphi_dn(l)*kernel3(x,y,loc_bndry_pts,lambda,mu);
            else
              % Otherwise: Entire element in Mach cone or out of it
              phibar_prime(j) = phibar_prime(j) - dphi_dn(l)*kernel(x,y,loc_bndry_pts,lambda,mu);
            end
          end
        end

      end % end loop
    end  % end loop
    phibar_prime = (1/pi) .* phibar_prime .* exp(-i * lambda * (1/(cos(mu)*sin(mu))) .* col_pts(:,1)); % physical coords
  end

  % Compute pressure fluctuation via finite differences -- physical pbar
  for l = 1:Nx
    for m = 1:Ny
      j = (m-1)*Nx + l;
      if ((l>1) && (l<Nx)) % interior point in x-directional sense
        %pbar_prime(j) = -rhoref*( Uinf*((phibar_prime(j+1) - phibar_prime(j-1))/(2*dx)) + i*omega*phibar_prime(j) );
        pbar_prime(j) = -rhoref*( Uinf*((phibar_prime(j) - phibar_prime(j-1))/(dx)) + i*omega*phibar_prime(j) );
      elseif (l==1) % Leading edge
        pbar_prime(j) = -rhoref*( Uinf*((phibar_prime(j+1) - phibar_prime(j))/dx) + i*omega*phibar_prime(j) );
      else % Training edge
        pbar_prime(j) = -rhoref*( Uinf*((phibar_prime(j) - phibar_prime(j-1))/dx) + i*omega*phibar_prime(j) );
      end
    end
  end

end
