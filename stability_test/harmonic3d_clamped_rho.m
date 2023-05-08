
close all;
clear all;
clc;

fs = 18;
afs = 14;
parM = 4;

% Add utilities
addpath('../utilities/3d/');
addpath('../utilities/stability/');

% Define beam attributes
h = 0.00101854;
Lx_phys = 0.2286;
Ly_phys = 2*Lx_phys;
Lx = Lx_phys / h;
Ly = Ly_phys / h;

flow.La = Lx_phys; flow.Lb = Ly_phys;

E = 38610640842;
nu = 0.35;
rho_m = 1762.035;
Dw = E*h^3 / (12*(1-nu^2));

rhoref = 0.497743562433333; flow.rhoref = rhoref;
flow.gamma = 1.4; R = 287;
M = 1.3; lambda = 217;

Uinf = sqrt(lambda*Dw/((Lx_phys^3)*rhoref));
aref = Uinf / M;
Tref = aref^2 / (flow.gamma * R);
D = Dw / ((aref^2) * rho_m * h^3);

% Define trial modes
num_modes_x = 2; num_modes_y = 2; num_modes = num_modes_x * num_modes_y;
nmx = num2str(num_modes_x); nmy = num2str(num_modes_y);
Nx = 31; Ny = 62; N = (Nx-1)*(Ny-1);
x = linspace(0,Lx,Nx); y = linspace(0,Ly,Ny);
col_pts = zeros(N, 2);  bndry_pts = zeros(4*N, 2);
bnd_x = linspace(0,Lx,Nx); bnd_y = linspace(0,Ly,Ny);
col_x = (bnd_x(1:end-1) + bnd_x(2:end)) ./ 2;
col_y = (bnd_y(1:end-1) + bnd_y(2:end)) ./ 2;
for k = 1:Ny-1
  for j = 1:Nx-1
    m = (k-1)*(Nx-1) + j;
    col_pts(m,:) = [col_x(j), col_y(k)];
    bndry_pts(4*(m-1)+1, :) = [bnd_x(j), bnd_y(k)];
    bndry_pts(4*(m-1)+2, :) = [bnd_x(j+1), bnd_y(k)];
    bndry_pts(4*(m-1)+3, :) = [bnd_x(j), bnd_y(k+1)];
    bndry_pts(4*(m-1)+4, :) = [bnd_x(j+1), bnd_y(k+1)];
  end
end

chis = [4.73, 7.859, pi*(2*3+1)/2, pi*(2*4+1)/2, pi*(2*5+1)/2, pi*(2*6+1)/2];
modeshape = @(n,x_in,L) (-1/sqrt(2)) * (cos(chis(n).*x_in./L) - cosh(chis(n).*x_in./L) - ((cos(chis(n))-cosh(chis(n)))/(sin(chis(n)) - sinh(chis(n))))*(sin(chis(n).*x_in./L) - sinh(chis(n).*x_in./L)) );
mode_deriv = @(n,x_in,L) (-1/sqrt(2)) * (chis(n)./L) .* (-sin(chis(n).*x_in./L) - sinh(chis(n).*x_in./L) - ((cos(chis(n))-cosh(chis(n)))/(sin(chis(n)) - sinh(chis(n))))*(cos(chis(n).*x_in./L) - cosh(chis(n).*x_in./L)) );


for k = 1:num_modes_y
  for j = 1:num_modes_x

    % Construct mode shapes
    ind = (k-1)*num_modes_x + j;
    psi = modeshape(j,col_x,Lx)' * modeshape(k,col_y,Ly);
    psi = reshape(psi,1,N); bem_psi(ind,:) = psi;
    dpsi_dx = mode_deriv(j,col_x,Lx)' * modeshape(k,col_y,Ly);
    dpsi_dx = reshape(dpsi_dx,1,N); bem_dpsi_dx(ind,:) = dpsi_dx;

    % Construct stiffness matrix
    K(ind,ind) = (Lx*Ly*D/4) * ((chis(j)/Lx)^4 + 2*((chis(j)/Lx)^2)*((chis(k)/Ly)^2) + (chis(k)/Ly)^4);
    omega_nat(ind) = sqrt(K(ind,ind)*(4/(Lx*Ly)));

    % figure;
    % subplot(121);
    % scatter(col_pts(:,1), col_pts(:,2), 24, bem_psi(ind,:)); grid on; hold on;
    % colorbar; colormap('jet'); ax = gca; ax.FontSize = afs;
    % ylabel('$y / h$', 'interpreter', 'latex', 'fontsize', fs);
    % xlabel('$x / h$', 'interpreter', 'latex', 'fontsize', fs);
    % title('$\psi (x/h, y/h)$', 'interpreter', 'latex', 'fontsize', fs);
    % subplot(122);
    % scatter(col_pts(:,1), col_pts(:,2), 24, bem_dpsi_dx(ind,:)); grid on; hold on;
    % colorbar; colormap('jet'); ax = gca; ax.FontSize = afs;
    % ylabel('$y / h$', 'interpreter', 'latex', 'fontsize', fs);
    % xlabel('$x / h$', 'interpreter', 'latex', 'fontsize', fs);
    % title('$\partial \psi / \partial x(x/h, y/h)$', 'interpreter', 'latex', 'fontsize', fs);
    % sgtitle('Trial Functions for Clamped Panel', 'interpreter', 'latex', 'fontsize', fs);
    % pause;

  end
end

disp(K);
disp(omega_nat.*(aref/(2*pi*h)));

% Rearrange to march format of col_pts and bndry_pts
bem_psi = bem_psi'; bem_dpsi_dx = bem_dpsi_dx';

%% Begin stability search
chi = 0.5; eps1 = 10^(-3); eps2 = 1e-2;
M = 1.3; flow.M = M; Tinf = 227.11111; flow.Tref = Tinf;
lambda_min = 100; lambda_max = 200; lnum = 4;

lambda_list = linspace(lambda_min, lambda_max, lnum);
omega_soln = zeros(lnum, num_modes);
p_bem_orig = zeros(num_modes, (Nx-1)*(Ny-1), lnum);
p_bem_soln = zeros(num_modes, (Nx-1)*(Ny-1), lnum);
p_bem_loc = zeros((Nx-1)*(Ny-1), num_modes);

detA = zeros(1,num_modes);
for m = 1:lnum
  lambda = lambda_list(m);

  conv = 0;

  disp(strcat(['Solving flow regime M = ', num2str(M), ', lambda = ', num2str(lambda)]));

  % Compute flow quantities keeping lambda fixed
  aref = sqrt(flow.gamma * R * Tinf);
  Uinf = M*aref;
  rho = lambda * Dw / ((Lx_phys^3) * Uinf^2); flow.rhoref = rho;
  pref = rho * R * Tinf;
  D = Dw / ((aref^2) * rho_m * h^3);
  flow.pref = pref;

  % Recompute non-dim stiffness matrix based on flow quantities
  omega = zeros(num_modes, 1);
  for k = 1:num_modes_y
    for j = 1:num_modes_x
      ind = (k-1)*num_modes_x + j;
      K(ind,ind) = (Lx*Ly*D/4) * ((chis(j)/Lx)^4 + 2*((chis(j)/Lx)^2)*((chis(k)/Ly)^2) + (chis(k)/Ly)^4);
      omega(ind) = sqrt(K(ind,ind)*(4/(Lx*Ly)));
    end
  end

  iter_num = 1;
  while (conv < 1)
  
    omegap = zeros(num_modes, 1);
    for j = 1:num_modes

      % Pull current iteration
      omega0 = omega(j);
  
      % Compute pressure matrix with current omega
      P = zeros(num_modes, num_modes);
      % parfor (l = 1:num_modes, parM)
      for l = 1:num_modes
        jy = floor((l-1)/num_modes_x) + 1; jx = l - num_modes_x*(jy-1);
        flow.psi = @(xin,yin,Lx_in,Ly_in) h .* modeshape(jx,xin,Lx_in) .* modeshape(jy,yin,Ly_in);
        flow.dpsi = @(xin,yin,Lx_in,Ly_in) h .* mode_deriv(jx,xin,Lx_in) .* modeshape(jy,yin,Ly_in);

        [phi_bem, p_bem] = bem3d(h.*col_pts, h.*bndry_pts, Nx-1, h.*bem_psi(:,l), bem_dpsi_dx(:,l), 2*i, -(aref/h)*omega0, flow);
        p_bem = (1/(rho_m*aref^2)).*p_bem; if (l==j) p_bem_loc(:,j) = p_bem; end

        if ((iter_num == 1) && (l==j))
          p_bem_orig(j,:,m) = p_bem;
        end

        %figure;
        %subplot(121);
        %scatter(col_pts(:,1), col_pts(:,2), 24, (rho_m*aref^2).*real(p_bem));
        %colorbar; colormap('jet');
        %ax = gca; ax.FontSize = afs;
        %xlabel('$x / h$', 'interpreter', 'latex', 'fontsize', fs);
        %ylabel('$y / h$', 'interpreter', 'latex', 'fontsize', fs);
        %title('Re( $p^{\prime} ) / \rho_m a^2$ ', 'interpreter', 'latex', 'fontsize', fs);
        %subplot(122);
        %scatter(col_pts(:,1), col_pts(:,2), 24, (rho_m*aref^2).*imag(p_bem));
        %colorbar; colormap('jet');
        %ax = gca; ax.FontSize = afs;
        %xlabel('$x / h$', 'interpreter', 'latex', 'fontsize', fs);
        %ylabel('$y / h$', 'interpreter', 'latex', 'fontsize', fs);
        %title('Im( $p^{\prime} ) / \rho_m a^2$', 'interpreter', 'latex', 'fontsize', fs);
        %sgtitle(strcat(['BEM Solution, $M = ' num2str(M), '$']), 'interpreter', 'latex', 'fontsize', fs);

        for k = 1:num_modes
          pjk = gauss_1pt_quad_3d(bndry_pts, p_bem.', bem_psi(:,k));
          P(l,k) = pjk;
        end

      end
      % Compute stability matrix
      Lmat = (Lx*Ly*(omega0^2)/4).*eye(num_modes);
      A = K + P - Lmat;
      detA(j) = abs(det(A)); detA_hist(iter_num,j) = detA(j);
      if (iter_num == 1)
        detA_init(j) = abs(detA(j));
      end
      A = sym(A);
      
      % Change relevant entry
      syms wpp
      A(j,j) = K(j,j) + P(j,j) - (Lx*Ly*wpp^2)/4;
      eqn = det(A) == 0;
      B = double(solve(eqn));
      
      omegap(j) = B(real(B)>0);
      disp('freq solved, waiting.');
    end
  
    % Check convergence
    omega_new = (1-chi).*omegap + chi.*omega; omega_hist(:,iter_num) = omega_new;
    diff = abs((omega_new - omega) ./ omega)
    if (max(diff) < eps1)
      if (max(abs(detA ./ detA_init)) < eps2)
        omega = omega_new;
        conv = 1;
        omega_soln(m,:) = omega;
        for mm = 1:num_modes
          p_bem_soln(mm,:,m) = p_bem_loc(:,mm);
        end
      else
        omega = omega_new
      end
    else
      omega = omega_new
    end
    iter_num = iter_num + 1;

  end

  disp('Frequency solutions:');
  disp(omega);

end

if (lnum == 1)

  figure; Ns = size(detA_hist); Np = Ns(1);
  for j = 1:4
    semilogy(1:Np, abs(detA_hist(:,j))); grid on; hold on;
  end
  legend('$\omega_1$', '$\omega_2$', '$\omega_3$', '$\omega_4$', 'interpreter', 'latex', 'fontsize', fs);
  xlabel('Iterations', 'interpreter', 'latex', 'fontsize', fs);
  ylabel('| det A($\omega$) |', 'interpreter', 'latex', 'fontsize', fs);
  title(strcat(['Determinant History, $M = ', num2str(M_min), ' N_x = ', num2str(Nx), ', N_y = ', num2str(Ny), ', L_x = ', num2str(Lx), ', L_y = ', num2str(Ly), '$']), 'interpreter', 'latex', 'fontsize', fs);
  
  figure;
  subplot(121);
  for j = 1:4
    plot(1:Np, real(omega_hist(j,:))); grid on; hold on;
  end
  legend('$\omega_1$', '$\omega_2$', '$\omega_3$', '$\omega_4$', 'interpreter', 'latex', 'fontsize', fs, 'location', 'northeast');
  xlabel('Iterations', 'interpreter', 'latex', 'fontsize', fs);
  ylabel('Re( $\omega$ )', 'interpreter', 'latex', 'fontsize', fs);
  subplot(122);
  for j = 1:4
    plot(1:Np, imag(omega_hist(j,:))); grid on; hold on;
  end
  xlabel('Iterations', 'interpreter', 'latex', 'fontsize', fs);
  ylabel('Im( $\omega$ )', 'interpreter', 'latex', 'fontsize', fs);
  sgtitle(strcat(['Frequency History, $N_x = ', num2str(Nx), ', N_y = ', num2str(Ny), ', L_x = ', num2str(Lx), ', L_y = ', num2str(Ly), ', N_{mx} = ', num2str(num_modes_x), ', N_{my} = ', num2str(num_modes_y), '$']), 'interpreter', 'latex', 'fontsize', fs);

end

mode_num = 1; mach_ind = 1;
%ret = plot_soln_3d(col_pts, p_bem_orig, p_bem_soln, mode_num, mach_ind, M_list);

figure;
for j = 1:num_modes_x
  sp = str2num(['2', nmx, num2str(j)]);
  subplot(sp);
  plot(lambda_list, real(omega_soln(:,j)), '-or'); grid on;
  ax = gca; ax.FontSize = afs;
  if (j==1) ylabel('Re( $\omega$ )', 'interpreter', 'latex', 'fontsize', fs); end
  title(strcat(['$\omega_', num2str(j), '$']), 'interpreter', 'latex', 'fontsize', fs);

  sp = str2num(['2', nmx, num2str(j+num_modes_x)]);
  subplot(sp);
  plot(lambda_list, imag(omega_soln(:,j)), '-or'); grid on;
  ax = gca; ax.FontSize = afs;
  if (j==1) ylabel('Im( $\omega$ )', 'interpreter', 'latex', 'fontsize', fs); end
  xlabel('$\lambda$', 'interpreter', 'latex', 'fontsize', fs);
end
sgtitle(strcat(['Eigenfrequencies, $N_x = ', num2str(Nx), ', N_y = ', num2str(Ny), ', L_x = ', num2str(Lx), ', L_y = ', num2str(Ly), '$']), 'interpreter', 'latex', 'fontsize', fs);

