
close all;
clear all;
clc;

fs = 18;
afs = 14;

% Add utilities
addpath('../utilities/2d/');
addpath('../utilities/stability/');

% Define beam attributes
L = 250;   % L = L_dim / h
Mw = 0;

E = 2.11e11; % Pa
nu = 0.3;
rho_m = 7500; % kg

% Define dimensional attributes
h = 0.001; L_panel = L*h;

% Define flow attributes
M = 1.2;
a = 328.6; % m/s
rho = 0.91; % kg / m^3
gamma = 1.4; R = 287;
T = a^2 / (gamma*R);
p = rho*R*T;

mu = rho / rho_m;
D = E / (12*(1-nu^2)*(a^2)*rho_m);

flow.pref = p;
flow.Tref = T;
flow.rhoref = rho;
flow.gamma = gamma;
flow.L = L_panel;

% Define trial modes
num_modes = 5;
N = 241; x = linspace(0,L,N);
modeshapes = zeros(num_modes,N);

bndry_pts = x;
col_pts = (bndry_pts(1:end-1) + bndry_pts(2:end)) ./ 2;

modeshape = @(n,x_ind,L) sin(n*pi.*x_ind ./ L);
mode_deriv = @(n,x_ind,L) (n*pi/L).*cos(n*pi.*x_ind ./ L);
for j = 1:num_modes
  modeshapes(j,:) = modeshape(j,x,L);
  mode_derivs(j,:) = mode_deriv(j,x,L);

  bem_psi(j,:) = modeshape(j,col_pts,L);
  bem_dpsi_dx(j,:) = mode_deriv(j,col_pts,L);
end

% figure;
% subplot(211);
% plot(col_pts, bem_psi(1,:), 'b', 'LineWidth', 2); grid on; hold on;
% ax = gca; ax.FontSize = afs;
% ylabel('$\psi_1 (x/h)$', 'interpreter', 'latex', 'fontsize', fs);
% subplot(212);
% plot(col_pts, bem_psi(2,:), 'b', 'LineWidth', 2); grid on; hold on;
% ax = gca; ax.FontSize = afs;
% ylabel('$\psi_2 (x/h)$', 'interpreter', 'latex', 'fontsize', fs);
% xlabel('$x / h$', 'interpreter', 'latex', 'fontsize', fs);
% sgtitle('Trial Modes to Compute Stability of Simply-Supported Beam', 'interpreter', 'latex', 'fontsize', fs);
% pause;

mu = 12e-5; D = 23.9;

%% Construct stability matrix
K = eye(num_modes);
for j = 1:num_modes
  K(j,j) = (D*(j*pi/L)^4 + (Mw^2)*(j*pi/L)^2)*(L/2);
end

% Begin stability search
chi = 0.5; eps1 = 10^(-5); eps2 = 1e-3;
M_min = 1.055; M_max = 1.6; Mn = 51;

M_list = linspace(M_min, M_max, Mn);
omega_soln = zeros(Mn, num_modes);
p_bem_orig = zeros(num_modes, N-1, Mn);
p_bem_soln = zeros(num_modes, N-1, Mn);
p_bem_loc = zeros(N-1,num_modes);

detA = zeros(1,num_modes);
for m = 1:Mn
  M = M_list(m);
  flow.M = M;
  conv = 0;

  disp(strcat(['Solving flow regime M = ', num2str(M)]));

  % Compute initial guesses
  omega = zeros(num_modes,1);
  for j = 1:num_modes
    omega(j) = sqrt(D*(j*pi/L)^4 + (Mw^2)*(j*pi/L)^2);
    %disp(j); disp(omega(j)*a/h);
  end

  % for j = 1:num_modes
  %   for k = 1:num_modes
  %     p_nondim = p_upt(mu, M, omega(k), col_pts, bndry_pts, bem_psi(j,:), bem_dpsi_dx(j,:));
  %     [phi_bem, p_bem] = bem2d(h.*col_pts, h.*bndry_pts, h.*bem_psi(j,:), bem_dpsi_dx(j,:), 2*i, -(a/h)*omega(k), flow);
  %     p_bem_nondim = p_bem ./ (rho_m*a^2);
  %     f = figure; f.Position = [560 431 704 516];
  %     subplot(211);
  %     plot(col_pts, real(p_nondim), 'r', 'LineWidth', 2); grid on; hold on;
  %     plot(col_pts, real(p_bem_nondim), 'b', 'LineWidth', 2);
  %     ax = gca; ax.FontSize = 16;
  %     legend('UPT', 'BEM', 'interpreter', 'latex', 'fontsize', fs, 'location', 'southeast');
  %     ylabel('Re( $p^{\prime} / (\rho_m a^2)$ )', 'interpreter', 'latex', 'fontsize', fs);
  %     subplot(212);
  %     plot(col_pts, imag(p_nondim), 'r', 'LineWidth', 2); grid on; hold on;
  %     plot(col_pts, imag(p_bem_nondim), 'b', 'LineWidth', 2);
  %     ax = gca; ax.FontSize = 16;
  %     ylabel('Im( $p^{\prime} / (\rho_m a^2)$ )', 'interpreter', 'latex', 'fontsize', fs);
  %     xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
  %     sgtitle(strcat(['UPT (Ved.) - BEM Comparison, $\psi_', num2str(j), ', \omega_', num2str(k), ', N_e = ', num2str(N-1), '$']), 'interpreter', 'latex', 'fontsize', fs);
  %   end
  % end
  % pause;

  %omega = [0.0010 + 0.0006i   0.0028 - 0.0005i   0.0067 - 0.0004i   0.0122 - 0.0002i   0.0192 - 0.0002i   0.0277 - 0.0002i   0.0378 - 0.0001i]';
  % omega_ved = dlmread('./stability_data/vedeneev_eigenfreq_solns.csv');
  % omega = omega_ved(Mn,:); omega = omega';

  iter_num = 1;
  while (conv < 1)
  
    omegap = zeros(num_modes, 1);
    for j = 1:num_modes

      % Pull current iteration
      omega0 = omega(j);
  
      % Compute pressure matrix with current omega
      P = zeros(num_modes, num_modes);
      for l = 1:num_modes
        % p_nondim = p_upt(mu, M, omega0, col_pts, bndry_pts, bem_psi(l,:), bem_dpsi_dx(l,:));

        flow.L = L*h;
        [phi_bem, p_bem] = bem2d(h.*col_pts, h.*bndry_pts, h.*bem_psi(l,:), bem_dpsi_dx(l,:), 2*i, -(a/h)*omega0, flow);
        p_bem = (1/(rho_m*a^2)).*p_bem; if (l==j) p_bem_loc(:,j) = p_bem; end

        if ((iter_num == 1) && (l==j))
          p_bem_orig(j,:,m) = p_bem;
        end

        % if (l == 2)
        % figure;
        % subplot(211);
        % plot(col_pts, real(p_nondim), 'b'); grid on; hold on;
        % plot(col_pts, real(p_bem), 'r');
        % ax = gca; ax.FontSize = afs;
        % legend('$p_{V}$', '$p_{BEM}$', 'interpreter', 'latex', 'fontsize', fs, 'location', 'northeast');
        % ylabel('Re( $p^{\prime} ) / \rho_m a^2$ ', 'interpreter', 'latex', 'fontsize', fs);
        % subplot(212);
        % plot(col_pts, imag(p_nondim), 'b'); grid on; hold on;
        % plot(col_pts, imag(p_bem), 'r');
        % ax = gca; ax.FontSize = afs;
        % ylabel('Im( $p^{\prime} ) / \rho_m a^2$', 'interpreter', 'latex', 'fontsize', fs);
        % xlabel('$x / h$', 'interpreter', 'latex', 'fontsize', fs);
        % sgtitle(strcat(['Comparison of BEM and Vedeneev Unsteady Potential Flow, $M = ' num2str(M), '$']), 'interpreter', 'latex', 'fontsize', fs);
        % pause;
        % end

        for k = 1:num_modes
          %pjk = gauss_1pt_quad(bndry_pts, p_nondim, bem_psi(k,:));
          pjk = gauss_1pt_quad(bndry_pts, p_bem, bem_psi(k,:));
          P(l,k) = pjk;
        end
      end
      % Compute stability matrix
      Lmat = (L*(omega0^2)/2).*eye(num_modes);
      A = K + P - Lmat;
      detA(j) = abs(det(A)); detA_hist(iter_num,j) = detA(j);
      if (iter_num == 1)
        detA_init(j) = abs(detA(j));
      end
      A = sym(A);
      
      % Change relevant entry
      syms wpp
      A(j,j) = K(j,j) + P(j,j) - (L*wpp^2)/2;
      eqn = det(A) == 0;
      B = double(solve(eqn));
      
      omegap(j) = B(real(B)>0);
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
        omega = omega_new;
      end
    else
      omega = omega_new
    end
    iter_num = iter_num + 1;

  end

  disp('Frequency solutions:');
  disp(omega);

end

figure; Np = max(size(detA_hist));
for j = 1:2
  semilogy(1:Np, abs(detA_hist(:,j))); grid on; hold on;
end
legend('$\omega_1$', '$\omega_2$', 'interpreter', 'latex', 'fontsize', fs);
xlabel('Iterations', 'interpreter', 'latex', 'fontsize', fs);
ylabel('| det A($\omega$) |', 'interpreter', 'latex', 'fontsize', fs);
title(strcat(['Determinant History, $M = ', num2str(M_min), ' N_e = ', num2str(N), '$']), 'interpreter', 'latex', 'fontsize', fs);

figure;
subplot(121);
for j = 1:2
  plot(1:Np, real(omega_hist(j,:))); grid on; hold on;
end
legend('$\omega_1$', '$\omega_2$', 'interpreter', 'latex', 'fontsize', fs, 'location', 'northeast');
xlabel('Iterations', 'interpreter', 'latex', 'fontsize', fs);
ylabel('Re( $\omega$ )', 'interpreter', 'latex', 'fontsize', fs);
subplot(122);
for j = 1:2
  plot(1:Np, imag(omega_hist(j,:))); grid on; hold on;
end
xlabel('Iterations', 'interpreter', 'latex', 'fontsize', fs);
ylabel('Im( $\omega$ )', 'interpreter', 'latex', 'fontsize', fs);
sgtitle(strcat(['Frequency History, $M = ', num2str(M_min), ' N_e = ', num2str(N), '$']), 'interpreter', 'latex', 'fontsize', fs);

check_vedeneev = 1;
if (check_vedeneev)

  % Read in Vedeneev eigenfrequencies
  eigenfreq1 = dlmread('./stability_data/eigenfreq1.csv');
  eigenfreq2 = dlmread('./stability_data/eigenfreq2.csv');
  eigenfreq3 = dlmread('./stability_data/eigenfreq3.csv');
  eigenfreq4 = dlmread('./stability_data/eigenfreq4.csv');
  jj=1:4; omega_nat = sqrt(D.*(jj.*pi./L).^4);

  f_comp = figure;
  s(1) = scatter(eigenfreq1(:,1), eigenfreq1(:,2), 'xk'); grid on; hold on;
  s(2) = scatter(eigenfreq2(:,1), eigenfreq2(:,2), 'xk');
  s(3) = scatter(eigenfreq3(:,1), eigenfreq3(:,2), 'xk');
  s(4) = scatter(eigenfreq4(:,1), eigenfreq4(:,2), 'xk');

  if (num_modes > 4)

    s(5) = plot(real(omega_soln(:,1)) .* 10^3, imag(omega_soln(:,1)) .* 10^4, '--k', 'DisplayName', '$\omega_1$');
    s(6) = plot(real(omega_soln(:,2)) .* 10^3, imag(omega_soln(:,2)) .* 10^4, '--k', 'DisplayName', '$\omega_2$');
    s(7) = plot(real(omega_soln(:,3)) .* 10^3, imag(omega_soln(:,3)) .* 10^4, '--k', 'DisplayName', '$\omega_3$');
    s(8) = plot(real(omega_soln(:,4)) .* 10^3, imag(omega_soln(:,4)) .* 10^4, '--k', 'DisplayName', '$\omega_4$');
    s(9) = scatter(real(omega_soln(:,1)) .* 10^3, imag(omega_soln(:,1)) .* 10^4, 24, M_list); 
    colorbar; colormap('jet'); caxis([M_list(1), M_list(end)]);
    s(10) = scatter(real(omega_soln(:,2)) .* 10^3, imag(omega_soln(:,2)) .* 10^4, 24, M_list); 
    s(11) = scatter(real(omega_soln(:,3)) .* 10^3, imag(omega_soln(:,3)) .* 10^4, 24, M_list); 
    s(12) = scatter(real(omega_soln(:,4)) .* 10^3, imag(omega_soln(:,4)) .* 10^4, 24, M_list); 
    s(13) = scatter(omega_nat .* 10^3, zeros(1,4), 'filled');
    % legend(s(5:8), 'interpreter', 'latex', 'fontsize', fs, 'location', 'northwest');
  else
    for j = 1:num_modes
      s(j+4) = plot(real(omega_soln(:,1)).*10^3, imag(omega_soln(:,j)).*10^4, '--k', 'DisplayName', strcat(['$\omega_', num2str(j), '$']));
      s(j+4+num_modes) = scatter(real(omega_soln(:,j)).*10^3, imag(omega_soln(:,j)).*10^4, 24, M_list);
    end
  end
  cb = colorbar; colormap('jet'); caxis([M_list(1), M_list(end)]); ylabel(cb, '$M$', 'interpreter', 'latex', 'FontSize', fs);
  ax = gca; ax.FontSize = afs;
  xlabel('Re( $\omega$ )', 'interpreter', 'latex', 'fontsize', fs);
  ylabel('Im( $\omega$ )', 'interpreter', 'latex', 'fontsize', fs);
  title(strcat(['Solution to Eigenfrequencies, $M = ', num2str(M_list(1)), '$ to $M = ', num2str(M_list(end)), '$, $N_e = ', num2str(N), '$']), 'interpreter', 'latex', 'fontsize', fs);

end

mode_num = 1; mach_ind = 1;
ret = plot_soln(col_pts, p_bem_orig, p_bem_soln, mode_num, mach_ind, M_list);


function [p] = p_upt(mu, M, omega, x, bndry_pts, W, dWdx)

  N = length(x); p = zeros(N,1);
  for j = 1:N
    x_loc = x(j);
    p(j) = (mu*M/sqrt(M^2-1))*(-i*omega*W(j) + M*dWdx(j));
    running_int = 0;
    for k = 1:j-1
      xx_loc = x(k);
      intg = (-i*omega*W(k) + M*dWdx(k)) * exp(i*M*omega*(x_loc-xx_loc)/(M^2-1)) * (i*besselj(0,-omega*(x_loc-xx_loc)/(M^2-1)) + M*besselj(1,-omega*(x_loc-xx_loc)/(M^2-1)));
      running_int = running_int + ((bndry_pts(k+1)-bndry_pts(k))/2)*2*intg;
    end
    % Account for final sliver
    intg = (-i*omega*W(end) + M*dWdx(end)) * (i*besselj(0,0) + M*besselj(1,0));
    running_int = running_int + ((bndry_pts(end) - bndry_pts(end-1))/2)*intg;
    p(j) = p(j) + (mu*omega/((M^2 - 1)^(3/2)))*running_int;
  end
  
end

