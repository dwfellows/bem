
close all;
clear all;
clc;

fs = 24;
afs = 20;

% Add utilities
addpath('../utilities/2d/');
addpath('../utilities/stability/');

% Define beam attributes
L = 300;   % L = L_dim / h
Mw = 0;

E = 2.11e11; % Pa
nu = 0.3;
rho_m = 7500; % kg

% Define dimensional attributes
h = 0.001; L_panel = L*h;

% Define flow attributes
a = 328.6; % m/s
rho = 0.91;  % kg / m^3
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
num_modes = 2;
N = 501; x_bem = linspace(0,L,N);
modeshapes = zeros(num_modes,N);

bndry_pts = x_bem;
col_pts = (bndry_pts(1:end-1) + bndry_pts(2:end)) ./ 2;

% Define DMD quantities
x_dmd = linspace(0,L,61);
bndry_pts_dmd = x_dmd;
col_pts_dmd = (bndry_pts_dmd(1:end-1) + bndry_pts_dmd(2:end)) ./ 2;

chis = [4.73, 7.859, pi*(2*3+1)/2, pi*(2*4+1)/2, pi*(2*5+1)/2, pi*(2*6+1)/2];
modeshape = @(n,x_in,L) (-1/sqrt(2)) * (cos(chis(n).*x_in./L) - cosh(chis(n).*x_in./L) - ((cos(chis(n))-cosh(chis(n)))/(sin(chis(n)) - sinh(chis(n))))*(sin(chis(n).*x_in./L) - sinh(chis(n).*x_in./L)) );
mode_deriv = @(n,x_in,L) (-1/sqrt(2)) * (chis(n)./L) .* (-sin(chis(n).*x_in./L) - sinh(chis(n).*x_in./L) - ((cos(chis(n))-cosh(chis(n)))/(sin(chis(n)) - sinh(chis(n))))*(cos(chis(n).*x_in./L) - cosh(chis(n).*x_in./L)) );
for j = 1:num_modes

  modeshapes(j,:) = modeshape(j,x_bem,L);
  mode_derivs(j,:) = mode_deriv(j,x_bem,L);

  bem_psi(j,:) = modeshape(j,col_pts,L);
  bem_dpsi_dx(j,:) = mode_deriv(j,col_pts,L);

  modeshapes_dmd(j,:) = modeshape(j,x_dmd,L);
  mode_derivs_dmd(j,:) = mode_deriv(j,x_dmd,L);

  bem_psi_dmd(j,:) = modeshape(j,col_pts_dmd,L);
  bem_dpsi_dx_dmd(j,:) = mode_deriv(j,col_pts_dmd,L);

end

% figure;
% subplot(211);
% plot(col_pts, bem_psi(1,:), 'b', 'LineWidth', 2); grid on; hold on;
% plot(col_pts, bem_dpsi_dx(1,:), 'r', 'LineWidth', 2);
% ax = gca; ax.FontSize = afs;
% ylabel('$\psi_1 (x/h)$', 'interpreter', 'latex', 'fontsize', fs);
% subplot(212);
% plot(col_pts, bem_psi(2,:), 'b', 'LineWidth', 2); grid on; hold on;
% plot(col_pts, bem_dpsi_dx(2,:), 'r', 'LineWidth', 2);
% ax = gca; ax.FontSize = afs;
% ylabel('$\psi_2 (x/h)$', 'interpreter', 'latex', 'fontsize', fs);
% xlabel('$x / h$', 'interpreter', 'latex', 'fontsize', fs);
% sgtitle('Trial Modes to Compute Stability of Clamped Beam', 'interpreter', 'latex', 'fontsize', fs);


%% Construct stability matrix
K = eye(num_modes);
for j = 1:num_modes
  K(j,j) = D*((chis(j)/L)^4)*(L/2);
end

for j = 1:num_modes
  omega_nat(j) = sqrt(D*(chis(j)/L)^4);
end

% Begin stability search
chi = 0.5; eps1 = 10^(-5); eps2 = 1e-3;
M_min = 1.2; M_max = 1.5; Mn = 41;

interp_on = 1;
M_list = linspace(M_min, M_max, Mn);
omega_soln = zeros(Mn, num_modes);
p_bem_orig = zeros(num_modes, N-1, Mn);
p_bem_soln = zeros(num_modes, N-1, Mn);
p_bem_loc = zeros(N-1,num_modes);

detA = zeros(1,num_modes);
for m = 1:Mn
  M = M_list(m);
  m_choice = num2str(M); m_choice = strcat([m_choice(1), m_choice(3:end)]);
  flow.M = M;
  conv = 0;

  disp(strcat(['Solving flow regime M = ', num2str(M)]));

  % Compute initial guesses
  omega = zeros(num_modes,1);
  for j = 1:num_modes
    omega(j) = sqrt(D*(chis(j)/L)^4);
    %disp(j); disp(omega(j)*a/h);
  end

  omegap = zeros(num_modes, 1);
  for j = 1:num_modes

    % Pull current iteration
    omega0 = omega(j);
  
    % Compute pressure matrix with current omega
    P = zeros(num_modes, num_modes);
    for l = 1:num_modes
     % Compute unsteady and BEM-based response
      p_uns = p_upt_func(mu, M, omega0, x_bem, bndry_pts, modeshapes(l,:), mode_derivs(l,:));

      flow.L = L*h;
      [phi_bem, p_bem] = bem2d(h.*col_pts, h.*bndry_pts, h.*bem_psi(l,:), bem_dpsi_dx(l,:), 2*i, -(a/h)*omega0, flow);
      p_bem = (1/(rho_m*a^2)).*p_bem; if (l==j) p_bem_loc(:,j) = p_bem; end

      % Compute local piston theory response
      p_lpt_nondim = p_lpt_func(mu, M, omega0, x_dmd, bndry_pts_dmd, modeshapes_dmd(l,:), mode_derivs_dmd(l,:));
      p_lpt = p_lpt_nondim .* (rho_m*a^2);

      % Compute pressure response via reading existing mode
      if (interp_on == 0)
        dmd_mode_exist = dlmread(strcat(['./2d_dmd_modes_clamped/m', m_choice, '/mode', num2str(l), '_omega', num2str(j), '/dmd_modes.csv']));
        alpha_tilde_exist = dlmread(strcat(['./2d_dmd_modes_clamped/m', m_choice, '/mode', num2str(l), '_omega', num2str(j), '/alphas_tilde.csv']));
        dmd_choice = 2;
        dmd_mode_exist = dmd_mode_exist(:,dmd_choice); alpha_tilde_exist = alpha_tilde_exist(dmd_choice);
      else
        % Compute pressure response via interpolating mode
        dmd_mode_interp_real = dlmread(strcat(['./2d_dmd_modes_clamped/interp/mode', num2str(l), '_omega', num2str(j), '_mode_real_coeffs_interp.csv']));
        dmd_mode_interp_imag = dlmread(strcat(['./2d_dmd_modes_clamped/interp/mode', num2str(l), '_omega', num2str(j), '_mode_imag_coeffs_interp.csv']));
        atilde_interp_real = dlmread(strcat(['./2d_dmd_modes_clamped/interp/mode', num2str(l), '_omega', num2str(j), '_atilde_real_coeffs_interp.csv']));
        atilde_interp_imag = dlmread(strcat(['./2d_dmd_modes_clamped/interp/mode', num2str(l), '_omega', num2str(j), '_atilde_imag_coeffs_interp.csv']));
        alpha_tilde = polyval(atilde_interp_real, M) + i*polyval(atilde_interp_imag, M);
        src = size(dmd_mode_interp_real); srcp = src(2);
        sic = size(dmd_mode_interp_imag); sicp = sic(2);
        mrc = zeros(1, srcp); mic = zeros(1, sicp);
        for ll = 1:srcp
          loc_p_coeffs = dmd_mode_interp_real(:,ll);
          loc_p = polyval(loc_p_coeffs, M);
          mrc(ll) = loc_p;
        end
        for ll = 1:sicp
          loc_p_coeffs = dmd_mode_interp_imag(:,ll);
          loc_p = polyval(loc_p_coeffs, M);
          mic(ll) = loc_p;
        end
        dmd_mode = polyval(mrc, x_dmd.*h./L_panel) + i.*polyval(mic, x_dmd.*h./L_panel);
        dmd_mode = dmd_mode.';
      end

      if (interp_on == 0)
        % Compute pressure mode using DMD
        p_corr_exist = p_lpt + p.*(-i*2*h).*dmd_mode_exist.*alpha_tilde_exist;
        p_corr_nondim_exist = p_corr_exist ./ (rho_m*a^2);
      else
        p_corr = p_lpt + p.*(-i*2*h).*dmd_mode.*alpha_tilde;
        p_corr_nondim = p_corr ./ (rho_m*a^2);
      end

      % p_corr_nondim = p_corr_nondim_exist;

      %% Plot mode comparison to UPT and LPT
      %f = figure; f.Position = [356 324 774 545];
      %subplot(211); plot(x_dmd, real(p_lpt_nondim), 'LineWidth', 2); grid on; hold on; plot(x_dmd, real(p_corr_nondim_exist), 'LineWidth', 2); plot(x_dmd, real(p_corr_nondim), 'LineWidth', 2); plot(x_bem, real(p_uns), 'LineWidth', 2);
      %ax = gca; ax.FontSize = afs;
      %ylabel('Re( $\hat{p}^{\prime} / \rho_m a^2$ )', 'interpreter', 'latex', 'fontsize', fs);
      %legend('LPT', 'DMD', 'DMD-Interp', 'UPT', 'interpreter', 'latex', 'fontsize', fs);
      %subplot(212); plot(x_dmd, imag(p_lpt_nondim), 'LineWidth', 2); grid on; hold on; plot(x_dmd, imag(p_corr_nondim_exist), 'LineWidth', 2); plot(x_dmd, imag(p_corr_nondim), 'LineWidth', 2); plot(x_bem, imag(p_uns), 'LineWidth', 2);
      %ax = gca; ax.FontSize = afs;
      %ylabel('Im( $\hat{p}^{\prime} / \rho_m a^2$ )', 'interpreter', 'latex', 'fontsize', fs);
      %xlabel('$x / h$', 'interpreter', 'latex', 'fontsize', fs);
      %sgtitle(strcat(['Comparison of Pressure Fluctuation from Theories, $\psi_', num2str(l), '$, $\omega_', num2str(j), '$, $M = ', num2str(M), '$']), 'interpreter', 'latex', 'fontsize', fs);
      %pause;

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
        if (interp_on == 1)
          pjk = gauss_1pt_quad(bndry_pts_dmd, p_corr_nondim, modeshapes_dmd(k,:));
        else
          pjk = gauss_1pt_quad(bndry_pts_dmd, p_corr_nondim_exist, modeshapes_dmd(k,:));
        end
        P(l,k) = pjk;
      end
    end
    % Compute stability matrix
    Lmat = (L*(omega0^2)/2).*eye(num_modes);
    A = K + P - Lmat;
    A = sym(A);
    
    % Change relevant entry
    syms wpp
    A(j,j) = K(j,j) + P(j,j) - (L*wpp^2)/2;
    eqn = det(A) == 0;
    B = double(solve(eqn));
    
    omegap(j) = B(real(B)>0);
  end

  disp('Frequency solutions:');
  disp(omegap);
  omega_soln(m,:) = omegap;

end

f_results = figure; f_results.Position = [708 -807 1020 483];
for j = 1:num_modes
  s(j) = plot(real(omega_soln(:,j)), imag(omega_soln(:,j)), '--k', 'LineWidth', 2); grid on; hold on;
  s(j+num_modes) = scatter(real(omega_soln(:,j)), imag(omega_soln(:,j)), 72, M_list, 'filled');
end
scatter(real(omega_nat), imag(omega_nat), 96, 'filled');
cb = colorbar; colormap('jet'); caxis([M_list(1), M_list(end)]); 
title(cb,'$M$','interpreter', 'latex', 'fontsize',fs);
ax = gca; ax.FontSize = afs;
xlabel('Re( $\omega$ )', 'interpreter', 'latex', 'fontsize', fs);
ylabel('Im( $\omega$ )', 'interpreter', 'latex', 'fontsize', fs);
if (interp_on == 0)
  title(strcat(['Solution to Eigenfrequencies, $M = ', num2str(M_list(1)), '$ to $M = ', num2str(M_list(end)), '$, $N_m = ', num2str(num_modes), '$']), 'interpreter', 'latex', 'fontsize', fs);
else
  title(strcat(['Interpolated Solution, $M = ', num2str(M_list(1)), '$ to $M = ', num2str(M_list(end)), '$, $N_m = ', num2str(num_modes), '$']), 'interpreter', 'latex', 'fontsize', fs);
end


mode_num = 1; mach_ind = 1;
ret = plot_soln(col_pts, p_bem_orig, p_bem_soln, mode_num, mach_ind, M_list);

function [p] = p_lpt_func(mu, M, omega_in, x, bndry_pts, W, dWdx)

  N = length(x); p = zeros(N,1);
  for j = 1:N
    p(j) = mu*M*dWdx(j) - i*mu*omega_in*W(j);
  end

end

function [p] = p_upt_func(mu, M, omega, x, bndry_pts, W, dWdx)

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

function [pjk] = rect_1pt(x, p, psi)

  Np = length(p); dx = x(2) - x(1);
  pjk = 0;
  for j = 1:Np-1
    pjk = pjk + p(j)*psi(j)*dx;
  end

end

