
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
M = 1.5;
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

m_choice = num2str(M); m_choice = strcat([m_choice(1), m_choice(3:end)]);

% Define trial modes
num_modes = 2;
N = 241; x = linspace(0,L,N); x_dmd = linspace(0,L,61);
modeshapes = zeros(num_modes,N);

bndry_pts = x;
col_pts = (bndry_pts(1:end-1) + bndry_pts(2:end)) ./ 2;

modeshape = @(n,x_ind,L) sin(n*pi.*x_ind ./ L);
mode_deriv = @(n,x_ind,L) (n*pi/L).*cos(n*pi.*x_ind ./ L);
for j = 1:num_modes
  modeshapes(j,:) = modeshape(j,x,L);
  mode_derivs(j,:) = mode_deriv(j,x,L);
end

figure;
subplot(211);
plot(x, modeshapes(1,:), 'b', 'LineWidth', 2); grid on; hold on;
plot(x, mode_derivs(1,:), 'r', 'LineWidth', 2);
ax = gca; ax.FontSize = afs;
ylabel('$\psi_1 (x/h)$', 'interpreter', 'latex', 'fontsize', fs);
subplot(212);
plot(x, modeshapes(2,:), 'b', 'LineWidth', 2); grid on; hold on;
plot(x, mode_derivs(2,:), 'r', 'LineWidth', 2);
ax = gca; ax.FontSize = afs;
ylabel('$\psi_2 (x/h)$', 'interpreter', 'latex', 'fontsize', fs);
xlabel('$x / h$', 'interpreter', 'latex', 'fontsize', fs);
sgtitle('Trial Modes to Compute Stability of Simply-Supported Beam', 'interpreter', 'latex', 'fontsize', fs);

mu = 12e-5; D = 23.9;

%% Construct stability matrix
K = eye(num_modes);
for j = 1:num_modes
  K(j,j) = (D*(j*pi/L)^4 + (Mw^2)*(j*pi/L)^2)*(L/2);
end

% Compute natural frequencies
omega0 = zeros(1, num_modes);
for j = 1:num_modes
  omega0(j) = sqrt(D*(j*pi/L)^4 + (Mw^2)*(j*pi/L)^2);
end

% Compute pressure matrix using DMD
P = zeros(num_modes, num_modes);
for j = 1:num_modes
  omega = omega0(j);
  for l = 1:num_modes
    p_lpt_nondim = p_lpt_func(mu, M, omega, x, bndry_pts, modeshapes(l,:), mode_derivs(l,:));
    p_uns = p_upt_func(mu, M, omega, x, bndry_pts, modeshapes(l,:), mode_derivs(l,:));
    p_lpt = p_lpt_nondim .* (rho_m*a^2);
    dmd_mode = dlmread(strcat(['./2d_dmd_modes/m', m_choice, '/sin', num2str(l), '_omega', num2str(j), '/dmd_modes.csv']));
    alpha_tilde = dlmread(strcat(['./2d_dmd_modes/m', m_choice, '/sin', num2str(l), '_omega', num2str(j), '/alphas_tilde.csv']));
    alphas = dlmread(strcat(['./2d_dmd_modes/m', m_choice, '/sin', num2str(l), '_omega', num2str(j), '/alphas.csv']));
    dmd_mode = dmd_mode(:,2); alpha_tilde = alpha_tilde(2);
    p_corr = p_lpt + p.*(-i*2*h).*dmd_mode.*alpha_tilde;
    p_corr_nondim = p_corr ./ (rho_m*a^2);

    %  p_lpt_nondim = p_lpt_func(mu, M, omega, x, bndry_pts, modeshapes(l,:), mode_derivs(l,:));
    %  p_uns = p_upt_func(mu, M, omega, x, bndry_pts, modeshapes(l,:), mode_derivs(l,:));
    %  dmd_mode = dlmread(strcat(['./2d_dmd_modes/m', m_choice, '/sin', num2str(l), '_omega', num2str(j), '/dmd_modes.csv']));
    %  alpha_tilde = dlmread(strcat(['./2d_dmd_modes/m', m_choice, '/sin', num2str(l), '_omega', num2str(j), '/alphas_tilde.csv']));
    %  alphas = dlmread(strcat(['./2d_dmd_modes/m', m_choice, '/sin', num2str(l), '_omega', num2str(j), '/alphas.csv']));
    %  dmd_mode = dmd_mode(:,2); alpha_tilde = alpha_tilde(2);
    %  p_corr = p.*(-i*2*h).*dmd_mode.*alpha_tilde; p_corr_nondim = p_corr ./ (rho_m*a^2);

    f = figure; f.Position = [356 324 774 545];
    subplot(211); plot(x, real(p_lpt_nondim), 'LineWidth', 2); grid on; hold on; plot(x,real(p_corr_nondim), 'LineWidth', 2); plot(x, real(p_uns), 'LineWidth', 2);
    ax = gca; ax.FontSize = afs;
    ylabel('Re( $\hat{p}^{\prime} / \rho_m a^2$ )', 'interpreter', 'latex', 'fontsize', fs);
    legend('LPT', 'DMD', 'UPT', 'interpreter', 'latex', 'fontsize', fs);
    subplot(212); plot(x, imag(p_lpt_nondim), 'LineWidth', 2); grid on; hold on; plot(x,imag(p_corr_nondim), 'LineWidth', 2); plot(x, imag(p_uns), 'LineWidth', 2);
    ax = gca; ax.FontSize = afs;
    ylabel('Im( $\hat{p}^{\prime} / \rho_m a^2$ )', 'interpreter', 'latex', 'fontsize', fs);
    xlabel('$x / h$', 'interpreter', 'latex', 'fontsize', fs);
    sgtitle(strcat(['Comparison of Pressure Fluctuation from Theories, $\psi_', num2str(l), '$, $\omega_', num2str(j), '$, $M = ', num2str(M), '$']), 'interpreter', 'latex', 'fontsize', fs);

    for k = 1:num_modes
      ploc = rect_1pt(x, p_corr_nondim, modeshapes(k,:));
      P(l,k) = ploc;
    end
  end
  Lmat = (L*(omega^2)/2) .* eye(num_modes);
  A = K + P - Lmat;
  A = sym(A);
  syms wpp
  A(j,j) = K(j,j) + P(j,j) - (L*wpp^2)/2;
  eqn = det(A) == 0;
  B = double(solve(eqn));
  omegap(j) = B(real(B)>0);
end

disp('omegap: '); disp(omegap);

function [p] = p_lpt_func(mu, M, omega, x, bndry_pts, W, dWdx)

  N = length(x); p = zeros(N,1);
  for j = 1:N
    x_loc = x(j);
    p(j) = (mu*M/sqrt(M^2 - 1))*(-i*omega*W(j) + M*dWdx(j));
    p(j) = mu*M*dWdx(j) - i*mu*omega*W(j); % To match Zhang local piston theory
    %p(j) = mu*M*dWdx(j) + i*mu*omega*W(j);
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
