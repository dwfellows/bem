
close all;
clear all;
clc;

addpath('./utilities/2d/');

fs = 20;
afs = 18;
lfs = 20;

% Define problem parameters
case_choice = '22'; 
sin_choice = case_choice(1); sc = str2num(sin_choice); omega_choice = case_choice(2); oc = str2num(omega_choice);
m_choice = '14'; amp_choice = 'origamp';   % for c6: 05, 085, 15, 3; amps: twice/orig/half/fifth
M = str2num(strcat([m_choice(1), '.', m_choice(2:end)]));
amp = 0.00336757;
if (oc == 1)
  omega = 253.68;
else
  omega = 1014.7;
end

flow.pref = 70185.673999999;
flow.rhoref = 0.91;
flow.gamma = 1.4;
flow.Tref = 268.7355898456943;
flow.M = M;

% Define deformation function -- Vedeneev approximation of first natural beam mode
L_panel = 0.3;
flow.L = L_panel;

R = flow.pref ./ (flow.rhoref * flow.Tref);
k_wav = omega / sqrt(flow.gamma * R * flow.Tref);

psi = @(n,x_in,Lp) sin(n.*pi.*x_in./Lp);
dpsi_dx = @(n,x_in,Lp) (n*pi/Lp) .* cos(n.*pi.*x_in./Lp);

% Define discretization
N = 241; % Number of points on panel where velocity perturbation is computed

bndry_pts = linspace(0,L_panel,N+1); % locations of boundary element endpoints
dx = bndry_pts(2) - bndry_pts(1); 
col_pts = (bndry_pts(1:end-1) + bndry_pts(2:end)) ./ 2; % location of collocation points

psi_in = psi(sc,col_pts,L_panel);
dpsi_dx_in = dpsi_dx(sc,col_pts,L_panel);

figure;
subplot(211);
plot(col_pts, psi_in, 'b', 'LineWidth', 2); grid on;
ax = gca; ax.FontSize = afs;
ylabel('$\psi$', 'interpreter', 'latex', 'fontsize', fs);
subplot(212);
plot(col_pts, dpsi_dx_in, 'b', 'LineWidth', 2); grid on;
ax = gca; ax.FontSize = afs;
ylabel('$\psi^{\prime}$', 'interpreter', 'latex', 'fontsize', fs);
xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
sgtitle('C6 Continuous Bump', 'interpreter', 'latex', 'fontsize', fs);

[phi, p_bar] = bem2d(col_pts, bndry_pts, psi_in, dpsi_dx_in, amp, omega, flow);

f_soln = figure; f_soln.Position = [244 215 885 475];
subplot(211);
plot(col_pts ./ L_panel, real(phi), 'b', 'LineWidth', 2); hold on; grid on;
plot(col_pts ./ L_panel, imag(phi), '--k', 'LineWidth', 2);
ax = gca; ax.FontSize = afs;
legend('Real', 'Imag.', 'interpreter', 'latex', 'fontsize', lfs, 'location', 'northwest');
xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$\bar{\phi}^{\prime}$', 'interpreter', 'latex', 'fontsize', fs);
subplot(212);
plot(col_pts ./ L_panel, real(p_bar) ./ flow.pref, 'b', 'LineWidth', 2); hold on; grid on;
plot(col_pts ./ L_panel, imag(p_bar) ./ flow.pref, '--k', 'LineWidth', 2);
ax = gca; ax.FontSize = afs;
xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$\bar{p}^{\prime} / p_{in}$', 'interpreter', 'latex', 'fontsize', fs);
sgtitle(strcat(['Velocity Perturbation Potential and Spatial Pressure Mode, C6 Bump, $M = ', num2str(M), '$']), 'interpreter', 'latex', 'fontsize', fs);



% Set up video comparison against PC2 data
pc2_dt_nondim = 4.5970 * 10^(-5);
Uref = sqrt(flow.pref / flow.rhoref);
pref = flow.pref;
t_scale = 1 / Uref;
pc2_dt = pc2_dt_nondim * t_scale;
N_tsteps = 500000;
snap_level = 1;
snapshot = 100*snap_level;
dmd_dt = snapshot*pc2_dt;
snapshot_begin = 1;
snapshot_end = 5001;
t = snapshot_begin:snap_level:snapshot_end; 
t_phys = t.*dmd_dt;
Nt = length(t);

pf_bem_history = (p_bar.') * exp(i*omega*t_phys); pf_bem_history = (pf_bem_history + conj(pf_bem_history)) ./ flow.pref;
pf_bem_history = pf_bem_history';
yloc_history = (-i*(amp/2).*psi_in.') * exp(i*omega*t_phys); yloc_history = (yloc_history + conj(yloc_history));
yloc_history = yloc_history';
ydot_history = ((amp/2)*omega.*psi_in.') * exp(i*omega*t_phys); ydot_history = (ydot_history + conj(ydot_history)) ./ Uref;
ydot_history = ydot_history';

dlmwrite(strcat(['~/research/dmdc_beam/master_data_sinusoid/m', m_choice, '_bem/sin', sin_choice, '_omega', omega_choice, '/surface_pressure_history_sin', sin_choice, '_omega', omega_choice, '_m', m_choice, '.csv']), pf_bem_history+1, 'delimiter', ',', 'precision', 12);
dlmwrite(strcat(['~/research/dmdc_beam/master_data_sinusoid/m', m_choice, '_bem/sin', sin_choice, '_omega', omega_choice, '/surface_yloc_history_sin', sin_choice, '_omega', omega_choice, '_m', m_choice, '.csv']), yloc_history, 'delimiter', ',', 'precision', 12);
dlmwrite(strcat(['~/research/dmdc_beam/master_data_sinusoid/m', m_choice, '_bem/sin', sin_choice, '_omega', omega_choice, '/surface_ydot_history_sin', sin_choice, '_omega', omega_choice, '_m', m_choice, '.csv']), ydot_history, 'delimiter', ',', 'precision', 12);
dlmwrite(strcat(['~/research/dmdc_beam/master_data_sinusoid/m', m_choice, '_bem/sin', sin_choice, '_omega', omega_choice, '/surface_xloc_history_sin', sin_choice, '_omega', omega_choice, '_m', m_choice, '.csv']), yloc_history, 'delimiter', ',', 'precision', 12);
dlmwrite(strcat(['~/research/dmdc_beam/master_data_sinusoid/m', m_choice, '_bem/sin', sin_choice, '_omega', omega_choice, '/surface_xdot_history_sin', sin_choice, '_omega', omega_choice, '_m', m_choice, '.csv']), ydot_history, 'delimiter', ',', 'precision', 12)
disp('Data written. Pausing.'); pause;

max_pf_bem = max(max(pf_bem_history)); min_pf_bem = min(min(pf_bem_history));
f_vid = figure; f_vid.Position = [627 -851 908 602];
v_fname = strcat(['./video_results/2d/m', m_choice, '_sinusoid_Ne', num2str(N), '_x121_comp.csv']);
v_data = VideoWriter(v_fname, 'MPEG-4');
v_data.FrameRate = 60; v_data.Quality = 10;
open(v_data);
for m = 1:length(t)
  t_loc = t_phys(m);
  p_loc = pf_bem_history(:,m); 
  y_loc = yloc_history(:,m);

  % Plot pressure
  subplot(211);
  plot(col_pts ./ L_panel, p_loc, 'b', 'LineWidth', 2); hold on; grid on;
  ax = gca; ax.FontSize = afs;
  ylim([min_pf_bem, max_pf_bem]);
  ylabel('$p^{\prime} / p_{in}$', 'interpreter', 'latex', 'fontsize', fs);
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);

  subplot(212);
  plot(col_pts ./ L_panel, y_loc, 'b', 'LineWidth', 2); grid on;
  ax = gca; ax.FontSize = afs;
  ylim([-amp, amp]);
  ylabel('$w(x)$', 'interpreter', 'latex', 'fontsize', fs);
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);

  t_fname = strcat(['Comparison of CFD and Theories, C6 Bump Geometry, $M = ', m_choice(1), '.', m_choice(2:end), ', t = ', num2str(m), '/', num2str(length(t)), '$']);
  sgtitle(t_fname, 'interpreter', 'latex', 'fontsize', fs);

  frame = getframe(gcf);
  writeVideo(v_data, frame);

  clf;

end
close(v_data);

