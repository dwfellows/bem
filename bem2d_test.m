
close all;
clear all;
clc;

addpath('./utilities/2d/');

fs = 26;
afs = 18;
lfs = 20;
dmd_on = 0;

% Define problem parameters
m_choice = '15'; amp_choice = 'origamp';
M = str2num(strcat([m_choice(1), '.', m_choice(2:end)]));
if (M < 1)
  amp = 0.000722158;
elseif (M > 1)
  amp = 0.00336757;
end

if (strcmp(amp_choice, 'halfamp'))
  amp = amp / 2;
elseif (strcmp(amp_choice, 'twiceamp'))
  amp = amp*2;
elseif (strcmp(amp_choice, 'fifthamp'))
  amp = amp / 5;
elseif (strcmp(amp_choice, 'tenthamp'))
  amp = amp / 10;
end
omega = 376.99111843;

flow.pref = 97326.6883347;
flow.rhoref = 1.2792;
flow.gamma = 1.44;
flow.Tref = 273;
flow.M = M;
R = flow.pref / (flow.rhoref * flow.Tref);
k_wav = omega / sqrt(flow.rhoref * R * flow.Tref);

% Define deformation function -- Vedeneev approximation of first natural beam mode
c1 = 4.73; L_panel = 0.3; flow.L = L_panel;
psi = @(x,Lp) cos(c1.*x./Lp) - cosh(c1.*x./Lp) - ((cos(c1) - cosh(c1))/(sin(c1) - sinh(c1))).*(sin(c1.*x./Lp) - sinh(c1.*x./Lp));
psi_max = psi(L_panel/2, L_panel);
psi = @(x,Lp) (1/psi_max) .* ( cos(c1.*x./Lp) - cosh(c1.*x./Lp) - ((cos(c1) - cosh(c1))/(sin(c1) - sinh(c1))).*(sin(c1.*x./Lp) - sinh(c1.*x./Lp)) );
dpsi_dx = @(x,Lp) (1/psi_max) .* (-c1/Lp) .* ( ((cos(c1) - cosh(c1))/(sin(c1) - sinh(c1))).*(cos(c1.*x./Lp) - cosh(c1.*x./Lp)) + sin(c1.*x./Lp) + sinh(c1.*x./Lp) );
dpsi_dx2 = @(x,Lp) (1./psi_max).*((c1/Lp)^2).*( cos(c1.*x./Lp) + cosh(c1.*x./Lp) - ((cos(c1) - cosh(c1))/(sin(c1)-sinh(c1))).*(sin(c1.*x./Lp) + sinh(c1.*x./Lp)) );
n1 = @(q,x,Lp) -q.*dpsi_dx(x,Lp).*(q.^2 .* (dpsi_dx(x,Lp)).^2 + 1).^(-1/2);
n2 = @(q,x,Lp) ((q.^2) .* (dpsi_dx(x,Lp)).^2 + 1).^(-1/2);
dn1_dq = @(q,x,Lp) -1.*dpsi_dx(x,Lp).*((q.^2).*(dpsi_dx(x,Lp)).^2 + 1).^(-3/2);
dn2_dq = @(q,x,Lp) -q.*((dpsi_dx(x,Lp)).^2).*((q.^2).*((dpsi_dx(x,Lp)).^2) + 1).^(-3/2);

% Define discretization
N = 1001; % Number of points on panel where velocity perturbation is computed

bndry_pts = linspace(0,L_panel,N+1); % locations of boundary element endpoints
dx = bndry_pts(2) - bndry_pts(1); 
col_pts = (bndry_pts(1:end-1) + bndry_pts(2:end)) ./ 2; % location of collocation points

psi_in = psi(col_pts, L_panel);
dpsi_dx_in = dpsi_dx(col_pts, L_panel);

figure;
subplot(211); plot(col_pts, psi_in); grid on; ylabel('$\psi$', 'interpreter', 'latex', 'fontsize', fs);
subplot(212); plot(col_pts, dpsi_dx_in); grid on; ylabel('$\partial \psi / \partial x$', 'interpreter', 'latex', 'fontsize', fs);

figure; Uinf = flow.M*sqrt(flow.gamma*287*flow.Tref);
subplot(211); plot(col_pts, amp.*Uinf.*dpsi_dx_in); grid on; hold on; ylabel('$A \cdot U_{\infty} \cdot \frac{\partial \psi}{\partial x}$', 'interpreter', 'latex', 'fontsize', fs);
subplot(212); plot(col_pts, amp.*omega.*psi_in); grid on; ylabel('$A \cdot \omega \cdot \psi$', 'interpreter', 'latex', 'fontsize', fs);
xlabel('$x$', 'interpreter', 'latex', 'fontsize', fs);
sgtitle(strcat(['$M = ', num2str(M), '$']), 'interpreter', 'latex', 'fontsize', fs);

M_comparison = 0;
if (M_comparison == 1)
  M_list = [0.5, 0.65, 0.8, 0.85]; Nm = max(size(M_list));
  %M_list = [1.5, 2.0, 2.5, 3.0]; Nm = max(size(M_list));
  phi_list = zeros(Nm, N);
  pbar_list = zeros(Nm, N);
  for j = 1:Nm
  
    M = M_list(j);
    flow.M = M;
    [phi, p_bar] = bem2d(col_pts, bndry_pts, psi_in, dpsi_dx_in, amp, omega, flow);
    phi_list(j,:) = phi;
    pbar_list(j,:) = p_bar;
  
  end
  
  f_comp = figure; f_comp.Position = [244 215 885 475];
  subplot(221);
  colors = ['r', 'b', 'g', 'm'];
  for j = 1:Nm
    color = colors(j);
    l(2*(j-1)+1) = plot(col_pts ./ L_panel, real(phi_list(j,:)), color,  'LineWidth', 2, 'DisplayName', strcat(['$M = ', num2str(M_list(j)), '$'])); grid on; hold on;
  end
  ax = gca; ax.FontSize = afs;
  ylabel('Re( $\bar{\phi}^{\prime}$ )', 'interpreter', 'latex', 'fontsize', fs);
  subplot(222);
  for j = 1:Nm
    color = colors(j);
    l(2*(j-1)+2) = plot(col_pts ./ L_panel, imag(phi_list(j,:)), color, 'LineWidth', 2); grid on; hold on;
  end
  ax = gca; ax.FontSize = afs;
  legend([l(1), l(3), l(5), l(7)], 'interpreter', 'latex', 'fontsize', fs, 'location', 'northwest');
  ylabel('Im( $\bar{\phi}^{\prime}$ )', 'interpreter', 'latex', 'fontsize', fs);

  subplot(223);
  for j = 1:Nm
    color = colors(j);
    plot(col_pts ./ L_panel, real(pbar_list(j,:)) ./ flow.pref, color, 'LineWidth', 2); grid on; hold on;
  end
  ax = gca; ax.FontSize = afs;
  ylabel('Re( $\bar{p}^{\prime}$ ) / $p_{\infty}$', 'interpreter', 'latex', 'fontsize', fs);
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
  subplot(224);
  for j = 1:Nm
    color = colors(j);
    plot(col_pts ./ L_panel, imag(pbar_list(j,:)) ./ flow.pref, color, 'LineWidth', 2); grid on; hold on;
  end
  ax = gca; ax.FontSize = afs;
  ylabel('Im( $\bar{p}^{\prime}$ ) / $p_{\infty}$', 'interpreter', 'latex', 'fontsize', fs);
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);

  if (M_list(1) < 1)
    sgtitle('Reduced Velocity Perturbation Potential and Spatial Pressure Mode, 2D Configuration, Subsonic Flows', 'interpreter', 'latex', 'fontsize', fs);
  else
    sgtitle('Reduced Velocity Perturbation Potential and Spatial Pressure Mode, 2D Configuration, Supersonic Flows', 'interpreter', 'latex', 'fontsize', fs);
  end

  pause;
end

[phi, p_bar] = bem2d(col_pts, bndry_pts, psi_in, dpsi_dx_in, amp, omega, flow);

% Compute LPT modes
p_bar_lpt = (-i*(amp/2)*flow.gamma*flow.M*dpsi_dx_in + flow.gamma*k_wav*(amp/2)*psi_in);  % p'_{LPT} / p_{in}
figure;
plot(col_pts ./ L_panel, real(p_bar_lpt), 'b', 'LineWidth', 2); grid on; hold on;
plot(col_pts ./ L_panel, imag(p_bar_lpt), '--k', 'LineWidth', 2);

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
sgtitle(strcat(['Velocity Perturbation Potential and Spatial Pressure Mode, $M = ', num2str(M), '$']), 'interpreter', 'latex', 'fontsize', fs);

% Read in DMD mode
if (dmd_on == 1)
  p_bar_dmd = dlmread(strcat(['./mode1_modes_clamped/m', m_choice, '_modes_', amp_choice, '.csv']));
  amplitude = dlmread(strcat(['./mode1_modes_clamped/m', m_choice, '_amps_', amp_choice, '.csv']));
  p_bar_dmd = [p_bar_dmd(:,1)*amplitude(1), p_bar_dmd(:,2)*amplitude(2)];

  dmd_choice = 1;

  bem_mode = p_bar./flow.pref - p_bar_lpt;
  max_re_bem = max(real(bem_mode)); min_re_bem = min(real(bem_mode));
  max_re_dmd = max(real(p_bar_dmd(:,dmd_choice))); min_re_dmd = min(real(p_bar_dmd(:,dmd_choice)));
  max_im_bem = max(imag(bem_mode)); min_im_bem = min(imag(bem_mode));
  max_im_dmd = max(imag(p_bar_dmd(:,dmd_choice))); min_im_dmd = min(imag(p_bar_dmd(:,dmd_choice)));
  max_re = max([max_re_bem, max_re_dmd]); min_re = min([min_re_bem, min_re_dmd]); 
  max_im = max([max_im_bem, max_im_dmd]); min_im = min([min_im_bem, min_im_dmd]);
  max_lim = max([max_re, max_im]); min_lim = min([min_re, min_im]);

  f_full_comp = figure; f_full_comp.Position = [244 215 885 475];
  pdx = linspace(0,L_panel,61);
  subplot(121);
  plot(col_pts ./ L_panel, real(p_bar)./flow.pref - real(p_bar_lpt), 'b', 'LineWidth', 2); grid on; hold on;
  plot(pdx ./ L_panel, real(p_bar_dmd(:,1)), 'r', 'LineWidth', 2);
  ax = gca; ax.FontSize = afs; ylim([min_lim, max_lim]);
  legend('BEM', 'DMD', 'interpreter', 'latex', 'fontsize', fs, 'location', 'northwest');
  ylabel('Re ( $\bar{p}^{\prime} / p_{\infty}$ )', 'interpreter', 'latex', 'fontsize', fs);
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
  
  subplot(122);
  plot(col_pts ./ L_panel, imag(p_bar)./flow.pref - imag(p_bar_lpt), 'b', 'LineWidth', 2); grid on; hold on;
  plot(pdx ./ L_panel, imag(p_bar_dmd(:,1)), 'r', 'LineWidth', 2);
  ax = gca; ax.FontSize = afs; ylim([min_lim, max_lim]);
  ylabel('Im ( $\bar{p}^{\prime} / p_{\infty}$ )', 'interpreter', 'latex', 'fontsize', fs);
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
  sgtitle(strcat(['Comparison of DMD Mode to BEM, $M = ', num2str(flow.M), '$']), 'interpreter', 'latex', 'fontsize', fs);
  
  pause;
end

pc2_dt_nondim = 4.5970 * 10^(-5);
Uref = sqrt(flow.pref / flow.rhoref); 
R = flow.pref / (flow.rhoref * flow.Tref);
aref = sqrt(flow.gamma * R * flow.Tref);
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


%% Compare against PC2 simulation data
pressure_data = dlmread(strcat(['./comparison_data/2d/m', m_choice, '/', amp_choice, '/surface_pressure_history_m', m_choice, '_', amp_choice, '.csv'])); pressure_data = pressure_data';
y_data = dlmread(strcat(['./comparison_data/2d/m', m_choice, '/', amp_choice, '/surface_yloc_history_m', m_choice, '_', amp_choice, '.csv'])); y_data = y_data';

% pressure_data = dlmread('./comparison_data/2d/m085/x241/surface_pressure_history_m085_x241.csv'); pressure_data = pressure_data';
% y_data = dlmread('./comparison_data/2d/m085/x241/surface_yloc_history_m085_x241.csv'); y_data = y_data';

pc2_dt_nondim = 4.5970 * 10^(-5);
Uref = sqrt(flow.pref / flow.rhoref); 
R = flow.pref / (flow.rhoref * flow.Tref);
aref = sqrt(flow.gamma * R * flow.Tref);
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

p_size = size(pressure_data); pdNx = p_size(1); pdNt = p_size(2);
pdx = linspace(0,L_panel,pdNx);
pdx_dmd = linspace(0,L_panel,pdNx);

pf_bem = (p_bar').*exp(i*omega.*t_phys); pf_bem = (pf_bem + conj(pf_bem))./flow.pref; 
pf_data = pressure_data - 1;
pf_lpt = (p_bar_lpt').*exp(i*omega.*t_phys); pf_lpt = pf_lpt + conj(pf_lpt);
max_pf_bem = max(max(pf_bem)); min_pf_bem = min(min(pf_bem));
max_pf_cfd = max(max(pf_data(:,2:end))); min_pf_cfd = min(min(pf_data(:,2:end)));
max_pf_lpt = max(max(pf_lpt)); min_pf_lpt = min(min(pf_lpt));
max_pf = max([max_pf_bem, max_pf_cfd, max_pf_lpt]);
min_pf = min([min_pf_bem, min_pf_cfd, min_pf_lpt]);

diff = ((max_pf_bem - max_pf_cfd)/max_pf_cfd)*100 ; disp(diff);

disp('Preparing to begin video comparison to PC2. Press any button to begin.'); pause;

pf_bem_history = imag( (p_bar)' * exp(i*omega*t_phys) ) ./ flow.pref;
f_cfd_comp = figure; f_cfd_comp.Position = [521 -787 1093 583];
plot(col_pts ./ L_panel, pf_bem_history(:, 2251), 'b', 'LineWidth', 2); grid on; hold on;
plot(pdx ./ L_panel, pressure_data(:,2251) - 1, '--k', 'LineWidth', 2);
ax = gca; ax.FontSize = afs;
legend('BEM', 'CFD', 'interpreter', 'latex', 'fontsize', fs, 'location', 'southeast');
xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$p^{\prime} / p_{\infty}$', 'interpreter', 'latex', 'fontsize', fs);
title(strcat(['Pressure Comparison, 2D Subsonic BEM and CFD, $M = ', num2str(M), '$']), 'interpreter', 'latex', 'fontsize', fs);
pause;

% Construct video of subsonic pressure behavior and discrepancy from local piston theory
pf_err_history = zeros(N, Nt);

max_pf = 1.05*max_pf;
min_pf = 1.05*min_pf;

f_vid = figure; f_vid.Position = [523 -891 919 590];
if (dmd_on == 0) v_fname = strcat(['./video_results/2d/m', m_choice, '_', amp_choice, '_Ne', num2str(N), '_x241.mp4']);
else v_fname = strcat(['./video_results/2d/m', m_choice, '_', amp_choice, '_Ne', num2str(N), '_DMD.mp4']); end
v_data = VideoWriter(v_fname, 'MPEG-4');
v_data.FrameRate = 20; v_data.Quality = 12;
open(v_data);
for m = 1:length(t)
  t_loc = t_phys(m);
  p_loc = p_bar.*exp(i*omega*t_loc); 
  p_loc = p_loc + conj(p_loc);  % complex-conjugate method
  if (dmd_on == 1)
    p_loc_dmd = p_bar_dmd(:,1).*exp(i*omega*t_loc);
    p_loc_dmd = p_loc_dmd + conj(p_loc_dmd);
  end

  % Compute LPT prediction
  q_loc = amp*sin(omega*t_loc); qdot_loc = amp*omega*cos(omega*t_loc); qddot_loc = -(omega^2)*q_loc;
  w_mis = (M*sqrt(flow.gamma*R*flow.Tref)).*-n1(q_loc,pdx,L_panel);
  w_mot = qdot_loc.*psi(pdx,L_panel).*n2(q_loc,pdx,L_panel);
  w = w_mis + w_mot;
  p = flow.pref.*ones(1,pdNx) + (flow.rhoref*aref).*w;
  p_lpt = p ./ flow.pref; p_lpt = p_lpt - 1;

  % Plot pressure
  subplot(211);
  if (dmd_on == 0)
    plot(col_pts ./ L_panel, p_loc ./ flow.pref, 'b', 'LineWidth', 2); hold on; grid on;
    plot(col_pts ./ L_panel, p_lpt, 'm', 'LineWidth', 2);
    plot(pdx ./ L_panel, pressure_data(:,m) - 1, '--k', 'LineWidth', 2);
    legend('BEM', 'LPT', 'CFD', 'interpreter', 'latex', 'fontsize', fs, 'location', 'eastoutside');
    ylim([min_pf, max_pf]); ax = gca; ax.FontSize = afs;
    ylabel('$p^{\prime} / p_{\infty}$', 'interpreter', 'latex', 'fontsize', fs);
  else
    plot(col_pts ./ L_panel, p_loc ./ flow.pref, 'b', 'LineWidth', 2); hold on; grid on;
    plot(pdx_dmd ./ L_panel, p_loc_dmd + p_lpt', 'g', 'LineWidth', 2);
    %plot(col_pts ./ L_panel, p_lpt, 'm', 'LineWidth', 2);
    plot(pdx ./ L_panel, pressure_data(:,m) - 1, '--k', 'LineWidth', 2);
    legend('BEM', 'DMD', 'CFD', 'interpreter', 'latex', 'fontsize', fs, 'location', 'eastoutside');
    ylim([min_pf, max_pf]); ax = gca; ax.FontSize = afs;
    ylabel('$p^{\prime} / p_{\infty}$', 'interpreter', 'latex', 'fontsize', fs);
  end

  subplot(212);
  plot(pdx ./ L_panel, y_data(:,m), 'b', 'LineWidth', 2); grid on;
  ylim([-amp, amp]); ax = gca; ax.FontSize = afs;
  ylabel('$w(x)$', 'interpreter', 'latex', 'fontsize', fs);

  t_fname = strcat(['Comparison of CFD and Theories, $N_e = ', num2str(N), ', M = ', m_choice(1), '.', m_choice(2:end), ', t = ', num2str(m), '/', num2str(length(t)), '$']);
  sgtitle(t_fname, 'interpreter', 'latex', 'fontsize', fs);

  frame = getframe(gcf);
  writeVideo(v_data, frame);

  p_history(:,m) = p_loc';

  clf;

end
close(v_data);


