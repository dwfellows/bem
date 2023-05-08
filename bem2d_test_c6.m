
close all;
clear all;
clc;

addpath('./utilities/2d/');

fs = 20;
afs = 18;
lfs = 20;

dmd_on = 1;
Ndmd = 61; read_dmd = (Ndmd-1)*2 + 1;  % 61 or 121

% Define problem parameters
m_choice = '3'; amp_choice = 'origamp';   % for c6: 05, 085, 15, 3; amps: twice/orig/half/fifth
M = str2num(strcat([m_choice(1), '.', m_choice(2:end)]));
amp = 0.000722158;
omega = 376.991;

flow.pref = 97326.6883347;
flow.rhoref = 1.2792;
flow.gamma = 1.44;
flow.Tref = 273;
flow.M = M;

% Define deformation function -- Vedeneev approximation of first natural beam mode
L_panel = 0.3;
flow.L = L_panel;

R = flow.pref ./ (flow.rhoref * flow.Tref);
k_wav = omega / sqrt(flow.gamma * R * flow.Tref);

% Construct matrix
A = zeros(8,8);
A(1,:) = ones(1,8);
A(2,:) = 7:14;

for ii = 7:14
  A(3,ii-6) = ii*(ii-1);
  A(4,ii-6) = ii*(ii-1)*(ii-2);
  A(5,ii-6) = ii*(ii-1)*(ii-2)*(ii-3);
  A(6,ii-6) = ii*(ii-1)*(ii-2)*(ii-3)*(ii-4);
  A(7,ii-6) = ii*(ii-1)*(ii-2)*(ii-3)*(ii-4)*(ii-5);
  A(8,ii-6) = 0.5^ii;
end

b = zeros(8,1); b(8) = 1;

p = A\b;
p = [flip(p); zeros(7,1)];

pd = 14:-1:1; pd = pd.*p(1:end-1)';

% Define discretization
N = 801; % Number of points on panel where velocity perturbation is computed

bndry_pts = linspace(0,L_panel,N+1); % locations of boundary element endpoints
dx = bndry_pts(2) - bndry_pts(1); 
col_pts = (bndry_pts(1:end-1) + bndry_pts(2:end)) ./ 2; % location of collocation points

psi_in = polyval(p, col_pts ./ L_panel);
dpsi_dx_in = polyval(pd, col_pts ./ L_panel) ./ L_panel;

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

M_comparison = 0;
if (M_comparison == 1)
  %M_list = [0.5, 0.65, 0.8, 0.85]; Nm = max(size(M_list));
  M_list = [1.5, 2.0, 2.5, 3.0]; Nm = max(size(M_list));
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


% Compute spatial mode associated with LPT
p_bar_lpt = (-i*(amp/2)*flow.gamma*flow.M*dpsi_dx_in + flow.gamma*k_wav*(amp/2)*psi_in);  % p'_{LPT} / p_{in}
bem_mode = p_bar ./ flow.pref - p_bar_lpt;

% Read in DMD mode
if (dmd_on == 1)
  p_bar_dmd = dlmread(strcat(['./mode1_modes/x', num2str(read_dmd), '/m', m_choice, '_modes_', amp_choice, '_x', num2str(read_dmd), '_c6.csv']));
  amplitude = dlmread(strcat(['./mode1_modes/x', num2str(read_dmd), '/m', m_choice, '_amps_', amp_choice, '_x', num2str(read_dmd), '_c6.csv']));

  p_bar_dmd = [p_bar_dmd(:,1)*amplitude(1), p_bar_dmd(:,2)*amplitude(2)];

  p_bar_dmd_size = size(p_bar_dmd); pdmdNx = p_bar_dmd_size(1); pdmd_x = linspace(0,1,pdmdNx);
  psi_in_dmd = polyval(p, pdmdNx); dpsi_dx_in_dmd = polyval(pd, pdmdNx) ./ L_panel;
  p_bar_lpt_dmd = (-i*(amp/2)*flow.gamma*flow.M*dpsi_dx_in_dmd + flow.gamma*k_wav*(amp/2)*psi_in_dmd);

  dmd_choice = 1;
  max_re_bem = max(real(bem_mode)); min_re_bem = min(real(bem_mode));
  max_re_dmd = max(real(p_bar_dmd(:,dmd_choice))); min_re_dmd = min(real(p_bar_dmd(:,dmd_choice)));
  max_im_bem = max(imag(bem_mode)); min_im_bem = min(imag(bem_mode));
  max_im_dmd = max(imag(p_bar_dmd(:,dmd_choice))); min_im_dmd = min(imag(p_bar_dmd(:,dmd_choice)));
  max_re = max([max_re_bem, max_re_dmd]); min_re = min([min_re_bem, min_re_dmd]); 
  max_im = max([max_im_bem, max_im_dmd]); min_im = min([min_im_bem, min_im_dmd]);
  max_lim = max([max_re, max_im]); min_lim = min([min_re, min_im]);

  pdx = linspace(0,L_panel,Ndmd);

  f_full_comp = figure; f_full_comp.Position = [244 215 885 475];
  subplot(121);
  plot(col_pts ./ L_panel, real(bem_mode), 'b', 'LineWidth', 2); grid on; hold on;
  plot(pdx ./ L_panel, real(p_bar_dmd(:,dmd_choice)), 'r', 'LineWidth', 2);
  ax = gca; ax.FontSize = afs; ylim([min_lim, max_lim]);
  legend('BEM', 'DMD', 'interpreter', 'latex', 'fontsize', fs, 'location', 'northwest');
  ylabel('Re ( $\bar{p}^{\prime} / p_{\infty}$ )', 'interpreter', 'latex', 'fontsize', fs);
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
  
  subplot(122);
  plot(col_pts ./ L_panel, imag(bem_mode), 'b', 'LineWidth', 2); grid on; hold on;
  plot(pdx ./ L_panel, imag(p_bar_dmd(:,dmd_choice)), 'r', 'LineWidth', 2);
  ax = gca; ax.FontSize = afs; ylim([min_lim, max_lim]);
  ylabel('Im ( $\bar{p}^{\prime} / p_{\infty}$ )', 'interpreter', 'latex', 'fontsize', fs);
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
  sgtitle(strcat(['Comparison of DMD Mode to BEM, C6-Continuous Geometry, $M = ', num2str(flow.M), '$']), 'interpreter', 'latex', 'fontsize', fs);
end

disp('Pausing before video comparison. Press anything to begin process.'); pause;

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

if (read_dmd == 121)
  pressure_data = dlmread(strcat(['./comparison_data/2d/c6/x', num2str(read_dmd), '/m', m_choice, '/surface_pressure_history_m', m_choice, '_c6.csv'])); pressure_data = pressure_data';
  y_data = dlmread(strcat(['./comparison_data/2d/c6/x', num2str(read_dmd), '/m', m_choice, '/surface_yloc_history_m', m_choice, '_c6.csv'])); y_data = y_data';
elseif (read_dmd == 241)
  pressure_data = dlmread(strcat(['./comparison_data/2d/c6/x', num2str(read_dmd), '/m', m_choice, '/surface_pressure_history_m', m_choice, '_x241_c6.csv'])); pressure_data = pressure_data';
  y_data = dlmread(strcat(['./comparison_data/2d/c6/x', num2str(read_dmd), '/m', m_choice, '/surface_yloc_history_m', m_choice, '_x241_c6.csv'])); y_data = y_data';
end
p_size = size(pressure_data); pdNx = p_size(1); pdNt = p_size(2);
pdx = linspace(0,L_panel,pdNx);

pf_bem_history = (p_bar.') * exp(i*omega*t_phys); pf_bem_history = (pf_bem_history + conj(pf_bem_history)) ./ flow.pref;
f_cfd_comp = figure; f_cfd_comp.Position = [521 -787 1093 583];
plot(col_pts ./ L_panel, pf_bem_history(:, 2251), 'b', 'LineWidth', 2); grid on; hold on;
plot(pdx ./ L_panel, pressure_data(:,2251) - 1, '--k', 'LineWidth', 2);
ax = gca; ax.FontSize = afs;
legend('BEM', 'CFD', 'interpreter', 'latex', 'fontsize', fs, 'location', 'southeast');
xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$p^{\prime} / p_{\infty}$', 'interpreter', 'latex', 'fontsize', fs);
if (M < 1)
  title(strcat(['Pressure Comparison, 2D Subsonic BEM and CFD, $M = ', num2str(M), '$']), 'interpreter', 'latex', 'fontsize', fs);
else
  title(strcat(['Pressure Comparison, 2D Supersonic BEM and CFD, $M = ', num2str(M), '$']), 'interpreter', 'latex', 'fontsize', fs);
end

pf_bem = (p_bar').*exp(i*omega.*t_phys); pf_bem = (pf_bem + conj(pf_bem))./flow.pref;
pf_data = pressure_data - 1;
pf_lpt = (p_bar_lpt').*exp(i*omega.*t_phys); pf_lpt = pf_lpt + conj(pf_lpt);
max_pf_bem = max(max(pf_bem)); min_pf_bem = min(min(pf_bem));
max_pf_cfd = max(max(pf_data(:,2:end))); min_pf_cfd = min(min(pf_data(:,2:end)));
max_pf_lpt = max(max(pf_lpt)); min_pf_lpt = min(min(pf_lpt));
max_pf = max([max_pf_bem, max_pf_cfd, max_pf_lpt]);
min_pf = min([min_pf_bem, min_pf_cfd, min_pf_lpt]);

psi_in_lpt = polyval(p, pdx ./ L_panel); dpsi_dx_in_lpt = polyval(pd, pdx ./ L_panel) ./ L_panel;
p_bar_lpt = (-i*(amp/2)*flow.gamma*flow.M*dpsi_dx_in_lpt + flow.gamma*k_wav*(amp/2)*psi_in_lpt);  % p'_{LPT} / p_{in}

f_vid = figure; f_vid.Position = [627 -851 908 602];
if (read_dmd == 121)
  v_fname = strcat(['./video_results/2d/m', m_choice, '_c6geom_Ne', num2str(N), '_x121_comp.csv']);
elseif (read_dmd == 241)
  v_fname = strcat(['./video_results/2d/m', m_choice, '_c6geom_Ne', num2str(N), '_x241_comp.csv']);
end
v_data = VideoWriter(v_fname, 'MPEG-4');
v_data.FrameRate = 60; v_data.Quality = 10;
open(v_data);
for m = 1:length(t)
  t_loc = t_phys(m);
  p_loc = p_bar.*exp(i*omega*t_loc); p_loc = p_loc + conj(p_loc);
  p_lpt = p_bar_lpt.*exp(i*omega*t_loc); p_lpt = p_lpt + conj(p_lpt);

  if (dmd_on == 1)
    p_loc_dmd = p_bar_dmd(:,dmd_choice).*exp(i*omega*t_loc); p_loc_dmd = p_loc_dmd + conj(p_loc_dmd);
  end

  % Plot pressure
  subplot(211);
  plot(col_pts ./ L_panel, p_loc ./ pref, 'b', 'LineWidth', 2); hold on; grid on;
  if (dmd_on == 1) plot(pdx ./ L_panel, p_loc_dmd + p_lpt', 'g', 'LineWidth', 2); end
  plot(pdx ./ L_panel, pressure_data(:,m) - 1, '--k', 'LineWidth', 2);
  ax = gca; ax.FontSize = afs;
  if (dmd_on == 1)
    legend('BEM', 'DMD', 'CFD', 'interpreter', 'latex', 'fontsize', fs, 'location', 'eastoutside');
  else
    legend('BEM', 'CFD', 'interpreter', 'latex', 'fontsize', fs, 'location', 'eastoutside');
  end
  % ylim([-0.05, 0.05]);
  ylim([-0.002, 0.002]);
  ylabel('$p^{\prime} / p_{in}$', 'interpreter', 'latex', 'fontsize', fs);
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);

  subplot(212);
  plot(pdx ./ L_panel, y_data(:,m), 'b', 'LineWidth', 2); grid on;
  ax = gca; ax.FontSize = afs;
  ylim([-amp, amp]);
  ylabel('$w(x)$', 'interpreter', 'latex', 'fontsize', fs);
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);

  t_fname = strcat(['Comparison of CFD and Theories, C6 Bump Geometry, $M = ', m_choice(1), '.', m_choice(2:end), ', t = ', num2str(m), '/', num2str(length(t)), '$']);
  sgtitle(t_fname, 'interpreter', 'latex', 'fontsize', fs);

  frame = getframe(gcf);
  writeVideo(v_data, frame);

  if (m == 1251) pause; end

  p_history(:,m) = p_loc';

  clf;

end
close(v_data);

