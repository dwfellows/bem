
close all;
clear all;
clc;

addpath('./utilities/3d/');

fs = 20;
afs = 18;
lfs = 20;
sz = 36;

dmd_on = 1;

% Define problem parameters
m_choice = '14'; amp_choice = 'origamp';
M = str2num(strcat([m_choice(1), '.', m_choice(2:end)]));
if (strcmp(amp_choice, 'origamp'))
  amp = 0.0000508;
end
omega = 691.15038379;

% Define flow properties based on choice of Mach number
if (M == 1.1)
  flow.pref = 34601.89642944;
  flow.rhoref = 0.49774215;
  flow.Tref = 242.222;
elseif (M == 1.2)
  flow.pref = 28872.53992512;
  flow.rhoref = 0.43729082;
  flow.Tref = 230.05556;
elseif (M == 1.3)
  flow.pref = 35836.7243159;
  flow.rhoref = 0.55034266;
  flow.Tref = 226.88889;
else % M \geq 1.4
  flow.pref = 41596.29322233;
  flow.rhoref = 0.64847794;
  flow.Tref = 223.5;
end
flow.gamma = 1.4;
flow.M = M;
R = flow.pref / (flow.rhoref * flow.Tref);
k_wav = omega / sqrt(flow.gamma*R*flow.Tref);

%% Define BEM collocation and boundary points
La = 0.2286; Lb = 2*La; flow.La = La;
Nx = 61; Ny = 61; N = (Nx-1)*(Ny-1);
col_pts = zeros(N, 2);  bndry_pts = zeros(4*N, 2);
bnd_x = linspace(0,La,Nx); bnd_y = linspace(0,Lb,Ny);
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

figure;
scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, 'r'); hold on; grid on;
scatter(bndry_pts(:,1) ./ La, bndry_pts(:,2) ./ La, 'b');
pbaspect([1 2 1]);
xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
title('Collocation and Element Boundary Points', 'interpreter', 'latex', 'fontsize', fs);


%% Define deformation function -- FEM solution of panel first natural mode
% Define 2D data
c1 = 4.73; L_panel = La;
psi2d = @(x,Lp) cos(c1.*x./Lp) - cosh(c1.*x./Lp) - ((cos(c1) - cosh(c1))/(sin(c1) - sinh(c1))).*(sin(c1.*x./Lp) - sinh(c1.*x./Lp));
psi_max_2d = psi2d(L_panel/2, L_panel);
psi2d = @(x,Lp) (1/psi_max_2d) .* ( cos(c1.*x./Lp) - cosh(c1.*x./Lp) - ((cos(c1) - cosh(c1))/(sin(c1) - sinh(c1))).*(sin(c1.*x./Lp) - sinh(c1.*x./Lp)) );
dpsi_dx2d = @(x,Lp) (1/psi_max_2d) .* (-c1/Lp) .* ( ((cos(c1) - cosh(c1))/(sin(c1) - sinh(c1))).*(cos(c1.*x./Lp) - cosh(c1.*x./Lp)) + sin(c1.*x./Lp) + sinh(c1.*x./Lp) );

psi = zeros(N,1); dpsi_dx = zeros(N,1);
for j = 1:(Ny-1)

  strip_indices = (j-1)*(Nx-1)+1:(j-1)*(Nx-1)+(Nx-1);
  loc_amp = psi2d(col_y(j), Lb);
  strip_psi_anal = loc_amp.*psi2d(col_x, La);
  strip_dpsi_dx_anal = loc_amp.*dpsi_dx2d(col_x, La);

  psi(strip_indices) = strip_psi_anal;
  dpsi_dx(strip_indices) = strip_dpsi_dx_anal;

end

% %% test -- use C6 data
% % Construct matrix
% A = zeros(8,8);
% A(1,:) = ones(1,8);
% A(2,:) = 7:14;
% 
% for ii = 7:14
%   A(3,ii-6) = ii*(ii-1);
%   A(4,ii-6) = ii*(ii-1)*(ii-2);
%   A(5,ii-6) = ii*(ii-1)*(ii-2)*(ii-3);
%   A(6,ii-6) = ii*(ii-1)*(ii-2)*(ii-3)*(ii-4);
%   A(7,ii-6) = ii*(ii-1)*(ii-2)*(ii-3)*(ii-4)*(ii-5);
%   A(8,ii-6) = 0.5^ii;
% end
% 
% b = zeros(8,1); b(8) = 1;
% 
% p = A\b;
% p = [flip(p); zeros(7,1)];
% 
% pd = 14:-1:1; pd = pd.*p(1:end-1)';
% for j = 1:Ny-1
%   strip_indices = (j-1)*(Nx-1)+1:(j-1)*(Nx-1)+(Nx-1);
%   loc_amp = polyval(p, col_pts(strip_indices(1), 2)./Lb);
%   strip_psi_anal = loc_amp.*polyval(p, col_pts(strip_indices, 1)./La);
%   strip_dpsi_dx_anal = loc_amp.*polyval(pd, col_pts(strip_indices,1)./La);
%   psi(strip_indices) = strip_psi_anal;
%   dpsi_dx(strip_indices) = strip_dpsi_dx_anal;
% end

figure;
scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, psi); grid on;
colorbar; colormap('jet'); pbaspect([1 2 1]);
ax = gca; ax.FontSize = afs;
xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
title('$\psi (x,y)$', 'interpreter', 'latex', 'fontsize', fs);

figure;
scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, dpsi_dx); grid on;
colorbar; colormap('jet'); pbaspect([1 2 1]);
ax = gca; ax.FontSize = afs;
xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
title('$\partial \psi(x,y) / \partial x$', 'interpreter', 'latex', 'fontsize', fs);



%% Compute spatial solution of reduced velocity perturbation potential and reduced pressure fluctuation
[phibar_prime, pbar_prime] = bem3d(col_pts, bndry_pts, Nx-1, psi, dpsi_dx, amp, omega, flow);

% Compute LPT mode
p_bar_lpt = (-i*(amp/2)*flow.gamma*flow.M.*dpsi_dx' + flow.gamma*k_wav*(amp/2).*psi');

figure;
subplot(121);
scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, real(p_bar_lpt)); grid on;
ax = gca; ax.FontSize = afs; colorbar; colormap('jet'); pbaspect([1 2 1]);
ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
title('Re( $\bar{p}^{\prime} / p_{\infty}$ )', 'interpreter', 'latex', 'fontsize', fs);
subplot(122);
scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, imag(p_bar_lpt)); grid on;
ax = gca; ax.FontSize = afs; colorbar; colormap('jet'); pbaspect([1 2 1]);
ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
title('Im( $\bar{p}^{\prime} / p_{\infty}$ )', 'interpreter', 'latex', 'fontsize', fs);
sgtitle('LPT Mode', 'interpreter', 'latex', 'fontsize', fs);

if (dmd_on == 1)
  p_bar_dmd = dlmread(strcat(['./mode1_modes_3d/m', m_choice, '_modes_', amp_choice, '.csv']));
  amplitude = dlmread(strcat(['./mode1_modes_3d/m', m_choice, '_amps_', amp_choice, '.csv']));
  p_bar_dmd = [p_bar_dmd(:,1)*amplitude(1), p_bar_dmd(:,2)*amplitude(2)];
end

f_soln = figure; f_soln.Position = [236 58 1196 889];
subplot(221);
scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, real(phibar_prime), 'filled');
ax = gca; ax.FontSize = afs;
colorbar; colormap('jet'); pbaspect([1 2 1]);
xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
title('$Re \big( \bar{\phi}^{\prime} \big)$', 'interpreter', 'latex', 'fontsize', fs);

subplot(222);
scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, imag(phibar_prime), 'filled');
ax = gca; ax.FontSize = afs;
colorbar; colormap('jet'); pbaspect([1 2 1]);
xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
title('$Im \big( \bar{\phi}^{\prime} \big)$', 'interpreter', 'latex', 'fontsize', fs);

subplot(223);
scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, real(pbar_prime) ./ flow.pref, 'filled');
ax = gca; ax.FontSize = afs;
colorbar; colormap('jet'); pbaspect([1 2 1]);
xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
title('$Re \big( \bar{p}^{\prime} / p_{\infty} \big)$', 'interpreter', 'latex', 'fontsize', fs);
sgtitle(strcat(['Reduced Velocity Perturbation Potential and Pressure Fluctuation, $M = ', num2str(M), '$']), 'interpreter', 'latex', 'fontsize', fs);

subplot(224);
scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, imag(pbar_prime) ./ flow.pref, 'filled');
ax = gca; ax.FontSize = afs;
colorbar; colormap('jet'); pbaspect([1 2 1]);
xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
title('$Im \big( \bar{p}^{\prime} / p_{\infty} \big)$', 'interpreter', 'latex', 'fontsize', fs);
sgtitle(strcat(['Reduced Velocity Perturbation Potential and Pressure Fluctuation, $M = ', num2str(M), '$']), 'interpreter', 'latex', 'fontsize', fs);

if (dmd_on == 1)

  x_plot_dmd = dlmread('./mode1_modes_3d/x_plot_dmd.csv');
  y_plot_dmd = dlmread('./mode1_modes_3d/y_plot_dmd.csv');

  dmd_choice = 1;

  bem_mode = pbar_prime./flow.pref - p_bar_lpt;
  max_re_bem = max(real(bem_mode)); max_im_bem = max(imag(bem_mode));
  min_re_bem = min(real(bem_mode)); min_im_bem = min(imag(bem_mode));
  max_re_dmd = max(real(p_bar_dmd(:,dmd_choice))); max_im_dmd = max(imag(p_bar_dmd(:,dmd_choice)));
  min_re_dmd = min(real(p_bar_dmd(:,dmd_choice))); min_im_dmd = min(imag(p_bar_dmd(:,dmd_choice)));
  max_re = max([max_re_bem, max_re_dmd]); min_re = min([min_re_bem, min_re_dmd]);
  max_im = max([max_im_bem, max_im_dmd]); min_im = min([min_im_bem, min_im_dmd]);

  f_full_comp = figure; f_full_comp.Position = [518 -965 1110 816];
  subplot(221);
  scatter(col_pts(:,1) ./ L_panel, col_pts(:,2) ./ L_panel, sz, real(bem_mode)); grid on;
  ax = gca; ax.FontSize = afs; colorbar; colormap('jet'); caxis([min_re, max_re]); pbaspect([1 2 1]);
  ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
  title('Re( $\bar{p}^{\prime}_{BEM}$ ) / $p_{\infty}$ ', 'interpreter', 'latex', 'fontsize', fs);

  subplot(222);
  scatter(col_pts(:,1) ./ L_panel, col_pts(:,2) ./ L_panel, sz, imag(bem_mode)); grid on;
  ax = gca; ax.FontSize = afs; colorbar; colormap('jet'); caxis([min_im, max_im]); pbaspect([1 2 1]);
  ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
  title('Im( $\bar{p}^{\prime}_{BEM}$ ) / $p_{\infty}$ ', 'interpreter', 'latex', 'fontsize', fs);

  subplot(223);
  scatter(x_plot_dmd, y_plot_dmd, sz, real(p_bar_dmd(:,dmd_choice))); grid on;
  ax = gca; ax.FontSize = afs; colorbar; colormap('jet'); caxis([min_re, max_re]); pbaspect([1 2 1]);
  xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
  ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
  title('Re( $\bar{p}^{\prime}_{DMD}$ ) / $p_{\infty}$', 'interpreter', 'latex', 'fontsize', fs);

  subplot(224);
  scatter(x_plot_dmd, y_plot_dmd, sz, imag(p_bar_dmd(:,dmd_choice))); grid on;
  ax = gca; ax.FontSize = afs; colorbar; colormap('jet'); caxis([min_im, max_im]); pbaspect([1 2 1]);
  ylabel('Im( $\bar{p}^{\prime}_{DMD}$ ) / $p_{\infty}$', 'interpreter', 'latex', 'fontsize', fs);
  xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
  title('Im( $\bar{p}^{\prime}_{DMD}$ ) / $p_{\infty}$', 'interpreter', 'latex', 'fontsize', fs);

  sgtitle(strcat(['Comparison of DMD Mode to BEM, $M = ', num2str(flow.M), '$']), 'interpreter', 'latex', 'fontsize', fs);
  
  % Plot mid-span results
  range = (floor((Ny-1)/2)*(Nx-1)+1):(floor((Ny-1)/2)*(Nx-1) + (Nx-1));
  dmd_range = (floor(41/2)*21+1):(floor(41/2)*21+21);
  figure;
  subplot(211);
  plot(col_pts(range,1) ./ L_panel, real(bem_mode(range)), 'b', 'LineWidth', 2); grid on; hold on;
  plot(x_plot_dmd(dmd_range), real(p_bar_dmd(dmd_range,dmd_choice)), 'r', 'LineWidth', 2);
  ax = gca; ax.FontSize = afs;
  legend('BEM', 'DMD', 'interpreter', 'latex', 'fontsize', fs, 'location', 'northwest');
  ylabel('Re( $\bar{p}^{\prime} / p_{\infty}$ )', 'interpreter', 'latex', 'fontsize', fs);
  subplot(212);
  plot(col_pts(range,1) ./ L_panel, imag(bem_mode(range)), 'b', 'LineWidth', 2); grid on; hold on;
  plot(x_plot_dmd(dmd_range), imag(p_bar_dmd(dmd_range,dmd_choice)), 'r', 'LineWidth', 2);
  ax = gca; ax.FontSize = afs;
  ylabel('Im( $\bar{p}^{\prime} / p_{\infty}$ )', 'interpreter', 'latex', 'fontsize', fs);
  xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
  sgtitle(strcat(['BEM-DMD Comparison, Half-Span Plane, $M = ', num2str(flow.M), '$']), 'interpreter', 'latex', 'fontsize', fs);

  pause;
end


% %% Temp -- plot just BEM
% snap_level = 1; snapshot_begin = 1; snapshot_end = 401;
% snapshot_vec = snapshot_begin:snap_level:snapshot_end;
% pc2_dt_nondim = 0.001;
% t_scale = 1 / sqrt(flow.pref / flow.rhoref);
% snapshot = 20*snap_level;
% dmd_dt = (pc2_dt_nondim*t_scale) * snapshot;
% 
% t_phys = (1:(snapshot_end-snapshot_begin+1)) .* dmd_dt; t_phys = t_phys - dmd_dt;
% maxT = max(size(t_phys));
% 
% max_pf_bem = max(max( imag( (pbar_prime') * exp(i*omega*t_phys) )./flow.pref ));
% min_pf_bem = min(min( imag( (pbar_prime') * exp(i*omega*t_phys) )./flow.pref ));
% 
% f_comp = figure; f_comp.Position = [189 152 1263 630];
% v_fname = strcat(['./video_results/3d/m', m_choice, '.mp4']);
% v_comp = VideoWriter(v_fname, 'MPEG-4');
% v_comp.FrameRate = 20; v_comp.Quality = 10;
% open(v_comp);
% for l = 1:length(t_phys)
%   t_loc = t_phys(l);
%   p_loc = imag( pbar_prime*exp(i*omega*t_loc) );
%   z_loc = imag( psi.*(amp*exp(i*omega*t_loc)) );
% 
%   subplot(121);
%   scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, p_loc ./ flow.pref, 'filled'); grid on;
%   colorbar; colormap('jet'); caxis([min_pf_bem, max_pf_bem]); pbaspect([1 2 1]);
%   ax = gca; ax.FontSize = afs;
%   xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
%   ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
%   title('BEM', 'interpreter', 'latex', 'fontsize', fs);
% 
%   subplot(122);
%   scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, z_loc, 'filled'); grid on;
%   colorbar; colormap('jet'); caxis([-amp,amp]); pbaspect([1 2 1]);
%   ax = gca; ax.FontSize = afs;
%   xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
%   ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
%   title('$w_{BEM} (x,y)$', 'interpreter', 'latex', 'fontsize', fs);
% 
%   sgtitle(strcat(['BEM, $M = ', num2str(M), ', t = ', num2str(l), '/', num2str(maxT), '$']), 'interpreter', 'latex', 'fontsize', fs);
%   
%   frame = getframe(gcf);
%   writeVideo(v_comp, frame);
%   clf;
% 
% end
% close(v_comp);


%% Comparison to PC2 simulation
disp('Preparing to begin video comparison to PC2. Press any button to begin.'); pause;

snap_level = 1; snapshot_begin = 1; snapshot_end = 401;
snapshot_vec = snapshot_begin:snap_level:snapshot_end;
pc2_dt_nondim = 0.001;
t_scale = 1 / sqrt(flow.pref / flow.rhoref);
snapshot = 20*snap_level;
dmd_dt = (pc2_dt_nondim*t_scale) * snapshot;

t_phys = (1:(snapshot_end-snapshot_begin+1)) .* dmd_dt; t_phys = t_phys - dmd_dt;
maxT = max(size(t_phys));

% Just single mode 1.1 \leq M \leq 1.4
pressure_data = dlmread(strcat(['./comparison_data/3d/panel_runs_supersonic/m', m_choice, '/surface_pressure_history_m', m_choice, '.csv']));
y_data = dlmread(strcat(['./comparison_data/3d/panel_runs_supersonic/m', m_choice, '/surface_yloc_history_m', m_choice, '.csv']));
x_data = dlmread(strcat(['./comparison_data/3d/panel_runs_supersonic/m', m_choice, '/surface_xloc_history_m', m_choice, '.csv']));
z_data = dlmread(strcat(['./comparison_data/3d/panel_runs_supersonic/m', m_choice, '/surface_zloc_history_m', m_choice, '.csv']));

% Full set of supersonic data
% pressure_data = dlmread(strcat(['./comparison_data/3d/supersonic_data/m', m_choice, '/surface_pressure_history_m', m_choice, '.csv']));
% y_data = dlmread(strcat(['./comparison_data/3d/supersonic_data/m', m_choice, '/surface_yloc_history_m', m_choice, '.csv']));
% x_data = dlmread(strcat(['./comparison_data/3d/supersonic_data/m', m_choice, '/surface_xloc_history_m', m_choice, '.csv']));
% z_data = dlmread(strcat(['./comparison_data/3d/supersonic_data/m', m_choice, '/surface_zloc_history_m', m_choice, '.csv']));

% M=3 single
% pressure_data = dlmread(strcat(['./comparison_data/3d/supersonic_data/m', m_choice, '_single/surface_pressure_history_m', m_choice, '_single.csv']));
% y_data = dlmread(strcat(['./comparison_data/3d/supersonic_data/m', m_choice, '_single/surface_yloc_history_m', m_choice, '_single.csv']));
% x_data = dlmread(strcat(['./comparison_data/3d/supersonic_data/m', m_choice, '_single/surface_xloc_history_m', m_choice, '_single.csv']));
% z_data = dlmread(strcat(['./comparison_data/3d/supersonic_data/m', m_choice, '_single/surface_zloc_history_m', m_choice, '_single.csv']));


pressure_data = pressure_data';
x_data = x_data';
y_data = y_data'; pd_y_locs = y_data(:,1); pd_y_locs = pd_y_locs - pd_y_locs(1);
z_data = z_data'; pd_x_locs = x_data(:,1); pd_x_locs = pd_x_locs - pd_x_locs(1);

max_pf_bem = max(max( ((pbar_prime')*exp(i*omega*t_phys) + conj((pbar_prime')*exp(i*omega*t_phys)))./flow.pref ));
min_pf_bem = min(min( ((pbar_prime')*exp(i*omega*t_phys) + conj((pbar_prime')*exp(i*omega*t_phys)))./flow.pref ));

% max_pf_bem = max(max( imag( (pbar_prime') * exp(i*omega*t_phys) )./flow.pref ));
% min_pf_bem = min(min( imag( (pbar_prime') * exp(i*omega*t_phys) )./flow.pref ));
max_pf_cfd = max(max( pressure_data - 1 ));
min_pf_cfd = min(min( pressure_data - 1 ));

disp(strcat(['max pf bem: ', num2str(max_pf_bem)]));
disp(strcat(['max pf cfd: ', num2str(max_pf_cfd)]));
diff = ((max_pf_bem - max_pf_cfd)/max_pf_cfd)*100;
disp(strcat(['difference: ', num2str(diff)]));
pause;

f_comp = figure; f_comp.Position = [189 152 1263 630];
v_fname = strcat(['./video_results/3d/m', m_choice, '_comparison.mp4']);
v_comp = VideoWriter(v_fname, 'MPEG-4');
v_comp.FrameRate = 20; v_comp.Quality = 10;
open(v_comp);
for l = 1:length(t_phys)
  t_loc = t_phys(l);
  p_loc = pbar_prime*exp(i*omega*t_loc) + conj(pbar_prime*exp(i*omega*t_loc));
  % p_loc = imag( pbar_prime*exp(i*omega*t_loc) );
  z_loc = imag( psi.*(amp*exp(i*omega*t_loc)) );

  subplot(141);
  scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, p_loc ./ flow.pref, 'filled'); grid on;
  colorbar; colormap('jet'); caxis([min_pf_bem, max_pf_bem]); pbaspect([1 2 1]);
  ax = gca; ax.FontSize = afs;
  xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
  ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
  title('BEM', 'interpreter', 'latex', 'fontsize', fs);

  subplot(142);
  scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, z_loc, 'filled'); grid on;
  colorbar; colormap('jet'); caxis([-amp,amp]); pbaspect([1 2 1]);
  ax = gca; ax.FontSize = afs;
  xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
  ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
  title('$w_{BEM} (x,y)$', 'interpreter', 'latex', 'fontsize', fs);

  subplot(143);
  scatter(pd_x_locs ./ La, pd_y_locs ./ La, sz, pressure_data(:,l) - 1, 'filled'); grid on;
  colorbar; colormap('jet'); caxis([min_pf_cfd, max_pf_cfd]); pbaspect([1 2 1]);
  ax = gca; ax.FontSize = afs;
  xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
  ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
  title('CFD', 'interpreter', 'latex', 'fontsize', fs);

  subplot(144);
  scatter(pd_x_locs ./ La, pd_y_locs ./ La, sz, z_data(:,l), 'filled'); grid on;
  colorbar; colormap('jet'); caxis([-amp,amp]); pbaspect([1 2 1]);
  ax = gca; ax.FontSize = afs;
  xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
  ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
  title('$w_{CFD} (x,y)$', 'interpreter', 'latex', 'fontsize', fs);

  sgtitle(strcat(['BEM vs. CFD, $M = ', num2str(M), ', t = ', num2str(l), '/', num2str(maxT), '$']), 'interpreter', 'latex', 'fontsize', fs);
  
  frame = getframe(gcf);
  writeVideo(v_comp, frame);
  clf;

end
close(v_comp);






%% Internal functions
function [belong_elem, fluid_node_bary] = find_struct_face(fluid_node_loc, surface_node_locs, structural_node_loc, structural_node_elem, structural_elem_node, structural_faces, neg_epsilon)

  % find clsest structural node
  node_dist = fluid_node_loc' - surface_node_locs(:,2:4)';
  node_dist = vecnorm(node_dist);
  [node_dist, sort_ind] = sort(node_dist);

  done = 0;
  for k = 1:length(node_dist)
    closest_node_ind = surface_node_locs(sort_ind(k), 1);

    % go through faces that share node
    elems_containing_node = structural_node_elem(closest_node_ind, :);
    %disp(elems_containing_node)
    for j = 1:length(elems_containing_node)
      % grab an element that contains closest node
      loc_elem = elems_containing_node(j);
  
      % find surface face that belongs to this element
      face_elem = find(structural_faces(:,1) == loc_elem);
      if isempty(face_elem)
        continue;
      end
  
      % get the nodes that comprise the actual surface face belonging to this element
      face_nodes = structural_faces(face_elem, 2:11);
  
      % Identify which corner nodes are on the surface face
      loc_elem_surf_num = find(face_nodes(1:4) == 0);
      loc_elem_surf = [1 2 3 4]; loc_elem_surf(loc_elem_surf_num) = [];
  
      % Backsolve for area coordinates -- in quadratic space
      coords = bary_coords(fluid_node_loc, loc_elem, structural_elem_node, structural_node_loc, loc_elem_surf_num);
      %disp(coords) 
      if ~((coords(1) < neg_epsilon) || (coords(2) < neg_epsilon) || (coords(3) < neg_epsilon) || (coords(4) < neg_epsilon))
        belong_elem = loc_elem;
        fluid_node_bary = coords;
        done = 1;
        break;
      end
    end

    if (done == 1)
      break;
    end
  end

end

function [coordinates] = bary_coords(fluid_node_loc, loc_struc_elem, elem_node, structural_node_loc, loc_elem_surf_num)

  loc_elem_nodes = elem_node(loc_struc_elem, :);
  loc_elem_nodelocs = structural_node_loc(loc_elem_nodes, :);
  fun = @(L) bary_eqns(L, loc_elem_nodelocs, fluid_node_loc);
  if (loc_elem_surf_num == 1)
    L0 = [0, 1/3, 1/3, 1/3];
  elseif (loc_elem_surf_num == 2)
    L0 = [1/3, 0, 1/3, 1/3];
  elseif (loc_elem_surf_num == 3)
    L0 = [1/3, 1/3, 0, 1/3];
  else
    L0 = [1/3, 1/3, 1/3, 0];
  end
  options = optimoptions('fsolve','Display','off');
  coordinates = fsolve(fun, L0, options);

end

function F = bary_eqns(L, loc_elem_nodelocs, fluid_node_loc)

  loc_elem_x = loc_elem_nodelocs(:,1); loc_elem_y = loc_elem_nodelocs(:,2); loc_elem_z = loc_elem_nodelocs(:,3);

  F(1) = loc_elem_x(1)*(2*L(1)-1)*L(1) + loc_elem_x(2)*(2*L(2)-1)*L(2) + loc_elem_x(3)*(2*L(3)-1)*L(3) + loc_elem_x(4)*(2*L(4)-1)*L(4) + ...
               4*loc_elem_x(5)*L(1)*L(2) + 4*loc_elem_x(6)*L(2)*L(3) + 4*loc_elem_x(7)*L(1)*L(3) + ...
               4*loc_elem_x(8)*L(1)*L(4) + 4*loc_elem_x(9)*L(2)*L(4) + 4*loc_elem_x(10)*L(3)*L(4) - fluid_node_loc(1);

  F(2) = loc_elem_y(1)*(2*L(1)-1)*L(1) + loc_elem_y(2)*(2*L(2)-1)*L(2) + loc_elem_y(3)*(2*L(3)-1)*L(3) + loc_elem_y(4)*(2*L(4)-1)*L(4) + ...
               4*loc_elem_y(5)*L(1)*L(2) + 4*loc_elem_y(6)*L(2)*L(3) + 4*loc_elem_y(7)*L(1)*L(3) + ...
               4*loc_elem_y(8)*L(1)*L(4) + 4*loc_elem_y(9)*L(2)*L(4) + 4*loc_elem_y(10)*L(3)*L(4) - fluid_node_loc(2);

  F(3) = loc_elem_z(1)*(2*L(1)-1)*L(1) + loc_elem_z(2)*(2*L(2)-1)*L(2) + loc_elem_z(3)*(2*L(3)-1)*L(3) + loc_elem_z(4)*(2*L(4)-1)*L(4) + ...
                4*loc_elem_z(5)*L(1)*L(2) + 4*loc_elem_z(6)*L(2)*L(3) + 4*loc_elem_z(7)*L(1)*L(3) + ...
                4*loc_elem_z(8)*L(1)*L(4) + 4*loc_elem_z(9)*L(2)*L(4) + 4*loc_elem_z(10)*L(3)*L(4) - fluid_node_loc(3);

  F(4) = L(1) + L(2) + L(3) + L(4) - 1;


end

function [interp_val] = bary_tet_interp(L, val)

  interp_val = val(1)*(2*L(1)-1)*L(1) + val(2)*(2*L(2)-1)*L(2) + val(3)*(2*L(3)-1)*L(3) + val(4)*(2*L(4)-1)*L(4) + ...
               4*val(5)*L(1)*L(2) + 4*val(6)*L(2)*L(3) + 4*val(7)*L(1)*L(3) + ...
               4*val(8)*L(1)*L(4) + 4*val(9)*L(2)*L(4) + 4*val(10)*L(3)*L(4);

end

