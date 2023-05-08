
close all;
clear all;
clc;

addpath('./utilities/3d/');
addpath('./utilities/2d/');

fs = 20;
afs = 18;
lfs = 20;
sz = 36;

% Define problem parameters
m_choice = '11'; amp_choice = 'origamp';
M = str2num(strcat([m_choice(1), '.', m_choice(2:end)]));
amp = 0.0000508; % amp = 0.000722157;
omega = 691.15038379;
dmd_on = 1;

run_3danal = 1;
run_3dstanal = 0;

flow.pref = 34601.89642944;
flow.rhoref = 0.49774215;
flow.gamma = 1.4;
flow.Tref = 242.222;
flow.M = M;

if (M >= 1.4)
  flow.pref = 41596.29322233;
  flow.rhoref = 0.64847794;
  flow.Tref = 223.5;
end

R = flow.pref / (flow.rhoref * flow.Tref);
k_wav = omega / sqrt(flow.gamma * R * flow.Tref);

%% Define deformation function -- C6 geometry
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

p = A\b; p = inv(A)*b;
p = [flip(p); zeros(7,1)];

pd = 14:-1:1; pd = pd.*p(1:end-1)';

%% Define BEM collocation and boundary points
La = 0.2286; Lb = 2*La;
flow.L = La; flow.La = La;
Nx = 72; Ny = 502;
%Nx = 62; Ny = 452; % resolution for thesis images
N = (Nx-1)*(Ny-1);
col_pts = zeros(N, 2);  bndry_pts = zeros(4*N, 2);
psi_in = zeros(N,1); dpsi_dx_in = zeros(N,1);
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

    psi_loc_y = polyval(p, col_y(k) ./ Lb);
    psi_in(m) = psi_loc_y * polyval(p, col_x(j) ./ La);
    dpsi_dx_in(m) = psi_loc_y * polyval(pd, col_x(j) ./ La) ./ La;

  end
end

% Define geometry function for pass to BEM 3D suite
psi_func = @(x,y,La,Lb) ...
(p(1).*(y./Lb).^14 + p(2).*(y./Lb).^13 + p(3).*(y./Lb).^12 + p(4).*(y./Lb).^11 + p(5).*(y./Lb).^10 + p(6).*(y./Lb).^9 + p(7).*(y./Lb).^8 + p(8).*(y./Lb).^7) .* ...
(p(1).*(x./La).^14 + p(2).*(x./La).^13 + p(3).*(x./La).^12 + p(4).*(x./La).^11 + p(5).*(x./La).^10 + p(6).*(x./La).^9 + p(7).*(x./La).^8 + p(8).*(x./La).^7);
%dpsi_func = @(x,y,La,Lb) ...
%(p(1).*(y./Lb).^14 + p(2).*(y./Lb).^13 + p(3).*(y./Lb).^12 + p(4).*(y./Lb).^11 + p(5).*(y./Lb).^10 + p(6).*(y./Lb).^9 + p(7).*(y./Lb).^8 + p(8).*(y./Lb).^7) .* ...
%(pd(1).*(x./La).^13 + pd(2).*(x./La).^12 + pd(3).*(x./La).^11 + pd(4).*(x./La).^10 + pd(5).*(x./La).^9 + pd(6).*(x./La).^8 + pd(7).*(x./La).^7 + pd(8).*(x./La).^6) ./ La;

dpsi_func = @(x,y,La,Lb) ...
(p(1).*(y./Lb).^14 + p(2).*(y./Lb).^13 + p(3).*(y./Lb).^12 + p(4).*(y./Lb).^11 + p(5).*(y./Lb).^10 + p(6).*(y./Lb).^9 + p(7).*(y./Lb).^8 + p(8).*(y./Lb).^7) .* ...
((14./(La^14)).*p(1).*x.^13 + (13./(La.^13)).*p(2).*x.^12 + (12./(La^12)).*p(3).*x.^11 + (11./(La^11)).*p(4).*x.^10 + (10./(La.^10)).*p(5).*x.^9 + (9./(La^9)).*p(6).*x.^8 + (8./(La^8)).*p(7).*x.^7 + (7./(La^7)).*p(8).*x.^6);

%(pd(1).*(x./La).^13 + pd(2).*(x./La).^12 + pd(3).*(x./La).^11 + pd(4).*(x./La).^10 + pd(5).*(x./La).^9 + pd(6).*(x./La).^8 + pd(7).*(x./La).^7 + pd(8).*(x./La).^6) ./ La;


flow.Lb = Lb; flow.psi = psi_func; flow.dpsi = dpsi_func;

figure;
scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, 'r'); hold on; grid on;
scatter(bndry_pts(:,1) ./ La, bndry_pts(:,2) ./ La, 'b');
pbaspect([1 2 1]);
xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
title('Collocation and Element Boundary Points', 'interpreter', 'latex', 'fontsize', fs);

figure;
scatter3(col_pts(:,1) ./ La, col_pts(:,2) ./ La, flow.psi(col_pts(:,1), col_pts(:,2), La, Lb), sz, flow.psi(col_pts(:,1), col_pts(:,2), La, Lb)); grid on;
colorbar; colormap('jet'); pbaspect([1 2 1]); view(2);
ax = gca; ax.FontSize = afs;
xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
zlabel('$\psi$', 'interpreter', 'latex', 'fontsize', fs);
title('$\psi (x,y)$', 'interpreter', 'latex', 'fontsize', fs);

figure;
scatter3(col_pts(:,1) ./ La, col_pts(:,2) ./ La, flow.dpsi(col_pts(:,1), col_pts(:,2), La, Lb), sz, flow.dpsi(col_pts(:,1), col_pts(:,2), La, Lb)); grid on;
colorbar; colormap('jet'); pbaspect([1 2 1]); view(2);
ax = gca; ax.FontSize = afs;
xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
zlabel('$\partial \psi / \partial x$', 'interpreter', 'latex', 'fontsize', fs);
title('$\partial \psi(x,y) / \partial x$', 'interpreter', 'latex', 'fontsize', fs);

%% Compute spatial solution of reduced velocity perturbation potential and reduced pressure fluctuation
%% using 3D BEM
[phibar_prime_3d_anal, pbar_prime_3d_anal] = bem3d(col_pts, bndry_pts, Nx-1, psi_in, dpsi_dx_in, amp, omega, flow);

% pbar_prime_3d_anal_mid = zeros(Nx-1, 1); col_pts_mid = zeros(Nx-1,1);
% for j = 1:Nx-1
%   ind = (61-1)*(Nx-1) + j;
%   pbar_prime_3d_anal_mid(j) = pbar_prime_3d_anal(ind);
%   col_pts_mid(j) = col_pts(ind);
% end
% figure;
% plot(col_pts_mid ./ La, real(pbar_prime_3d_anal_mid) ./ flow.pref, 'b', 'LineWidth', 2); grid on; hold on;
% plot(col_pts_mid ./ La, imag(pbar_prime_3d_anal_mid) ./ flow.pref, '--k', 'LineWidth', 2);
% ylim([-0.01, 0.02]);
% ax = gca; ax.FontSize = afs;
% xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
% ylabel('$\bar{p}^{\prime}_{BEM}$', 'interpreter', 'latex', 'fontsize', fs);
% title('Pressure Solution, 3D C6-Continuous Geometry', 'interpreter', 'latex', 'fontsize', fs);
% pause;
% dlmwrite('./m085_c6_temp.csv', pbase_prime_3d_anal_mid ./ flow.pref, 'delimiter', ',', 'precision', 12);

%% Plot all solutions

%% 3D Analytical Geometry
if (run_3danal == 1)
  f_soln_3d_anal = figure; f_soln_3d_anal.Position = [236 58 1196 889];
  subplot(221);
  scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, real(phibar_prime_3d_anal), 'filled');
  ax = gca; ax.FontSize = afs;
  colorbar; colormap('jet'); pbaspect([1 2 1]);
  xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
  ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
  title('$Re \big( \bar{\phi}^{\prime} \big)$', 'interpreter', 'latex', 'fontsize', fs);
  
  subplot(222);
  scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, imag(phibar_prime_3d_anal), 'filled');
  ax = gca; ax.FontSize = afs;
  colorbar; colormap('jet'); pbaspect([1 2 1]);
  xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
  ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
  title('$Im \big( \bar{\phi}^{\prime} \big)$', 'interpreter', 'latex', 'fontsize', fs);
  
  subplot(223);
  scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, real(pbar_prime_3d_anal) ./ flow.pref, 'filled');
  ax = gca; ax.FontSize = afs;
  colorbar; colormap('jet'); pbaspect([1 2 1]);
  xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
  ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
  title('$Re \big( \bar{p}^{\prime} / p_{\infty} \big)$', 'interpreter', 'latex', 'fontsize', fs);
  
  subplot(224);
  scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, imag(pbar_prime_3d_anal) ./ flow.pref, 'filled');
  ax = gca; ax.FontSize = afs;
  colorbar; colormap('jet'); pbaspect([1 2 1]);
  xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
  ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
  title('$Im \big( \bar{p}^{\prime} / p_{\infty} \big)$', 'interpreter', 'latex', 'fontsize', fs);
  sgtitle(strcat(['Reduced Velocity Perturbation Potential and Pressure Fluctuation, 3D Analytical Geom, $M = ', num2str(M), '$']), 'interpreter', 'latex', 'fontsize', fs);
end



%% Comparison to DMD

% Compute LPT mode
p_bar_lpt = (-i*(amp/2)*flow.gamma*flow.M.*dpsi_dx_in + flow.gamma*k_wav*(amp/2).*psi_in);

figure;
subplot(121);
scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, real(p_bar_lpt)); grid on;
colorbar; colormap('jet');
subplot(122);
scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, imag(p_bar_lpt)); grid on;
colorbar; colormap('jet');

if (dmd_on == 1)
  p_bar_dmd = dlmread(strcat(['./mode1_modes_3d/c6/m', m_choice, '_modes_', amp_choice, '.csv']));
  amplitude = dlmread(strcat(['./mode1_modes_3d/c6/m', m_choice, '_amps_', amp_choice, '.csv']));
  p_bar_dmd = [p_bar_dmd(:,1)*amplitude(1), p_bar_dmd(:,2)*amplitude(2)];
end

if (dmd_on == 1)

  x_plot_dmd = dlmread('./mode1_modes_3d/x_plot_dmd.csv');
  y_plot_dmd = dlmread('./mode1_modes_3d/y_plot_dmd.csv');

  bem_mode = (pbar_prime_3d_anal.')./flow.pref - p_bar_lpt;
  %bem_mode = (pbar_st.') ./ flow.pref - p_bar_lpt;
  dmd_choice = 1;

  f_full_comp = figure; f_full_comp.Position = [518 -965 1110 816];
  subplot(221);
  scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, real(bem_mode), 'filled'); grid on;
  ax = gca; ax.FontSize = afs; colorbar; colormap('jet'); pbaspect([1 2 1]);
  ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
  title('Re( $\bar{p}^{\prime}_{BEM}$ ) / $p_{\infty}$ ', 'interpreter', 'latex', 'fontsize', fs);

  subplot(222);
  scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, imag(bem_mode), 'filled'); grid on;
  ax = gca; ax.FontSize = afs; colorbar; colormap('jet'); pbaspect([1 2 1]);
  ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
  title('Im( $\bar{p}^{\prime}_{BEM}$ ) / $p_{\infty}$ ', 'interpreter', 'latex', 'fontsize', fs);

  subplot(223);
  scatter(x_plot_dmd, y_plot_dmd, sz, real(p_bar_dmd(:,dmd_choice)), 'filled'); grid on;
  ax = gca; ax.FontSize = afs; colorbar; colormap('jet'); pbaspect([1 2 1]);
  xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
  ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
  title('Re( $\bar{p}^{\prime}_{DMD}$ ) / $p_{\infty}$', 'interpreter', 'latex', 'fontsize', fs);

  subplot(224);
  scatter(x_plot_dmd, y_plot_dmd, sz, imag(p_bar_dmd(:,dmd_choice)), 'filled'); grid on;
  ax = gca; ax.FontSize = afs; colorbar; colormap('jet'); pbaspect([1 2 1]);
  ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
  xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
  title('Im( $\bar{p}^{\prime}_{DMD}$ ) / $p_{\infty}$', 'interpreter', 'latex', 'fontsize', fs);

  sgtitle(strcat(['Comparison of DMD Mode to BEM, $M = ', num2str(flow.M), '$']), 'interpreter', 'latex', 'fontsize', fs);
  

  % Plot mid-span results
  range = (floor((Ny-1)/2)*(Nx-1)+1):(floor((Ny-1)/2)*(Nx-1) + (Nx-1));
  %dmd_range = (floor(41/2)*21+1):(floor(41/2)*21+21);
  dmd_range = (floor(41/2)*41+1):(floor(41/2)*41+41);
  f_span = figure; f_span.Position = [560 478 847 469];
  subplot(211);
  plot(col_pts(range,1) ./ La, real(bem_mode(range)), 'b', 'LineWidth', 2); grid on; hold on;
  plot(x_plot_dmd(dmd_range), real(p_bar_dmd(dmd_range,dmd_choice)), 'r', 'LineWidth', 2);
  ax = gca; ax.FontSize = afs;
  legend('BEM', 'DMD', 'interpreter', 'latex', 'fontsize', fs, 'location', 'southwest');
  ylabel('Re( $\bar{p}^{\prime} / p_{\infty}$ )', 'interpreter', 'latex', 'fontsize', fs);
  subplot(212);
  plot(col_pts(range,1) ./ La, imag(bem_mode(range)), 'b', 'LineWidth', 2); grid on; hold on;
  plot(x_plot_dmd(dmd_range), imag(p_bar_dmd(dmd_range,dmd_choice)), 'r', 'LineWidth', 2);
  ax = gca; ax.FontSize = afs;
  ylabel('Im( $\bar{p}^{\prime} / p_{\infty}$ )', 'interpreter', 'latex', 'fontsize', fs);
  xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
  sgtitle(strcat(['BEM-DMD Comparison, Half-Span Plane, $M = ', num2str(flow.M), '$']), 'interpreter', 'latex', 'fontsize', fs);

  % Temporarily write results
  dlmwrite('~/research/dmdc_panel/pc2_varying_cases/dmd_results/BEM/xrange.csv', col_pts(range,1), 'delimiter', ',', 'precision', 12);
  dlmwrite('~/research/dmdc_panel/pc2_varying_cases/dmd_results/BEM/bem_mode.csv', bem_mode(range), 'delimiter', ',', 'precision', 12);

end


%% Comparisons to PC2 simulation
disp('Preparing to begin video comparison to PC2. Press any button to begin.'); %pause;

snap_level = 1; snapshot_begin = 1; snapshot_end = 361;
snapshot_vec = snapshot_begin:snap_level:snapshot_end;
if (M == 1.1)
  pc2_dt_nondim = 0.0009987201549;
elseif (M == 1.5)
  pc2_dt_nondim = 0.0009593471401;
else
  pc2_dt_nondim = 0.001;
end
t_scale = 1 / sqrt(flow.pref / flow.rhoref);
snapshot = 20*snap_level;
dmd_dt = (pc2_dt_nondim*t_scale) * snapshot;

t_phys = (1:(snapshot_end-snapshot_begin+1)) .* dmd_dt; t_phys = t_phys - dmd_dt;
maxT = max(size(t_phys));

if (flow.M > 1)
  pressure_data = dlmread(strcat(['./comparison_data/c6_3d/c6_3d_supersonic_data/m', m_choice, '/', amp_choice, '/surface_pressure_history_c6_m', m_choice, '_', amp_choice, '.csv']));
  x_data = dlmread(strcat(['./comparison_data/c6_3d/c6_3d_supersonic_data/m', m_choice, '/', amp_choice, '/surface_xloc_history_c6_m', m_choice, '_', amp_choice, '.csv']));
  y_data = dlmread(strcat(['./comparison_data/c6_3d/c6_3d_supersonic_data/m', m_choice, '/', amp_choice, '/surface_yloc_history_c6_m', m_choice, '_', amp_choice, '.csv']));
  z_data = dlmread(strcat(['./comparison_data/c6_3d/c6_3d_supersonic_data/m', m_choice, '/', amp_choice, '/surface_zloc_history_c6_m', m_choice, '_', amp_choice, '.csv']));
else
  pressure_data = dlmread(strcat(['./comparison_data/c6_3d/c6_3d_subsonic_data/m', m_choice, '/', amp_choice, '/surface_pressure_history_c6_m', m_choice, '_', amp_choice, '.csv']));
  x_data = dlmread(strcat(['./comparison_data/c6_3d/c6_3d_subsonic_data/m', m_choice, '/', amp_choice, '/surface_xloc_history_c6_m', m_choice, '_', amp_choice, '.csv']));
  y_data = dlmread(strcat(['./comparison_data/c6_3d/c6_3d_subsonic_data/m', m_choice, '/', amp_choice, '/surface_yloc_history_c6_m', m_choice, '_', amp_choice, '.csv']));
  z_data = dlmread(strcat(['./comparison_data/c6_3d/c6_3d_subsonic_data/m', m_choice, '/', amp_choice, '/surface_zloc_history_c6_m', m_choice, '_', amp_choice, '.csv']));
end
pressure_data = pressure_data';
x_data = x_data';
y_data = y_data'; pd_y_locs = y_data(:,1); pd_y_locs = pd_y_locs - pd_y_locs(1);
z_data = z_data'; pd_x_locs = x_data(:,1); pd_x_locs = pd_x_locs - pd_x_locs(1);

% Compute max/min pressure fluctuations from BEM approaches

pf_3d_bem_history = (pbar_prime_3d_anal.')*exp(i*omega*t_phys); pf_3d_bem_history = (pf_3d_bem_history + conj(pf_3d_bem_history)) ./ flow.pref;

max_pf_bem_3dgeom = max(max(  pf_3d_bem_history ));
min_pf_bem_3dgeom = min(min(  pf_3d_bem_history ));

max_pf_cfd = max(max( pressure_data - 1 ));
min_pf_cfd = min(min( pressure_data - 1 ));

f_cfd_comp = figure; f_cfd_comp.Position = [639 -815 722 475];
max_pf = max([max_pf_bem_3dgeom, max_pf_cfd]); min_pf = min([min_pf_bem_3dgeom, min_pf_cfd]);
pf_bem_history = (pbar_prime_3d_anal.') * exp(i*omega*t_phys); pf_bem_history = pf_bem_history + conj(pf_bem_history);
subplot(121);
scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, pf_bem_history(:,141) ./ flow.pref, 'filled'); grid on;
colorbar; colormap('jet'); caxis([min_pf, max_pf]); pbaspect([1 2 1]);
ax = gca; ax.FontSize = afs;
xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
title('BEM', 'interpreter', 'latex', 'fontsize', fs);
subplot(122);
scatter(pd_x_locs ./ La, pd_y_locs ./ La, sz, pressure_data(:,141) - 1, 'filled'); grid on;
colorbar; colormap('jet'); caxis([min_pf, max_pf]); pbaspect([1 2 1]);
ax = gca; ax.FontSize = afs;
xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
title('CFD', 'interpreter', 'latex', 'fontsize', fs);
sgtitle(strcat(['BEM-CFD Comparison, M = $', num2str(flow.M), '$']), 'interpreter', 'latex', 'fontsize', fs);

pause;
if (run_3danal == 1)
  f_comp_3danal = figure; f_comp_3danal.Position = [189 152 1263 630];
  v_fname = strcat(['./video_results/3d/m', m_choice, '_3danal_comparison_81.mp4']);
  v_comp_3danal = VideoWriter(v_fname, 'MPEG-4');
  v_comp_3danal.FrameRate = 20; v_comp_3danal.Quality = 10;
  open(v_comp_3danal);
  for l = 1:length(t_phys)
    t_loc = t_phys(l);
    p_loc = pbar_prime_3d_anal*exp(i*omega*t_loc); p_loc = p_loc + conj(p_loc);
    z_loc = (-i*amp/2).*psi_in.*exp(i*omega*t_loc); z_loc = z_loc + conj(z_loc);
  
    subplot(141);
    scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, p_loc ./ flow.pref, 'filled'); grid on;
    colorbar; colormap('jet'); %caxis([min_pf_bem_3dgeom, max_pf_bem_3dgeom]); pbaspect([1 2 1]);
    caxis([min_pf_bem_3dgeom, max_pf_bem_3dgeom]); pbaspect([1 2 1]);
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
  
    sgtitle(strcat(['BEM vs. CFD, 3D Analytical Geom., $M = ', num2str(M), ', t = ', num2str(l), '/', num2str(maxT), '$']), 'interpreter', 'latex', 'fontsize', fs);
    
    frame = getframe(gcf);
    writeVideo(v_comp_3danal, frame);
    clf;
  
  end
  close(v_comp_3danal);


  f_comp_3danal_halfspan = figure; f_comp_3danal_halfspan.Position = [189 152 1263 630];
  v_fname = strcat(['./video_results/3d/m', m_choice, '_3danal_halfspan_comparison_81.mp4']);
  v_comp_3danal_halfspan = VideoWriter(v_fname, 'MPEG-4');
  v_comp_3danal_halfspan.FrameRate = 20; v_comp_3danal_halfspan.Quality = 10;

  bem_range = (floor((Ny-1)/2)*(Nx-1)+1):(floor((Ny-1)/2)*(Nx-1) + (Nx-1));
  dmd_range = (floor(41/2)*41+1):(floor(41/2)*41+41);

  open(v_comp_3danal_halfspan);
  for l = 1:length(t_phys)
    t_loc = t_phys(l);
    p_loc = pbar_prime_3d_anal*exp(i*omega*t_loc); p_loc = p_loc + conj(p_loc);
    z_loc = (-i*amp/2).*psi_in(bem_range).*exp(i*omega*t_loc); z_loc = z_loc + conj(z_loc);

    subplot(211);
    plot(col_pts(bem_range,1) ./ La, p_loc(bem_range) ./ flow.pref, 'b', 'LineWidth', 2); grid on; hold on;
    plot(linspace(0,1,41), pressure_data(dmd_range,l) - 1, 'r', 'LineWidth', 2);
    ylim([min_pf_cfd,max_pf_cfd]);
    ax = gca; ax.FontSize = afs; 
    legend('BEM', 'CFD', 'interpreter', 'latex', 'fontsize', fs, 'location', 'eastoutside');  
    ylabel('$p^{\prime} / p_{\infty}$', 'interpreter', 'latex', 'fontsize', fs);

    subplot(212);
    plot(col_pts(bem_range,1) ./ La, z_loc, 'b', 'LineWidth', 2); grid on; hold on;
    plot(linspace(0,1,41), z_data(dmd_range,l), 'r', 'LineWidth', 2);
    ylim([-amp, amp]);
    ax = gca; ax.FontSize = afs;
    ylabel('$z$', 'interpreter', 'latex', 'fontsize', fs);
    xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
    sgtitle(strcat(['Half-Span Comparison of BEM and CFD, $M = ', num2str(M), ', t = ' num2str(l), '/', num2str(maxT), '$']), 'interpreter', 'latex', 'fontsize', fs);
   
    frame = getframe(gcf);
    writeVideo(v_comp_3danal_halfspan, frame);
    clf;
  
  end
  close(v_comp_3danal_halfspan);

end

if (run_3dstanal == 1)
  f_comp_stanal = figure; f_comp_stanal.Position = [189 152 1263 630];
  v_fname = strcat(['./video_results/3d/m', m_choice, '_stanal_comparison.mp4']);
  v_comp_stanal = VideoWriter(v_fname, 'MPEG-4');
  v_comp_stanal.FrameRate = 20; v_comp_stanal.Quality = 10;
  open(v_comp_stanal);
  for l = 1:length(t_phys)
    t_loc = t_phys(l);
    p_loc = pbar_st*exp(i*omega*t_loc); p_loc = p_loc + conj(p_loc);
    z_loc = imag( psi_in.*(amp*exp(i*omega*t_loc)) );
  
    subplot(141);
    scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, p_loc ./ flow.pref, 'filled'); grid on;
    colorbar; colormap('jet'); caxis([min_pf_bem_st, max_pf_bem_st]); pbaspect([1 2 1]);
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
  
    sgtitle(strcat(['BEM vs. CFD, Strip Theory Analytical Geom., $M = ', num2str(M), ', t = ', num2str(l), '/', num2str(maxT), '$']), 'interpreter', 'latex', 'fontsize', fs);
    
    frame = getframe(gcf);
    writeVideo(v_comp_stanal, frame);
    clf;
  
  end
  close(v_comp_stanal);
end


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

