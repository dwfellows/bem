
close all;
clear all;
clc;

addpath('./utilities/3d/');

fs = 20;
afs = 18;
lfs = 20;
sz = 36;

% Define problem parameters
m_choice = '05';
M = str2num(strcat([m_choice(1), '.', m_choice(2:end)]));
amp = 0.0000508;
omega = 691.15038379;

flow.pref = 34601.89642944;
flow.rhoref = 0.49774215;
flow.gamma = 1.4;
flow.Tref = 242.222;
flow.M = M;


%% Define BEM collocation and boundary points
La = 0.2286; Lb = 2*La; flow.La = La;
Nx = 21; Ny = 41; N = (Nx-1)*(Ny-1);
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
neg_epsilon = -1e-10;
addpath('~/research/dmdc_panel/hashimoto_fem/struct_energies_utils');
addpath('~/research/dmdc_panel/hashimoto_fem/fluid_function_utils');

% Import structural quantities
elem_node = importdata('~/research/dmdc_panel/hashimoto_fem/elements_hires.txt'); elem_node = elem_node.data;
node_elem = dlmread('~/research/dmdc_panel/hashimoto_fem/node_elem_hires.txt');
nodes = dlmread('~/research/dmdc_panel/hashimoto_fem/nodes_hires.txt');
structural_node_locs = nodes(:,2:4);
surface_node_data = dlmread('~/research/dmdc_panel/hashimoto_fem/face_nodes_hires.txt');
struct_faces = dlmread('~/research/dmdc_panel/hashimoto_fem/surf_faces_trim_hires.txt');

mode1_zdir_deform = dlmread('~/research/dmdc_panel/hashimoto_fem/mode1_zdir_deform_hires.txt');
def_func = mode1_zdir_deform(:,5);
def_func = def_func ./ max(def_func);

% fix location of structural origin
N_nodes = max(size(structural_node_locs));
shift = [zeros(N_nodes, 2), 0.00101854.*ones(N_nodes,1)];
structural_node_locs = structural_node_locs - shift;
N_face = max(size(surface_node_data));
shift = [zeros(N_face, 2), 0.00101854.*ones(N_face,1)];
surface_node_data(:,2:4) = surface_node_data(:,2:4) - shift;

% Compute psi and dpsi_dx
psi = zeros(N,1);
for l = 1:N

  loc = [col_pts(l,:), 0];
  [belong_elem, fluid_node_bary] = find_struct_face(loc, surface_node_data, structural_node_locs, node_elem, elem_node, struct_faces, neg_epsilon);
  node_bary(l,:) = [fluid_node_bary, belong_elem];
  belong_elem_nodes = elem_node(belong_elem, :);
  belong_elem_deform = def_func(belong_elem_nodes);
  loc_deform = bary_tet_interp(fluid_node_bary, belong_elem_deform);
  psi(l) = loc_deform;

  % Worth a double check to be safe, if results look poor
  belong_elem_nodelocs = structural_node_locs(belong_elem_nodes, :);
  [detJ, J, Jinv, L] = jacobian(fluid_node_bary(2), fluid_node_bary(3), fluid_node_bary(4), belong_elem_nodelocs);
  q = reshape([zeros(10,2), belong_elem_deform]', [], 1);
  [eps_xx, eps_yy, eps_zz, eps_xy, eps_yz, eps_xz] = strain(L, Jinv, q);
  dpsi_dx(l) = 2*eps_xz;

end

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
scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, real(pbar_prime), 'filled');
ax = gca; ax.FontSize = afs;
colorbar; colormap('jet'); pbaspect([1 2 1]);
xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
title('$Re \big( \bar{p}^{\prime} \big)$', 'interpreter', 'latex', 'fontsize', fs);
sgtitle(strcat(['Reduced Velocity Perturbation Potential and Pressure Fluctuation, $M = ', num2str(M), '$']), 'interpreter', 'latex', 'fontsize', fs);

subplot(224);
scatter(col_pts(:,1) ./ La, col_pts(:,2) ./ La, sz, imag(pbar_prime), 'filled');
ax = gca; ax.FontSize = afs;
colorbar; colormap('jet'); pbaspect([1 2 1]);
xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
title('$Im \big( \bar{p}^{\prime} \big)$', 'interpreter', 'latex', 'fontsize', fs);
sgtitle(strcat(['Reduced Velocity Perturbation Potential and Pressure Fluctuation, $M = ', num2str(M), '$']), 'interpreter', 'latex', 'fontsize', fs);

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

pressure_data = dlmread(strcat(['./comparison_data/3d/panel_runs_subsonic/m', m_choice, '/surface_pressure_history_m', m_choice, '.csv']));
y_data = dlmread(strcat(['./comparison_data/3d/panel_runs_subsonic/m', m_choice, '/surface_yloc_history_m', m_choice, '.csv']));
x_data = dlmread(strcat(['./comparison_data/3d/panel_runs_subsonic/m', m_choice, '/surface_xloc_history_m', m_choice, '.csv']));
z_data = dlmread(strcat(['./comparison_data/3d/panel_runs_subsonic/m', m_choice, '/surface_zloc_history_m', m_choice, '.csv']));
pressure_data = pressure_data';
x_data = x_data';
y_data = y_data'; pd_y_locs = y_data(:,1); pd_y_locs = pd_y_locs - pd_y_locs(1);
z_data = z_data'; pd_x_locs = x_data(:,1); pd_x_locs = pd_x_locs - pd_x_locs(1);

max_pf_bem = max(max( imag( (pbar_prime') * exp(i*omega*t_phys) )./flow.pref ));
min_pf_bem = min(min( imag( (pbar_prime') * exp(i*omega*t_phys) )./flow.pref ));
max_pf_cfd = max(max( pressure_data - 1 ));
min_pf_cfd = min(min( pressure_data - 1 ));

f_comp = figure; f_comp.Position = [189 152 1263 630];
v_fname = strcat(['./video_results/3d/m', m_choice, '_comparison.mp4']);
v_comp = VideoWriter(v_fname, 'MPEG-4');
v_comp.FrameRate = 20; v_comp.Quality = 10;
open(v_comp);
for l = 1:length(t_phys)
  t_loc = t_phys(l);
  p_loc = imag( pbar_prime*exp(i*omega*t_loc) );
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

