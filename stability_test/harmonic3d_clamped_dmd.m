
close all;
clear all;
clc;

fs = 24;
afs = 20;
sz = 24;

% Add utilities
addpath('../utilities/3d/');
addpath('../utilities/stability/');
addpath('~/research/dmdc_panel/hashimoto_fem/struct_energies_utils');
addpath('~/research/dmdc_panel/hashimoto_fem/fluid_function_utils');

% Define beam attributes
Lx_phys = 0.2286; 
Ly_phys = 2*Lx_phys;
h = 0.00101854; 
Lx = Lx_phys;
Ly = Ly_phys;

E = 38610640824; % Pa
nu = 0.35;
rho_m = 1762.035; % kg

% Define constant flow attributes
rho = 0.497743562433333;
gamma = 1.4; R = 287;
flow.rhoref = rho;
flow.gamma = gamma;
flow.La = Lx_phys;

% Define constant structural properties
num_modes_x = 2; num_modes_y = 1; num_modes = num_modes_x * num_modes_y;
T1 = 0.031682495535872; U1 = 15308.97392989149;
T2 = 0.037153300330665; U2 = 122072.1959817805;
T2_span = 0.033908307236766; U2_span = 27458.33043857992;

T_vals = [T1, T2]; U_vals = [U1, U2];
%T_vals = [T1, T2_span]; U_vals = [U1, U2_span];

% Define trial modes
num_modes_x = 2; num_modes_y = 1; num_modes = num_modes_x * num_modes_y;
nmx = num2str(num_modes_x); nmy = num2str(num_modes_y);
Nx = 41; Ny = 41; N = (Nx-1)*(Ny-1);
x_dmd = linspace(0,Lx,Nx); y_dmd = linspace(0,Ly,Ny);
bnd_x_dmd = linspace(0,Lx,Nx); bnd_y_dmd = linspace(0,Ly,Ny);

Nx_bem = 31; Ny_bem = 62; N_bem = (Nx_bem-1)*(Ny_bem-1);
x = linspace(0,Lx,Nx); y = linspace(0,Ly,Ny);
col_pts = zeros(N_bem, 2);  bndry_pts = zeros(4*N_bem, 2);
bnd_x = linspace(0,Lx,Nx_bem); bnd_y = linspace(0,Ly,Ny_bem);
col_x = (bnd_x(1:end-1) + bnd_x(2:end)) ./ 2;
col_y = (bnd_y(1:end-1) + bnd_y(2:end)) ./ 2;
for k = 1:Ny_bem-1
  for j = 1:Nx_bem-1
    m = (k-1)*(Nx_bem-1) + j;
    col_pts(m,:) = [col_x(j), col_y(k)];
    bndry_pts(4*(m-1)+1, :) = [bnd_x(j), bnd_y(k)];
    bndry_pts(4*(m-1)+2, :) = [bnd_x(j+1), bnd_y(k)];
    bndry_pts(4*(m-1)+3, :) = [bnd_x(j), bnd_y(k+1)];
    bndry_pts(4*(m-1)+4, :) = [bnd_x(j+1), bnd_y(k+1)];
  end
end

% Read in FEM information -- for derivative computation
elem_node = importdata('./3d_dmd_modes_clamped/fem/elements.txt'); elem_node = elem_node.data;
node_elem = dlmread('./3d_dmd_modes_clamped/fem/node_elem.txt');
nodes = dlmread('./3d_dmd_modes_clamped/fem/nodes.txt');
structural_node_locs = nodes(:,2:4);
surface_node_data = dlmread('./3d_dmd_modes_clamped/fem/face_nodes.txt');
struct_faces = dlmread('./3d_dmd_modes_clamped/fem/surf_faces_trim.txt');
mode1_zdir_deform = dlmread('./3d_dmd_modes_clamped/fem/mode1_zdir_deform.txt');
mode2_zdir_deform = dlmread('./3d_dmd_modes_clamped/fem/mode2_zdir_deform.txt');
mode3_zdir_deform = dlmread('./3d_dmd_modes_clamped/fem/mode2_span_zdir_deform.txt');
def_func1 = mode1_zdir_deform(:,5); def_func1 = def_func1 ./ max(def_func1);
def_func2 = mode2_zdir_deform(:,5); def_func2 = def_func2 ./ max(def_func2);
def_func3 = mode3_zdir_deform(:,5); def_func3 = def_func3 ./ max(def_func3);

N_nodes = max(size(structural_node_locs));
shift = [zeros(N_nodes, 2), 0.00101854.*ones(N_nodes,1)];
structural_node_locs = structural_node_locs - shift;
N_face = max(size(surface_node_data));
shift = [zeros(N_face, 2), 0.00101854.*ones(N_face,1)];
surface_node_data(:,2:4) = surface_node_data(:,2:4) - shift;

modeshape_locs = zeros((Nx-1)*(Ny-1),2);

for k = 1:Ny_bem-1
  for j = 1:Nx_bem-1
    ind_col = (k-1)*(Nx_bem-1) + j;
    loc = [col_x(j), col_y(k)];
    col_locs(ind_col,:) = loc;
    [psi1_col_loc, dpsi1_dx_col_loc] = compute_mode_deriv(elem_node, node_elem, structural_node_locs, surface_node_data, struct_faces, def_func1, loc);
    [psi2_col_loc, dpsi2_dx_col_loc] = compute_mode_deriv(elem_node, node_elem, structural_node_locs, surface_node_data, struct_faces, def_func2, loc);

    mode_col1(ind_col) = psi1_col_loc; mode_deriv_col1(ind_col) = dpsi1_dx_col_loc;
    mode_col2(ind_col) = psi2_col_loc; mode_deriv_col2(ind_col) = dpsi2_dx_col_loc;

    %mode_deriv_col1(ind_col) = compute_mode_deriv(elem_node, node_elem, structural_node_locs, surface_node_data, struct_faces, def_func1, loc);
    %mode_deriv_col2(ind_col) = compute_mode_deriv(elem_node, node_elem, structural_node_locs, surface_node_data, struct_faces, def_func2, loc);
  end
end

figure;
subplot(121);
scatter(col_locs(:,1), col_locs(:,2), sz, mode_col1); hold on;
colorbar; colormap('jet');
title('$\psi_1^{BEM}$', 'interpreter', 'latex');
subplot(122); scatter(col_locs(:,1), col_locs(:,2), sz, mode_col2); hold on;
colorbar; colormap('jet');
title('$\psi_2^{BEM}$', 'interpreter', 'latex');

figure;
subplot(121);
scatter(col_locs(:,1), col_locs(:,2), sz, mode_deriv_col1); hold on;
colorbar; colormap('jet');
title('$\partial \psi_1^{BEM} / \partial x$', 'interpreter', 'latex');
subplot(122); scatter(col_locs(:,1), col_locs(:,2), sz, mode_deriv_col2); hold on;
colorbar; colormap('jet');
title('$\partial \psi_2^{BEM} / \partial x$', 'interpreter', 'latex'); %pause;

for k = 1:Ny
  for j = 1:Nx
    ind = (k-1)*Nx + j;
    loc = [bnd_x_dmd(j), bnd_y_dmd(k)];
    modeshape_locs(ind,:) = [bnd_x_dmd(j), bnd_y_dmd(k)];
    if ((j == 1) || (j == Nx) || (k == 1) || (k == Ny))
      mode_deriv1(ind) = 0; mode_deriv2(ind) = 0;
    else
      [psi1, dpsi1_dx] = compute_mode_deriv(elem_node, node_elem, structural_node_locs, surface_node_data, struct_faces, def_func1, loc);
      [psi2, dpsi2_dx] = compute_mode_deriv(elem_node, node_elem, structural_node_locs, surface_node_data, struct_faces, def_func2, loc);
      [psi3, dpsi3_dx] = compute_mode_deriv(elem_node, node_elem, structural_node_locs, surface_node_data, struct_faces, def_func3, loc);

      mode_deriv1(ind) = dpsi1_dx;
      mode_deriv2(ind) = dpsi2_dx;
      %mode_deriv2(ind) = dpsi3_dx;

      %mode_deriv1(ind) = compute_mode_deriv(elem_node, node_elem, structural_node_locs, surface_node_data, struct_faces, def_func1, loc);
      %mode_deriv2(ind) = compute_mode_deriv(elem_node, node_elem, structural_node_locs, surface_node_data, struct_faces, def_func2, loc);
    end
  end
end
mode_deriv(1,:) = mode_deriv1; mode_deriv(2,:) = mode_deriv2;

for k = 1:num_modes_y
  for j = 1:num_modes_x

    % Construct mode shapes
    ind = (k-1)*num_modes_x + j;
    %psi = modeshape(j,col_x,Lx)' * modeshape(k,col_y,Ly);
    %psi = reshape(psi,1,N); bem_psi(ind,:) = psi;
    %dpsi_dx = mode_deriv(j,col_x,Lx)' * modeshape(k,col_y,Ly);
    %dpsi_dx = reshape(dpsi_dx,1,N); bem_dpsi_dx(ind,:) = dpsi_dx;

    if (j==1)
      mode_loc = dlmread('./3d_dmd_modes_clamped/Nx41_Ny41_panel_deform.csv');
    else
      mode_loc = dlmread('./3d_dmd_modes_clamped/Nx41_Ny41_panel_deform_mode5.csv');
      %mode_loc = dlmread('./3d_dmd_modes_clamped/Nx41_Ny41_panel_deform_mode2_span.csv');
    end
    modeshapes(ind,:) = reshape(mode_loc', 1, Nx*Ny);

    % Construct stiffness matrix
    T(ind,ind) = T_vals(j);
    U(ind,ind) = U_vals(j);

  end
end

figure;
subplot(121);
scatter(modeshape_locs(:,1), modeshape_locs(:,2), sz, modeshapes(1,:)); grid on; hold on;
colorbar; colormap('jet'); ax = gca; ax.FontSize = afs;
ylabel('$y$', 'interpreter', 'latex', 'fontsize', fs);
xlabel('$x$', 'interpreter', 'latex', 'fontsize', fs);
title('$\psi_1$', 'interpreter', 'latex', 'fontsize', fs);
subplot(122);
scatter(modeshape_locs(:,1), modeshape_locs(:,2), sz, modeshapes(2,:)); grid on; hold on;
colorbar; colormap('jet'); ax = gca; ax.FontSize = afs;
xlabel('$$', 'interpreter', 'latex', 'fontsize', fs);
title('$\psi_2$', 'interpreter', 'latex', 'fontsize', fs);

figure;
subplot(121);
scatter(modeshape_locs(:,1), modeshape_locs(:,2), sz, mode_deriv(1,:)); grid on; hold on;
colorbar; colormap('jet'); ax = gca; ax.FontSize = afs;
ylabel('$y / h$', 'interpreter', 'latex', 'fontsize', fs);
xlabel('$x / h$', 'interpreter', 'latex', 'fontsize', fs);
title('$\partial \psi_1 / \partial x$', 'interpreter', 'latex', 'fontsize', fs);
subplot(122);
scatter(modeshape_locs(:,1), modeshape_locs(:,2), sz, mode_deriv(2,:)); grid on; hold on;
colorbar; colormap('jet'); ax = gca; ax.FontSize = afs;
xlabel('$x / h$', 'interpreter', 'latex', 'fontsize', fs);
title('$\partial \psi_2 / \partial x$', 'interpreter', 'latex', 'fontsize', fs);

% Compute natural frequencies
for j = 1:num_modes
  omega_nat(j) = sqrt( U_vals(j) / T_vals(j) );
end

%% Begin stability search
%M_min = 1.1; M_max = 1.5; Mn = 4;
M_min = 1.2; M_max = 1.5; Mn = 4;
%M_list = [1.4];
M_list = linspace(M_min, M_max, Mn);
omega_soln = zeros(Mn, num_modes);
p_bem_orig = zeros(num_modes, (Nx-1)*(Ny-1), Mn);
p_bem_soln = zeros(num_modes, (Nx-1)*(Ny-1), Mn);
p_bem_loc = zeros((Nx-1)*(Ny-1), num_modes);

disp('Beginning stability computations, pause.'); %pause;
for m = 1:Mn
  M = M_list(m);
  m_choice = num2str(M); m_choice = strcat([m_choice(1), m_choice(3:end)]);
  flow.M = M;

  if (M == 1.1)
    pref = 41696.99607460994;
    Tref = 291.8886582844682;
  elseif (M == 1.2)
    pref = 35037.05920158196;
    Tref = 245.2675531418101;
  elseif (M == 1.3)
    pref = 29854.06227827102;
    Tref = 208.9853707243825;
  elseif (M == 1.4)
    pref = 25741.51288279492;
    Tref = 180.1965696552075;
  else
    pref = 22423.7178890121911;
    Tref = 156.9712340107566;
  end
  flow.pref = pref; flow.Tref = Tref;

  disp(strcat(['Solving flow regime M = ', num2str(M)]));

  for j = 1:num_modes

    % Pull natural frequency
    omega0 = omega_nat(j);
  
    % Compute pressure matrix with current omega
    P = zeros(num_modes, num_modes);
    for l = 1:num_modes
      %[phi_bem, p_bem] = bem3d(h.*col_pts, h.*bndry_pts, Nx-1, h.*bem_psi(:,l), bem_dpsi_dx(:,l), 2*i, -(a/h)*omega0, flow);
      %[phi_bem, p_bem] = bem3d(col_pts, bndry_pts, Nx_bem-1, mode_col1, mode_deriv_col1, 1, omega0, flow);
      %p_bem = (1/(rho_m*a^2)).*p_bem; if (l==j) p_bem_loc(:,j) = p_bem; end

      aref = sqrt(gamma*R*Tref);
      p_lpt_nondim = p_lpt_func(M, gamma, omega0/aref, modeshapes(l,:), mode_deriv(l,:)); p_lpt = pref.*p_lpt_nondim;
      dmd_mode = dlmread(strcat(['./3d_dmd_modes_clamped/m', m_choice, '/mode', num2str(l), '_omega', num2str(j), '/dmd_modes.csv']));
      alpha_tilde = dlmread(strcat(['./3d_dmd_modes_clamped/m', m_choice, '/mode', num2str(l), '_omega', num2str(j), '/alphas_tilde.csv']));
      dmd_choice = 1;
      dmd_mode = dmd_mode(:,dmd_choice); alpha_tilde = alpha_tilde(dmd_choice);
      p_corr_nondim = p_lpt_nondim.' + (2*i).*dmd_mode.*alpha_tilde;
      p_corr = pref .* p_corr_nondim;
      %p_corr = p_lpt.' + pref.*dmd_mode.*alpha_tilde; p_corr_nondim = p_corr ./ pref;

      %figure;
      %subplot(121);
      %scatter(col_pts(:,1), col_pts(:,2), sz, real(p_bem));
      %ax = gca; ax.FontSize = afs; colorbar; colormap('jet');
      %xlabel('$x$', 'interpreter', 'latex', 'fontsize', fs);
      %ylabel('$y$', 'interpreter', 'latex', 'fontsize', fs);
      %title('Re( $p^{\prime}_{BEM}$ )', 'interpreter', 'latex', 'fontsize', fs);
      %subplot(122);
      %scatter(col_pts(:,1), col_pts(:,2), sz, imag(p_bem));
      %ax = gca; ax.FontSize = afs; colorbar; colormap('jet');
      %xlabel('$x$', 'interpreter', 'latex', 'fontsize', fs);
      %ylabel('$y$', 'interpreter', 'latex', 'fontsize', fs);
      %title('Im( $p^{\prime}_{BEM}$ )', 'interpreter', 'latex', 'fontsize', fs);

      %figure;
      %subplot(121);
      %scatter(modeshape_locs(:,1), modeshape_locs(:,2), sz, real(dmd_mode));
      %ax = gca; ax.FontSize = afs; colorbar; colormap('jet');
      %xlabel('$x$', 'interpreter', 'latex', 'fontsize', fs);
      %ylabel('$y$', 'interpreter', 'latex', 'fontsize', fs);
      %title('Re( $\phi$ )', 'interpreter', 'latex', 'fontsize', fs);
      %subplot(122);
      %scatter(modeshape_locs(:,1), modeshape_locs(:,2), sz, imag(dmd_mode));
      %ax = gca; ax.FontSize = afs; colorbar; colormap('jet');
      %xlabel('$x$', 'interpreter', 'latex', 'fontsize', fs);
      %ylabel('$y$', 'interpreter', 'latex', 'fontsize', fs);
      %title('Im( $\phi$ )', 'interpreter', 'latex', 'fontsize', fs);

      figure;
      subplot(121);
      scatter(modeshape_locs(:,1), modeshape_locs(:,2), sz, real(p_lpt_nondim));
      ax = gca; ax.FontSize = afs; colorbar; colormap('jet');
      xlabel('$x$', 'interpreter', 'latex', 'fontsize', fs);
      ylabel('$y$', 'interpreter', 'latex', 'fontsize', fs);
      title('Re( $p^{\prime}_{LPT}$ )', 'interpreter', 'latex', 'fontsize', fs);
      subplot(122);
      scatter(modeshape_locs(:,1), modeshape_locs(:,2), sz, imag(p_lpt_nondim));
      ax = gca; ax.FontSize = afs; colorbar; colormap('jet');
      xlabel('$x$', 'interpreter', 'latex', 'fontsize', fs);
      ylabel('$y$', 'interpreter', 'latex', 'fontsize', fs);
      title('Im( $p^{\prime}_{LPT}$ )', 'interpreter', 'latex', 'fontsize', fs);

      figure;
      subplot(121);
      scatter(modeshape_locs(:,1), modeshape_locs(:,2), sz, real(p_corr_nondim));
      ax = gca; ax.FontSize = afs; colorbar; colormap('jet');
      xlabel('$x$', 'interpreter', 'latex', 'fontsize', fs);
      ylabel('$y$', 'interpreter', 'latex', 'fontsize', fs);
      title('Re( $p^{\prime}_{DMD}$ )', 'interpreter', 'latex', 'fontsize', fs);
      subplot(122);
      scatter(modeshape_locs(:,1), modeshape_locs(:,2), sz, imag(p_corr_nondim));
      ax = gca; ax.FontSize = afs; colorbar; colormap('jet');
      xlabel('$x$', 'interpreter', 'latex', 'fontsize', fs);
      ylabel('$y$', 'interpreter', 'latex', 'fontsize', fs);
      title('Im( $p^{\prime}_{DMD}$ )', 'interpreter', 'latex', 'fontsize', fs);
      sgtitle(strcat(['DMD Solution, $M = ' num2str(M), '$']), 'interpreter', 'latex', 'fontsize', fs);
      %pause;

      for k = 1:num_modes
        pjk = gauss_1pt_quadrature_3d(modeshape_locs, p_corr, modeshapes(k,:), Nx, Ny);
        P(l,k) = pjk;
      end

    end

    % Compute stability matrix
    A = -T.*(omega0^2) + U + P;
    A = sym(A);
    syms wpp
    A(j,j) = -T(j,j)*wpp^2 + U(j,j) + P(j,j);
    eqn = det(A) == 0;
    B = double(solve(eqn));

    omegap(j) = B(real(B)>0);

  end

  disp('Frequency solutions:');
  disp(omegap);
  omega_soln(m,:) = omegap;

end

figure;
for j = 1:num_modes_x
  sp = str2num(['2', nmx, num2str(j)]);
  subplot(sp);
  plot(M_list, real(omega_soln(:,j)), '-or'); grid on;
  ax = gca; ax.FontSize = afs;
  if (j==1) ylabel('Re( $\omega$ )', 'interpreter', 'latex', 'fontsize', fs); end
  title(strcat(['$\omega_', num2str(j), '$']), 'interpreter', 'latex', 'fontsize', fs);

  sp = str2num(['2', nmx, num2str(j+num_modes_x)]);
  subplot(sp);
  plot(M_list, imag(omega_soln(:,j)), '-or'); grid on;
  ax = gca; ax.FontSize = afs;
  if (j==1) ylabel('Im( $\omega$ )', 'interpreter', 'latex', 'fontsize', fs); end
  xlabel('$M$', 'interpreter', 'latex', 'fontsize', fs);
end
sgtitle(strcat(['Eigenfrequencies, $N_x = ', num2str(Nx), ', N_y = ', num2str(Ny), ', L_x = ', num2str(Lx), ', L_y = ', num2str(Ly), '$']), 'interpreter', 'latex', 'fontsize', fs);

figure;
subplot(121);
plot(M_list, real(omega_soln(:,1)), '-or', 'LineWidth', 2); grid on;
ax = gca; ax.FontSize = afs;
ylabel('Re( $\omega_1$ )', 'interpreter', 'latex', 'fontsize', fs);
subplot(122);
plot(M_list, imag(omega_soln(:,1)), '-or', 'LineWidth', 2); grid on;
ax = gca; ax.FontSize = afs;
ylabel('Im( $\omega_1$ )', 'interpreter', 'latex', 'fontsize', fs);
xlabel('$M$', 'interpreter', 'latex', 'fontsize', fs);
sgtitle(strcat(['Eigenfrequency Solutions, $N_x = ', num2str(Nx), ', N_y = ', num2str(Ny), '$']), 'interpreter', 'latex', 'fontsize', fs);

function [p_lpt_nondim] = p_lpt_func(M, gamma, k, modeshape, mode_deriv)

  % exp( i * omega * t ) convection
  p_lpt_nondim = gamma*M.*mode_deriv + i*k.*modeshape;
  %p_lpt_nondim = 0.5.*(i*gamma*M.*mode_deriv + gamma*k.*modeshape);

end

function [pjk] = gauss_1pt_quadrature_3d(modeshape_locs, p, psi, Nx, Ny)

  dx = modeshape_locs(2,1) - modeshape_locs(1,1);
  dy = modeshape_locs(Nx+1,2) - modeshape_locs(1,2);

  pjk = 0;
  for k = 1:Ny
    for j = 1:Nx
      ind = (k-1)*Nx + j;
      if ( ((j==1)&&(k==1)) || ((j==Nx)&&(k==1)) || ((j==1)&&(k==Ny)) || ((j==Nx)&&(k==Ny)) )
        % Corner cases
        pjk = pjk + (dx/2)*(dy/2)*p(ind)*psi(ind);
      elseif ( (j==1) || (j==Nx) )
        pjk = pjk + dx*(dy/2)*p(ind)*psi(ind);
      elseif ( (k==1) || (k==Ny) )
        pjk = pjk + dy*(dx/2)*p(ind)*psi(ind);
      else
        pjk = pjk + dx*dy*p(ind)*psi(ind);
      end
    end
  end
  %pjk = -1*pjk;

end

function [psi, dpsi_dx] = compute_mode_deriv(elem_node, node_elem, structural_node_locs, surface_node_data, struct_faces, def_func, loc)

  neg_epsilon = -1e-10;

  loc = [loc, 0];
  [belong_elem, fluid_node_bary] = find_struct_face(loc, surface_node_data, structural_node_locs, node_elem, elem_node, struct_faces, neg_epsilon);
  node_bary = [fluid_node_bary, belong_elem];
  belong_elem_nodes = elem_node(belong_elem, :);
  belong_elem_deform = def_func(belong_elem_nodes);
  loc_deform = bary_tet_interp(fluid_node_bary, belong_elem_deform);
  psi = loc_deform;
  
  % Worth a double check to be safe, if results look poor
  belong_elem_nodelocs = structural_node_locs(belong_elem_nodes, :);
  [detJ, J, Jinv, L] = jacobian(fluid_node_bary(2), fluid_node_bary(3), fluid_node_bary(4), belong_elem_nodelocs);
  q = reshape([zeros(10,2), belong_elem_deform]', [], 1);
  [eps_xx, eps_yy, eps_zz, eps_xy, eps_yz, eps_xz] = strain(L, Jinv, q);
  dpsi_dx = 2*eps_xz;

end

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


