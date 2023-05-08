function [col_pts, bndry_pts] = construct_grid(La, Lb, Nx, Ny)

  N = Nx*Ny;

  col_pts = zeros(N,2); bndry_pts = zeros(4*N, 2);
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

end
