function [int] = gauss_1pt_quad_3d(bndry_pts, p_nondim, bem_psi)

  N = length(bndry_pts) / 4;
  int = 0;
  for j = 1:N

    a = bndry_pts(4*(j-1)+1,1); b = bndry_pts(4*(j-1)+2,1);
    c = bndry_pts(4*(j-1)+1,2); d = bndry_pts(4*(j-1)+3,2);
    loc_int = (b-a)*(d-c)*p_nondim(j)*bem_psi(j);
    int = int + loc_int;

  end

end
