function [int] = gauss_1pt_quad(bndry_pts, p_nondim, bem_psi)

  N = length(bndry_pts) - 1;
  int = 0;
  for j = 1:N

    a = bndry_pts(j); b = bndry_pts(j+1);
    loc_int = ((b-a)/2)*2*p_nondim(j)*bem_psi(j);
    int = int + loc_int;

  end

end
