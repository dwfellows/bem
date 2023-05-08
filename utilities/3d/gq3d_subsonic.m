function [int_loc] = gq3d_subsonic(x,y,z,bndry_pts,k,M)

  pts = [-1/sqrt(3), 1/sqrt(3)];

  a = bndry_pts(1,1); b = bndry_pts(2,1);
  c = bndry_pts(1,2); d = bndry_pts(3,2);

  int_loc = 0;
  for l = 1:2
    for m = 1:2
      [G] = greens3d_subsonic(x,y,z, ((b-a)/2)*pts(l) + ((a+b)/2), ((d-c)/2)*pts(m) + ((c+d)/2), 0, k, M);
      int_loc = int_loc + G;
    end
  end
  int_loc = ((b-a)/2)*((d-c)/2)*int_loc;

end
