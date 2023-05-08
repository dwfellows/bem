function [int_loc] = gq2d(x,y,bndry_pts,k,M)

  pts = [-1/sqrt(3), 1/sqrt(3)];
  a = bndry_pts(1); b = bndry_pts(2);
  int_loc = 0;

  if (M < 1) % If subsonic, utilize subsonic Green's function
    for j = 1:2
      [G,dGdx,dGdy] = greens2d(x,y, ((b-a)/2)*pts(j) + ((a+b)/2), 0, k, M);
      int_loc = int_loc + G;
    end
    int_loc = ((b-a)/2)*int_loc;
  else % If supersonic, utilize supersonic Green's function
    f = @(x1) besselj(0,k*(x - x1));
    if (b < x)
      %int_loc = integral(f,a,b);
      int_loc = 0;
      for j = 1:2
        int_loc = int_loc + f( ((b-a)/2)*pts(j) + ((a+b)/2) );
      end
      int_loc = ((b-a)/2)*int_loc;
    else
      %int_loc = integral(f,a,x);
      int_loc = f( (x-a)/2 );
      int_loc = (x-a)*int_loc;
    end
  end

end
