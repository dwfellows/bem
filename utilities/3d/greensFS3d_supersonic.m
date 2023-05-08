function [GFS] = greensFS3d_supersonic(x,y,z,xi,eta,zeta,k,M,La)

  % Compute subsonic Green's function, derivatives, and radius info from points
  beta = sqrt(M^2 - 1);

  if ( (x-xi) > beta*sqrt( (y-eta)^2 + (z-zeta)^2 ) )
    H = 1;
  else
    H = 0;
  end

  R = (1/La)*sqrt( (1/beta^2)*(x-xi)^2 + (y-eta)^2 + (z-zeta)^2 );

  GFS = (H/R) * cos( (k*La/beta) * R );

end

