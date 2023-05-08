function [GFS] = greensFS3d_subsonic(x,y,z,xi,eta,zeta,k,M)

  % Compute subsonic Green's function, derivatives, and radius info from points
  beta = sqrt(1 - M^2);

  % Compute full-space information
  R = sqrt( ((x-xi)^2)/(beta^2) + (y-eta)^2 + (z-zeta)^2 );
  %GFS = (1/(4*pi*beta*R)) * exp(-i*k*R/beta) * exp(-i*k*M*(x-xi)/(beta^2));

  %GFS = exp( -i*k*( (sqrt((x-xi)^2 + (beta^2)*((y-eta)^2 + (z-zeta)^2)) + M*(x-xi))/(beta^2)) ) ./ ...
        sqrt((x-xi)^2 + (beta^2)*((y-eta)^2 + (z-zeta)^2));

  GFS = exp( -i*k*( (sqrt((x-xi)^2 + (beta^2)*((y-eta)^2 + (z-zeta)^2)) + M*(xi-x))/(beta^2)) ) ./ ...
        sqrt((x-xi)^2 + (beta^2)*((y-eta)^2 + (z-zeta)^2));


end

