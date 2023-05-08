function [GFS, dGFSdx, dGFSdy, dRdx, dRdy] = greensFS2d(x,y,xi,eta,k,M)

  % Compute Green's function, derivatives, and radius info from points
  beta = sqrt(1 - M^2);

  % Compute full-space information
  R = sqrt( ((x-xi)^2)/(beta^2) + (y-eta)^2 );
  dRdx = (x-xi) / ( (beta^2) * sqrt( ( ((x-xi)^2)/(beta^2) ) + (y-eta)^2 ) );
  dRdy = (y-eta) / sqrt( ( ((x-xi)^2)/(beta^2) ) + (y-eta)^2 );

  GFS = (-i/(4*beta)) * besselh(0,2, k*M*R/beta) * exp( -(i*k*M*(x-xi))/(beta^2) );

  dGFSdx = (-i/(4*beta)) * (-besselh(1,2, k*M*R/beta) * (k*M/beta) * dRdx) * exp( -(i*k*M*(x-xi))/(beta^2) ) ...
           -i/(4*beta) * besselh(0,2, k*M*R/beta) * (-i*k*M/(beta^2)) * exp( -(i*k*M*(x-xi))/(beta^2) );

  dGFSdy = (-i/(4*beta)) * (-besselh(1,2, k*M*R/beta) * (k*M/beta) * dRdy) * exp( -(i*k*M*(x-xi))/(beta^2) );

  GFS = (-i/(4*beta)) * besselh(0,2, (k/(beta^2)) * sqrt( (x-xi)^2 + (beta^2)*(y-eta)^2 ) ) * exp( (-i*k*M*(x-xi))/(beta^2) );
  GFS = (-i/(4*beta)) * besselh(0,2, (k/(beta^2)) * sqrt( (x-xi)^2 + (beta^2)*(y-eta)^2 ) ) * exp( (i*k*M*(x-xi))/(beta^2) );

end

