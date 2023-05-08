function [G, dGdx, dGdy] = greensFS_hm2d(x,y,xi,eta,k)

  % Compute full-space Greens function to Helmholtz equation
  r = sqrt( (x-xi)^2 + (y-eta)^2 );
  G = -(i/4) * besselh(0,1,k*r)

  drdx = (x-xi) / (2*sqrt( (x-xi)^2 + (y-eta)^2 ));
  drdy = (y-eta) / (2*sqrt( (x-xi)^2 + (y-eta)^2 ));

  dGdx = (i*k/4) * besselh(1,1,k*r) * drdx;
  dGdy = (i*k/4) * besselh(1,1,k*r) * drdy;

end

