function [G] = greens3d_subsonic(x,y,z,xi,eta,zeta,k,M);

  % Compute half-space subsonic Greens function

  [GFS] = greensFS3d_subsonic(x,y,z,xi,eta,zeta,k,M);
  [GFSim] = greensFS3d_subsonic(x,y,z,xi,eta,-zeta,k,M);

  G = GFS + GFSim;

end

