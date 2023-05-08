function [G] = greens3d_supersonic(x,y,z,xi,eta,zeta,k,M,La);

  % Compute half-space supersonic Greens function

  [GFS] = greensFS3d_supersonic(x,y,z,xi,eta,zeta,k,M,La);
  [GFSim] = greensFS3d_supersonic(x,y,z,xi,eta,-zeta,k,M,La);

  G = GFS + GFSim;

end

