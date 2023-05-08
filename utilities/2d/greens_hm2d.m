function [G, dGdx, dGdy] = greens_hm2d(x,y,xi,eta,k)

  % Compute half-space Greens function to Helmholtz equation
  [GFS, dGFSdx, dGFSdy] = greensFS_hm2d(x,y,xi,eta,k);
  [GFSim, dGFSdxim, dGFSdyim] = greensFS_hm2d(x,y,xi,-eta,k);

  G = GFS + GFSim;
  dGdx = dGFSdx + dGFSdxim;
  dGdy = dGFSdy + dGFSdyim;

end

