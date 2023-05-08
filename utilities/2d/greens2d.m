function [G, dGdx, dGdy] = greens2d(x,y,xi,eta,k,M);

  % Compute half-space Greens function

  [GFS, dGFSdx, dGFSdy, dRdx, dRdy] = greensFS2d(x,y,xi,eta,k,M);
  [GFSim, dGFSdxim, dGFSdyim, dRdxim, dRdyim] = greensFS2d(x,y,xi,-eta,k,M);

  G = GFS + GFSim;
  dGdx = dGFSdx + dGFSdxim;
  dGdy = dGFSdy + dGFSdyim;

end

