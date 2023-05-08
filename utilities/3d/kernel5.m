function [int_loc] = kernel5(x,y,bndry_pts,lambda,mu,flow,amp,omega,Uinf)

  %pts = [-1/sqrt(3), 1/sqrt(3)];
  pts = [-sqrt(3/5), 0, sqrt(3/5)];
  weights = [5/9, 8/9, 5/9];

  % boundary element corner values
  a = bndry_pts(1,1); b = bndry_pts(2,1);
  c = bndry_pts(1,2); d = bndry_pts(3,2);

  cent = [(a+b)/2, (c+d)/2];  % boundary element center coordinates

  il = [2/3, 1/6, 1/6; 1/6, 2/3, 1/6; 1/6, 1/6, 2/3];
  La = flow.La; Lb = flow.Lb;

  int_loc = 0;

  yupper = y + (x-b);
  ylower = y - (x-b);
  xupper = x - (d-y);
  xlower = x - (y-c);

  %figure;
  %scatter(bndry_pts(:,1), bndry_pts(:,2), 'b'); hold on;
  %scatter(x,y,'m','filled');
  %scatter(xlower,c,'g'); scatter(xupper,d,'g'); scatter(b,ylower,'r'); scatter(b,yupper,'r');

  % Compute left rectangle
  for l = 1:3
    for m = 1:3
      x1 = ((xupper-a)/2)*pts(l) + ((xupper+a)/2);
      y1 = ((d-c)/2)*pts(m) + ((c+d)/2);
      weight = weights(l)*weights(m);
      %scatter(x1,y1,'kx');

      R = sqrt( (x-x1)^2 - (y-y1)^2 );
      kernel_loc = cos( lambda*R/cos(mu) ) / R;
      dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(x1*La/tan(mu),y1*La,La,Lb) + (1/2)*amp*omega*flow.psi(x1*La/tan(mu),y1*La,La,Lb);
      dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* x1 ./ Uinf) .* dphi_dn;
      int_loc = int_loc + dphi_dn*kernel_loc*weight;
    end
  end
  int_loc = ((xupper-a)/2)*((d-c)/2)*int_loc;

  % Compute right rectangle
  for l = 1:3
    for m = 1:3
      x1 = ((b-xupper)/2)*pts(l) + ((xupper+b)/2);
      y1 = ((yupper-ylower)/2)*pts(m) + ((yupper+ylower)/2);
      weight = weights(l)*weights(m);
      %scatter(x1,y1,'kx');

      R = sqrt( (x-x1)^2 - (y-y1)^2 );
      kernel_loc = cos( lambda*R/cos(mu) ) / R;
      dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(x1*La/tan(mu),y1*La,La,Lb) + (1/2)*amp*omega*flow.psi(x1*La/tan(mu),y1*La,La,Lb);
      dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* x1 ./ Uinf) .* dphi_dn;
      int_loc = int_loc + ((b-xupper)/2)*((yupper-ylower)/2)*weight*dphi_dn*kernel_loc;
    end
  end

  % Compute upper triangle
  base = (b-xupper); height = (d - yupper); area = base*height/2;

  xl = xupper + (1/3)*base; yl = yupper + (1/3)*height; %scatter(xl,yl,'kx');
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (-27/48)*area*dphi_dn*kernel_loc;

  xl = xupper + (3/5)*base; yl = yupper + (1/5)*height; %scatter(xl,yl,'kx');
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

  xl = xupper + (1/5)*base; yl = yupper + (1/5)*height; %scatter(xl,yl,'kx');
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

  xl = xupper + (1/5)*base; yl = yupper + (3/5)*height; %scatter(xl,yl,'kx');
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;



  % 3pt Gaussian quadrature for triangle
  %xl = xupper + (1/6)*(b-xupper); yl = yupper + (1/6)*(d - yupper); %scatter(xl,yl,'kx');
  %R = sqrt( (x-xl)^2 - (y-yl)^2 );
  %dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  %dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  %kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

  %xl = xupper + (1/6)*(b-xupper); yl = yupper + (2/3)*(d - yupper); %scatter(xl,yl,'kx');
  %R = sqrt( (x-xl)^2 - (y-yl)^2 );
  %dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  %dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  %kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

  %xl = xupper + (2/3)*(b-xupper); yl = yupper + (1/6)*(d - yupper); %scatter(xl,yl,'kx');
  %R = sqrt( (x-xl)^2 - (y-yl)^2 );
  %dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  %dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  %kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

  % Compute lower triangle
  base = (b-xlower); height = (ylower - c); area = base*height/2;

  xl = xlower + (1/3)*base; yl = ylower - (1/3)*height; %scatter(xl,yl,'kx'); pause;
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (-27/48)*area*dphi_dn*kernel_loc;

  xl = xlower + (3/5)*base; yl = ylower - (1/5)*height; %scatter(xl,yl,'kx'); pause;
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

  xl = xlower + (1/5)*base; yl = ylower - (1/5)*height; %scatter(xl,yl,'kx'); pause;
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

  xl = xlower + (1/5)*base; yl = ylower - (3/5)*height; %scatter(xl,yl,'kx'); pause;
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;



  % 3pt Gaussian quadrature for triangles
  %xl = xlower + (1/6)*(b-xlower); yl = ylower - (1/6)*(ylower - c); %scatter(xl,yl,'kx');
  %R = sqrt( (x-xl)^2 - (y-yl)^2 );
  %dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  %dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  %kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

  %xl = xlower + (1/6)*(b-xlower); yl = ylower - (2/3)*(ylower - c); %scatter(xl,yl,'kx');
  %R = sqrt( (x-xl)^2 - (y-yl)^2 );
  %dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  %dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  %kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

  %xl = xlower + (2/3)*(b-xlower); yl = ylower - (1/6)*(ylower - c); %scatter(xl,yl,'kx'); pause;
  %R = sqrt( (x-xl)^2 - (y-yl)^2 );
  %dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  %dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  %kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

end
