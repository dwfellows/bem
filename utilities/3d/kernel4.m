function [int_loc] = kernel4(x,y,bndry_pts,lambda,mu,flow,amp,omega,Uinf)

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

  yupper_left = y + (x-a);
  ylower_left = y - (x-a);
  yupper_right = y + (x-b);
  ylower_right = y - (x-b);

  %figure;
  %scatter(bndry_pts(:,1), bndry_pts(:,2), 'b'); hold on;
  %scatter(x,y,'m','filled');
  %scatter(a,yupper_left, 'r'); scatter(a,ylower_left,'r');
  %scatter(b,yupper_right, 'g'); scatter(b,ylower_right,'g');

  % Compute inner rectangle
  for l = 1:3
    for m = 1:3
      x1 = ((b-a)/2)*pts(l) + ((a+b)/2);
      y1 = ((yupper_right-ylower_right)/2)*pts(m) + ((yupper_right+ylower_right)/2);
      weight = weights(l)*weights(m);
      %scatter(x1,y1,'kx');

      R = sqrt( (x-x1)^2 - (y-y1)^2 );
      kernel_loc = cos( lambda*R/cos(mu) ) / R;
      dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(x1*La/tan(mu),y1*La,La,Lb) + (1/2)*amp*omega*flow.psi(x1*La/tan(mu),y1*La,La,Lb);
      dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* x1 ./ Uinf) .* dphi_dn;

      int_loc = int_loc + dphi_dn*kernel_loc*weight;
    end
  end
  int_loc = ((b-a)/2)*((yupper_right-ylower_right)/2)*int_loc;

  % Compute upper triangle
  base = (b-a); height = (yupper_left - yupper_right); area = base*height/2;

  xl = a + (1/3)*base; yl = yupper_right + (1/3)*height;
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (-27/48)*area*dphi_dn*kernel_loc;

  xl = a + (3/5)*base; yl = yupper_right + (1/5)*height;
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

  xl = a + (1/5)*base; yl = yupper_right + (1/5)*height;
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

  xl = a + (1/5)*base; yl = yupper_right + (3/5)*height;
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

  % 3 pt Gaussian quadrature
  % xl = a + (1/6)*(b-a); yl = yupper_right + (1/6)*(yupper_left - yupper_right); %scatter(xl,yl,'kx');
  % R = sqrt( (x-xl)^2 - (y-yl)^2 );
  % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

  % xl = a + (1/6)*(b-a); yl = yupper_right + (2/3)*(yupper_left - yupper_right); %scatter(xl,yl,'kx');
  % R = sqrt( (x-xl)^2 - (y-yl)^2 );
  % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

  % xl = a + (2/3)*(b-a); yl = yupper_right + (1/6)*(yupper_left - yupper_right); %scatter(xl,yl,'kx');
  % R = sqrt( (x-xl)^2 - (y-yl)^2 );
  % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

  % Compute lower triangle
  base = (b-a); height = (ylower_right - ylower_left); area = base*height/2;

  xl = a + (1/3)*base; yl = ylower_right - (1/3)*height;
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (-27/48)*area*dphi_dn*kernel_loc;

  xl = a + (3/5)*base; yl = ylower_right - (1/5)*height;
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

  xl = a + (1/5)*base; yl = ylower_right - (1/5)*height;
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

  xl = a + (1/5)*base; yl = ylower_right - (3/5)*height;
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

  % 3 pt Gaussian quadrature for triangles
  %xl = a + (1/6)*(b-a); yl = ylower_right - (1/6)*(ylower_right - ylower_left); %scatter(xl,yl,'kx');
  %R = sqrt( (x-xl)^2 - (y-yl)^2 );
  %dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  %dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  %kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

  %xl = a + (1/6)*(b-a); yl = ylower_right - (2/3)*(ylower_right - ylower_left); %scatter(xl,yl,'kx');
  %R = sqrt( (x-xl)^2 - (y-yl)^2 );
  %dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  %dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  %kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

  %xl = a + (2/3)*(b-a); yl = ylower_right - (1/6)*(ylower_right - ylower_left); %scatter(xl,yl,'kx'); pause;
  %R = sqrt( (x-xl)^2 - (y-yl)^2 );
  %dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  %dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  %kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

end
