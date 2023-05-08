function [int_loc] = kernel6(x,y,bndry_pts,lambda,mu,flow,amp,omega,Uinf)

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

  yupper = y + (x-a);
  ylower = y - (x-a);

  % figure;
  % scatter(bndry_pts(:,1), bndry_pts(:,2), 'b'); hold on;
  % scatter(x,y,'m','filled');
  % scatter(a,yupper,'r'); scatter(a,ylower,'r');

  % Compute upper triangle
  base = (x-a); height = (yupper - y); area = base*height/2;

  xl = a + (1/3)*base; yl = y + (1/3)*height;
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (-27/48)*area*dphi_dn*kernel_loc;

  xl = a + (3/5)*base; yl = y + (1/5)*height;
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

  xl = a + (1/5)*base; yl = y + (1/5)*height;
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

  xl = a + (1/5)*base; yl = y + (3/5)*height;
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;


  % 3pt Gaussian quadrature for triangles
  %xl = a + (1/6)*(x-a); yl = y + (1/6)*(yupper - y); %scatter(xl,yl,'kx');
  %R = sqrt( (x-xl)^2 - (y-yl)^2 );
  %dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  %dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  %kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

  %xl = a + (1/6)*(x-a); yl = y + (2/3)*(yupper - y); %scatter(xl,yl,'kx');
  %R = sqrt( (x-xl)^2 - (y-yl)^2 );
  %dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  %dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  %kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

  %xl = a + (2/3)*(x-a); yl = y + (1/6)*(yupper - y); %scatter(xl,yl,'kx');
  %R = sqrt( (x-xl)^2 - (y-yl)^2 );
  %dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  %dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  %kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;


  % Compute lower triangle
  base = (x-a); height = (y - ylower); area = base*height/2;

  xl = a + (1/3)*base; yl = y - (1/3)*height;
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (-27/48)*area*dphi_dn*kernel_loc;

  xl = a + (3/5)*base; yl = y - (1/5)*height;
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

  xl = a + (1/5)*base; yl = y - (1/5)*height;
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

  xl = a + (1/5)*base; yl = y - (3/5)*height;
  R = sqrt( (x-xl)^2 - (y-yl)^2 );
  dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;


  % 3pt Gaussian quadrature for triangles
  % xl = a + (1/6)*(x-a); yl = y - (1/6)*(y - ylower); %scatter(xl,yl,'kx');
  % R = sqrt( (x-xl)^2 - (y-yl)^2 );
  % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

  % xl = a + (1/6)*(x-a); yl = y - (2/3)*(y - ylower); %scatter(xl,yl,'kx');
  % R = sqrt( (x-xl)^2 - (y-yl)^2 );
  % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

  % xl = a + (2/3)*(x-a); yl = y - (1/6)*(y - ylower); %scatter(xl,yl,'kx'); pause;
  % R = sqrt( (x-xl)^2 - (y-yl)^2 );
  % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
  % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
  % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

end
