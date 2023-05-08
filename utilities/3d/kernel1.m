function [int_loc] = kernel1(x,y,bndry_pts,lambda,mu,flow,amp,omega,Uinf)

  % boundary element corner values
  a = bndry_pts(1,1); b = bndry_pts(2,1);
  c = bndry_pts(1,2); d = bndry_pts(3,2);

  cent = [(a+b)/2, (c+d)/2];  % boundary element center coordinates

  il = [2/3, 1/6, 1/6; 1/6, 2/3, 1/6; 1/6, 1/6, 2/3];
  La = flow.La; Lb = flow.Lb;

  int_loc = 0;
  if (c > y) % Element above collocation point

    xb = x - (c-y);
    yb = y + (x-a);

    % Use 4 pt Gaussian quadrature for triangles
    base = (xb-a); height = (yb-c); area = base*height/2;
    xl = a + (1/3)*base; yl = c + (1/3)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (-27/48)*area*dphi_dn*kernel_loc;

    xl = a + (3/5)*base; yl = c + (1/5)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

    xl = a + (1/5)*base; yl = c + (1/5)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

    xl = a + (1/5)*base; yl = c + (3/5)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

    % 3 pt Gaussian quadrature for triangles
    % xl = a + (1/6)*(xb-a); yl = c + (1/6)*(yb-c);
    % R = sqrt( (x-xl)^2 - (y-yl)^2 );
    % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + dphi_dn*kernel_loc;

    % xl = a + (2/3)*(xb-a); yl = c + (1/6)*(yb-c);
    % R = sqrt( (x-xl)^2 - (y-yl)^2 );
    % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + dphi_dn*kernel_loc;

    % xl = a + (1/6)*(xb-a); yl = c + (2/3)*(yb-c);
    % R = sqrt( (x-xl)^2 - (y-yl)^2 );
    % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + dphi_dn*kernel_loc;

    %int_loc = ((xb-a)*(yb-c)/6)*int_loc;

  else % Element below collocation point

    xb = x - (y-d);
    yb = y - (x-a);

    base = (xb-a); height = (d-yb); area = base*height/2;

    % Use 4pt Gaussian quadrature for triangles
    xl = a + (1/3)*base; yl = d - (1/3)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (-27/48)*area*dphi_dn*kernel_loc;

    xl = a + (3/5)*base; yl = d - (1/5)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

    xl = a + (1/5)*base; yl = d - (1/5)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

    xl = a + (1/5)*base; yl = d - (3/5)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;


    % xl = a + (1/6)*(xb-a); yl = d - (1/6)*(d-yb);
    % R = sqrt( (x-xl)^2 - (y-yl)^2 );
    % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + dphi_dn*kernel_loc;

    % xl = a + (2/3)*(xb-a); yl = d - (1/6)*(d-yb);
    % R = sqrt( (x-xl)^2 - (y-yl)^2 );
    % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + dphi_dn*kernel_loc;

    % xl = a + (1/6)*(xb-a); yl = d - (2/3)*(d-yb);
    % R = sqrt( (x-xl)^2 - (y-yl)^2 );
    % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + dphi_dn*kernel_loc;

    % int_loc = ((xb-a)*(d-yb)/6)*int_loc;

  end

end
