function [int_loc] = kernel3(x,y,bndry_pts,lambda,mu,flow,amp,omega,Uinf)

  %pts = [-1/sqrt(3), 1/sqrt(3)];
  pts = [-sqrt(3/5), 0, sqrt(3/5)];
  weights = [5/9, 8/9, 5/9];

  % boundary element corner values
  a = bndry_pts(1,1); b = bndry_pts(2,1);
  c = bndry_pts(1,2); d = bndry_pts(3,2);

  cent = [(a+b)/2, (c+d)/2];  % boundary element center coordinates
  La = flow.La; Lb = flow.Lb;

  int_loc = 0;
  if (c > y) % Element above collocation point

    yb = y + (x-b);
    xb = x - (d-y);

    % Compute left rectangle
    for l = 1:3
      for m = 1:3
        x1 = ((xb-a)/2)*pts(l) + ((a+xb)/2);
        y1 = ((d-c)/2)*pts(m) + ((c+d)/2);
        weight = weights(l)*weights(m);

        R = sqrt( (x-x1)^2 - (y-y1)^2 );
        kernel_loc = cos( lambda*R/cos(mu) ) / R;
        dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(x1*La/tan(mu),y1*La,La,Lb) + (1/2)*amp*omega*flow.psi(x1*La/tan(mu),y1*La,La,Lb);
        dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* x1 ./ Uinf) .* dphi_dn;

        int_loc = int_loc + dphi_dn*kernel_loc*weight;
      end
    end
    int_loc = ((xb-a)*(d-c)/4)*int_loc;

    % Compute right rectangle
    for l = 1:3
      for m = 1:3
        x1 = ((b-xb)/2)*pts(l) + ((b+xb)/2);
        y1 = ((yb-c)/2)*pts(m) + ((yb+c)/2);
        weight = weights(l)*weights(m);

        R = sqrt( (x-x1)^2 - (y-y1)^2 );
        kernel_loc = cos( lambda*R/cos(mu) ) / R;
        dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(x1*La/tan(mu),y1*La,La,Lb) + (1/2)*amp*omega*flow.psi(x1*La/tan(mu),y1*La,La,Lb);
        dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* x1 ./ Uinf) .* dphi_dn;

        int_loc = int_loc + ((b-xb)*(yb-c)/4)*dphi_dn*kernel_loc*weight;
      end
    end

    % Compute enclosed triangle
    base = (b-xb); height = (d-yb); area = base*height/2;

    % 4pt Gaussian quadrature for triangles
    xl = xb + (1/3)*base; yl = yb + (1/3)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (-27/48)*area*dphi_dn*kernel_loc;

    xl = xb + (3/5)*base; yl = yb + (1/5)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

    xl = xb + (1/5)*base; yl = yb + (1/5)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

    xl = xb + (1/5)*base; yl = yb + (3/5)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;


    % 3pt Gaussian quadrature for triangles
    % xl = xb + (1/6)*(b-xb); yl = yb + (1/6)*(d-yb);
    % R = sqrt( (x-xl)^2 - (y-yl)^2 );
    % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

    % xl = xb + (2/3)*(b-xb); yl = yb + (1/6)*(d-yb);
    % R = sqrt( (x-xl)^2 - (y-yl)^2 );
    % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

    % xl = xb + (1/6)*(b-xb); yl = yb + (2/3)*(d-yb);
    % R = sqrt( (x-xl)^2 - (y-yl)^2 );
    % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

  else % Element below collocation point

    xb = x - (y-c);
    yb = y - (x-b);

    % Compute left rectangle
    for l = 1:3
      for m = 1:3
        x1 = ((xb-a)/2)*pts(l) + ((a+xb)/2);
        y1 = ((d-c)/2)*pts(m) + ((c+d)/2); 
        weight = weights(l)*weights(m);

        R = sqrt( (x-x1)^2 - (y-y1)^2 );
        kernel_loc = cos( lambda*R/cos(mu) ) / R;
        dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(x1*La/tan(mu),y1*La,La,Lb) + (1/2)*amp*omega*flow.psi(x1*La/tan(mu),y1*La,La,Lb);
        dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* x1 ./ Uinf) .* dphi_dn;
        int_loc = int_loc + dphi_dn*kernel_loc*weight;
      end
    end
    int_loc = ((xb-a)*(d-c)/4)*int_loc;

    % Compute right rectangle
    for l = 1:3
      for m = 1:3
        x1 = ((b-xb)/2)*pts(l) + ((b+xb)/2);
        y1 = ((d-yb)/2)*pts(m) + ((yb+d)/2);
        weight = weights(l)*weights(m);

        R = sqrt( (x-x1)^2 - (y-y1)^2 );
        kernel_loc = cos( lambda*R/cos(mu) ) / R;
        dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(x1*La/tan(mu),y1*La,La,Lb) + (1/2)*amp*omega*flow.psi(x1*La/tan(mu),y1*La,La,Lb);
        dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* x1 ./ Uinf) .* dphi_dn;

        int_loc = int_loc + ((b-xb)*(d-yb)/4)*dphi_dn*weight*kernel_loc;
      end
    end

    base = (b-xb); height = (yb-c); area = base*height/2;

    xl = xb + (1/3)*base; yl = yb - (1/3)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (-27/48)*area*dphi_dn*kernel_loc;

    xl = xb + (3/5)*base; yl = yb - (1/5)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

    xl = xb + (1/5)*base; yl = yb - (1/5)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

    xl = xb + (1/5)*base; yl = yb - (3/5)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;


    % 3pt Gaussian quadrature for triangles
    % xl = xb + (1/6)*(b-xb); yl = yb - (1/6)*(yb-c);
    % R = sqrt( (x-xl)^2 - (y-yl)^2 );
    % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

    % xl = xb + (2/3)*(b-xb); yl = yb - (1/6)*(yb-c);
    % R = sqrt( (x-xl)^2 - (y-yl)^2 );
    % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

    % xl = xb + (1/6)*(b-xb); yl = yb - (2/3)*(yb-c);
    % R = sqrt( (x-xl)^2 - (y-yl)^2 );
    % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

  end

end
