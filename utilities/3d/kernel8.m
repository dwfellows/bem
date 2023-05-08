function [int_loc] = kernel8(x,y,bndry_pts,lambda,mu,flow,amp,omega,Uinf)

  %pts = [-1/sqrt(3), 1/sqrt(3)];
  pts = [-sqrt(3/5), 0, sqrt(3/5)];
  weights = [5/9, 8/9, 5/9];

  % boundary element corner values
  a = bndry_pts(1,1); b = bndry_pts(2,1);
  c = bndry_pts(1,2); d = bndry_pts(3,2);

  cent = [(a+b)/2, (c+d)/2];  % boundary element center coordinates

  La = flow.La; Lb = flow.Lb;

  % figure;
  % scatter(bndry_pts(:,1), bndry_pts(:,2), 'b'); hold on;
  % scatter(x,y,'m','filled');

  int_loc = 0;
  if (c > y)

    xleft = x - (d-y);
    xright = x - (c-y);
    %scatter(xleft,d,'r');
    %scatter(xright,c,'g');

    % Compute left rectangle
    for l = 1:3
      for m = 1:3
        x1 = ((xleft-a)/2)*pts(l) + ((xleft+a)/2);
        y1 = ((d-c)/2)*pts(m) + ((d+c)/2);
        weight = weights(l)*weights(m);
        %scatter(x1,y1,'kx');

        R = sqrt( (x-x1)^2 - (y-y1)^2 );
        kernel_loc = cos( lambda*R/cos(mu) ) / R;
        dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(x1*La/tan(mu),y1*La,La,Lb) + (1/2)*amp*omega*flow.psi(x1*La/tan(mu),y1*La,La,Lb);
        dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* x1 ./ Uinf) .* dphi_dn;
        int_loc = int_loc + dphi_dn*kernel_loc*weight;
      end
    end
    int_loc = ((xleft-a)*(d-c)/4)*int_loc;

    % Compute right triangle
    base = (xright-xleft); height = (d-c); area = base*height/2;

    xl = xleft + (1/3)*base; yl = c + (1/3)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (-27/48)*area*dphi_dn*kernel_loc;

    xl = xleft + (3/5)*base; yl = c + (1/5)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

    xl = xleft + (1/5)*base; yl = c + (1/5)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

    xl = xleft + (1/5)*base; yl = c + (3/5)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;


    % 3pt Gaussian quadrature rule for triangles
    %xl = xleft + (1/6)*(xright-xleft); yl = c + (1/6)*(d-c); %scatter(xl,yl,'kx');
    %R = sqrt( (x-xl)^2 - (y-yl)^2 );
    %dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    %dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    %kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

    %xl = xleft + (2/3)*(xright-xleft); yl = c + (1/6)*(d-c); %scatter(xl,yl,'kx');
    %R = sqrt( (x-xl)^2 - (y-yl)^2 );
    %dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    %dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    %kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

    %xl = xleft + (1/6)*(xright-xleft); yl = c + (2/3)*(d-c); %scatter(xl,yl,'kx'); pause;
    %R = sqrt( (x-xl)^2 - (y-yl)^2 );
    %dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    %dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    %kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

  else % Element below

    xright = x - (y-c);
    xleft = x - (y-d);
    %scatter(xright,c,'g'); scatter(xleft,d,'r');

    % Compute left rectangle
    for l = 1:3
      for m = 1:3
        x1 = ((xleft-a)/2)*pts(l) + ((xleft+a)/2);
        y1 = ((d-c)/2)*pts(m) + ((d+c)/2);
        weight = weights(l)*weights(m);
        %scatter(x1,y1,'kx');

        R = sqrt( (x-x1)^2 - (y-y1)^2 );
        kernel_loc = cos( lambda*R/cos(mu) ) / R;
        dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(x1*La/tan(mu),y1*La,La,Lb) + (1/2)*amp*omega*flow.psi(x1*La/tan(mu),y1*La,La,Lb);
        dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* x1 ./ Uinf) .* dphi_dn;
        int_loc = int_loc + dphi_dn*kernel_loc*weight;
      end
    end
    int_loc = ((xleft-a)*(d-c)/4)*int_loc;

    % Compute right triangle
    base = (xright-xleft); height = (d-c); area = base*height/2;

    xl = xleft + (1/3)*base; yl = c + (1/3)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (-27/48)*area*dphi_dn*kernel_loc;

    xl = xleft + (3/5)*base; yl = c + (1/5)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

    xl = xleft + (1/5)*base; yl = c + (1/5)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

    xl = xleft + (1/5)*base; yl = c + (3/5)*height;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;


    % % 3pt Gaussian quadrature for triangles
    % xl = xleft + (1/6)*(xright-xleft); yl = c + (1/6)*(d-c); %scatter(xl,yl,'kx');
    % R = sqrt( (x-xl)^2 - (y-yl)^2 );
    % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

    % xl = xleft + (2/3)*(xright-xleft); yl = c + (1/6)*(d-c); %scatter(xl,yl,'kx');
    % R = sqrt( (x-xl)^2 - (y-yl)^2 );
    % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

    % xl = xleft + (1/6)*(xright-xleft); yl = c + (2/3)*(d-c); %scatter(xl,yl,'kx'); pause;
    % R = sqrt( (x-xl)^2 - (y-yl)^2 );
    % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

  end

end
