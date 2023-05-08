function [int_loc] = kernel2(x,y,bndry_pts,lambda,mu,flow,amp,omega,Uinf)

  %pts = [-1/sqrt(3), 1/sqrt(3)];
  pts = [-sqrt(3/5), 0, sqrt(3/5)];
  weights = [5/9, 8/9, 5/9];

  % boundary element corner values
  a = bndry_pts(1,1); b = bndry_pts(2,1);
  c = bndry_pts(1,2); d = bndry_pts(3,2);

  cent = [(a+b)/2, (c+d)/2];  % boundary element center coordinates
  La = flow.La; Lb = flow.Lb;

  %figure;
  %scatter(bndry_pts(:,1), bndry_pts(:,2), 'b'); hold on;
  %scatter(x,y,'m','filled');

  int_loc = 0;
  if (c > y) % Element above collocation point

    ylower = y + (x-b); %scatter(b,ylower,'g');
    yupper = y + (x-a); %scatter(a,yupper,'r');

    %figure;
    %scatter(bndry_pts(:,1), bndry_pts(:,2), 'b'); hold on;
    %scatter(x,y,'m','filled');
    %scatter(b,ylower,'g'); scatter(a,yupper,'r');

    % Compute enclosed rectangle
    for l = 1:3
      for m = 1:3
        x1 = ((b-a)/2)*pts(l) + ((a+b)/2);
        y1 = ((ylower-c)/2)*pts(m) + ((c+ylower)/2); 
        weight = weights(l)*weights(m);
        %scatter(x1,y1,'kx');

        R = sqrt( (x-x1)^2 - (y-y1)^2 );
        kernel_loc = cos( lambda*R/cos(mu) ) / R;
        dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(x1*La/tan(mu),y1*La,La,Lb) + (1/2)*amp*omega*flow.psi(x1*La/tan(mu),y1*La,La,Lb);
        dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* x1 ./ Uinf) .* dphi_dn;
        int_loc = int_loc + dphi_dn*kernel_loc*weight;
      end
    end
    int_loc = ((b-a)*(ylower-c)/4)*int_loc;

    % Compute enclosed triangle
    base = (b-a); height = (yupper-ylower); area = base*height/2;

    xl = a + (1/3)*base; yl = ylower + (1/3)*height; %scatter(xl,yl,'kx');
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (-27/48)*area*dphi_dn*kernel_loc;

    xl = a + (3/5)*base; yl = ylower + (1/5)*height; %scatter(xl,yl,'kx');
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

    xl = a + (1/5)*base; yl = ylower + (1/5)*height; %scatter(xl,yl,'kx');
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

    xl = a + (1/5)*base; yl = ylower + (3/5)*height; %scatter(xl,yl,'kx'); pause;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;


    % 3pt Gaussian quadrature for triangles
    % xl = a + (1/6)*(b-a); yl = ylower + (1/6)*(yupper-ylower); %scatter(xl,yl,'kx');
    % R = sqrt( (x-xl)^2 - (y-yl)^2 );
    % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

    % xl = a + (2/3)*(b-a); yl = ylower + (1/6)*(yupper-ylower); %scatter(xl,yl,'kx');
    % R = sqrt( (x-xl)^2 - (y-yl)^2 );
    % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

    % xl = a + (1/6)*(b-a); yl = ylower + (2/3)*(yupper-ylower); %scatter(xl,yl,'kx'); pause;
    % R = sqrt( (x-xl)^2 - (y-yl)^2 );
    % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

  else % Element below collocation point

    yupper = y - (x-b);
    ylower = y - (x-a);

    %figure;
    %scatter(bndry_pts(:,1), bndry_pts(:,2), 'b'); hold on;
    %scatter(x,y,'m','filled');
    %scatter(b,yupper,'g'); scatter(a,ylower,'r');

    % Compute enclosed rectangle
    for l = 1:3
      for m = 1:3
        x1 = ((b-a)/2)*pts(l) + ((a+b)/2);
        y1 = ((d-yupper)/2)*pts(m) + ((d+yupper)/2); 
        weight = weights(l)*weights(m);
        %scatter(x1,y1,'kx');

        R = sqrt( (x-x1)^2 - (y-y1)^2 );
        kernel_loc = cos( lambda*R/cos(mu) ) / R;
        dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(x1*La/tan(mu),y1*La,La,Lb) + (1/2)*amp*omega*flow.psi(x1*La/tan(mu),y1*La,La,Lb);
        dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* x1 ./ Uinf) .* dphi_dn;

        int_loc = int_loc + dphi_dn*kernel_loc*weight;
      end
    end
    int_loc = ((b-a)*(d-yupper)/4)*int_loc;

    base = (b-a); height = (yupper-ylower); area = base*height/2;

    % 4pt Gaussian quadrature for triangles
    xl = a + (1/3)*base; yl = yupper - (1/3)*height; %scatter(xl,yl,'kx');
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (-27/48)*area*dphi_dn*kernel_loc;

    xl = a + (3/5)*base; yl = yupper - (1/5)*height; %scatter(xl,yl,'kx');
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

    xl = a + (1/5)*base; yl = yupper - (1/5)*height; %scatter(xl,yl,'kx');
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;

    xl = a + (1/5)*base; yl = yupper - (3/5)*height; %scatter(xl,yl,'kx'); pause;
    R = sqrt( (x-xl)^2 - (y-yl)^2 );
    dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (25/48)*area*dphi_dn*kernel_loc;



    % % 3pt Gaussian quadrature for triangles
    % xl = a + (1/6)*(b-a); yl = yupper - (1/6)*(yupper-ylower); %scatter(xl,yl,'kx');
    % R = sqrt( (x-xl)^2 - (y-yl)^2 );
    % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

    % xl = a + (2/3)*(b-a); yl = yupper - (1/6)*(yupper-ylower); %scatter(xl,yl,'kx');
    % R = sqrt( (x-xl)^2 - (y-yl)^2 );
    % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

    % xl = a + (1/6)*(b-a); yl = yupper - (2/3)*(yupper-ylower); %scatter(xl,yl,'kx'); pause;
    % R = sqrt( (x-xl)^2 - (y-yl)^2 );
    % dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(xl*La/tan(mu),yl*La,La,Lb) + (1/2)*amp*omega*flow.psi(xl*La/tan(mu),yl*La,La,Lb);
    % dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* xl ./ Uinf) .* dphi_dn;
    % kernel_loc = cos( lambda*R/cos(mu) ) / R; int_loc = int_loc + (area/3)*dphi_dn*kernel_loc;

  end

end
