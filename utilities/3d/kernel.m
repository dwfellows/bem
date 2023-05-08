function [int_loc] = kernel(x,y,bndry_pts,lambda,mu,flow,amp,omega,Uinf)

  pts = [-1/sqrt(3), 1/sqrt(3)];

  a = bndry_pts(1,1); b = bndry_pts(2,1);
  c = bndry_pts(1,2); d = bndry_pts(3,2);

  cent = [(a+b)/2, (c+d)/2];
  La = flow.La; Lb = flow.Lb;

  int_loc = 0; int_num = 0;
  for l = 1:2
    for m = 1:2
      x1 = ((b-a)/2)*pts(l) + ((a+b)/2);
      y1 = ((d-c)/2)*pts(m) + ((c+d)/2); 

      R = sqrt( (x-x1)^2 - (y-y1)^2 );
      kernel_loc = cos( lambda*R/cos(mu) ) / R;

      dphi_dn = -(i/2)*amp*Uinf*flow.dpsi(x1*La/tan(mu),y1*La,La,Lb) + (1/2)*amp*omega*flow.psi(x1*La/tan(mu),y1*La,La,Lb);
      dphi_dn = La .* exp(i*omega*(La/(cos(mu)*sin(mu))) .* x1 ./ Uinf) .* dphi_dn;


      % if ( (x > x1) && (y1 >= (y - sqrt((x1-x)^2))) && (y1 <= (y + sqrt((x1-x)^2)))  )
      %   R = sqrt( (x-x1)^2 - (y-y1)^2 );
      %   kernel_loc = cos( lambda*R/cos(mu) ) / R;
      %   int_num = int_num + 1;
      % else
      %   kernel_loc = 0;
      %   int_loc = 0; return
      % end

      int_loc = int_loc + dphi_dn*kernel_loc;
    end
  end

  %int_loc = ((b-a)/2)*((d-c)/2)*int_loc*(int_num/4);
  int_loc = ((b-a)/2)*((d-c)/2)*int_loc;

end
