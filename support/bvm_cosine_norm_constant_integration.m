function [] = bvm_cosine_norm_constant_integration(k1,k2,k3)

  %fun = @(x,y) 1./( sqrt(x + y) .* (1 + x + y).^2 )
  fun = @(x,y) exp( k1 * cos(x) + k2 * cos(y) - k3 * cos(x-y)); 

  xlow = 0;
  xhigh = 2*pi;
  ylow = 0;
  yhigh = 2*pi;

  q = integral2(fun,xlow,xhigh,ylow,yhigh)

end
