function err = err_h1(u,v,h)
% This function computes the reltive H1 error of the two matrices (same 
% size)with v being the reference

dux = (u(:,2:end)-u(:,1:end-1))/h;
duy = (u(2:end,:)-u(1:end-1,:))/h;

dvx = (v(:,2:end)-v(:,1:end-1))/h;
dvy = (v(2:end,:)-v(1:end-1,:))/h;

l2_err = sum(sum((u-v).^2));
grad_err = sum(sum((dux-dvx).^2)) + sum(sum((duy-dvy).^2));

h1_norm = sum(sum(v.^2));

err = sqrt((l2_err+grad_err)/h1_norm);

end