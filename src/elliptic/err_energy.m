function err = err_energy(u,v,ax,ay)
% This function computes the reltive energy error of the two matrices (same 
% size)with v being the reference. a is the permeability.

dux = (u(:,2:end)-u(:,1:end-1));
duy = (u(2:end,:)-u(1:end-1,:));

dvx = (v(:,2:end)-v(:,1:end-1));
dvy = (v(2:end,:)-v(1:end-1,:));

err = sqrt( ( sum(sum(ax.*(dux-dvx).^2)) + sum(sum(ay.*(duy-dvy).^2)) )...
    /( sum(sum(ax.*dvx.^2)) + sum(sum(ay.*dvy.^2)) ) );

end