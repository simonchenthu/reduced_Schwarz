function err = err_inf(u,v)
% This function computes the reltive L-infinity error of the two matrices u
% v (same size) with v being the reference

err = max(max(abs(u-v)))/max(max(abs(v)));

end