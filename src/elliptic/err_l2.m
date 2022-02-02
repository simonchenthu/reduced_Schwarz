function err = err_l2(u,v)
% This function computes the reltive L2 error of the two matrices (same 
% size)with v being the reference

err = sqrt(sum(sum((u-v).^2))/sum(sum(v.^2)));

end