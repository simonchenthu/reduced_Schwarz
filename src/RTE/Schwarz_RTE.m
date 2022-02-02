function [f,theta,iter_residual,q] = Schwarz_RTE(epsilon,sigma_x,...
                                                 f_bdy,theta_bdy,...
                                                 t,s,dx,v0,w0)

Mx = length(t); Nx = length(t(1):dx:s(end));

Nv = length(v0);

w_bdy = [w0';1;1];

i_update_l = round((s(1:end-1)-t(2:end))/dx)+1; i_update_l = [0,i_update_l];
i_update_r = round((s(1:end-1)-t(2:end))/dx); i_update_r = [i_update_r,0];

%%% Iteration

Patch_bdy = ones(Nv+2,Mx); 
Patch_bdy(Nv/2+1:Nv,1) = f_bdy(Nv/2+1:Nv); Patch_bdy(Nv+1,1) = theta_bdy(1);
Patch_bdy(1:Nv/2,Mx) = f_bdy(1:Nv/2); Patch_bdy(end,Mx) = theta_bdy(2);
Patch_bdy_temp = Patch_bdy;   % 1:Nv/2...v<0; Nv/2+1:Nv...v>0; Nv+1...x=t; Nv+2...x=s

iter_residual = 1; 
q = 0;


f_patch = cell(1,Mx); theta_patch = cell(1,Mx); 
while iter_residual>=1e-3 && q+1<=800
    q = q+1;
    
    for j = 1:Mx
        
        f_bdy_temp = Patch_bdy(1:Nv,j); theta_bdy_temp = Patch_bdy(Nv+1:end,j);
        sigma_x_patch = sigma_x(t(j)/dx+1:s(j)/dx+1);
        x_patch = t(j):dx:s(j);
        
        [f,theta] = nonlinearRTE_FD_aa_monotone(f_bdy_temp,theta_bdy_temp,epsilon,...
                                                        sigma_x_patch,x_patch,w0,v0,dx); 
        
        f_patch{1,j} = f; theta_patch{1,j} = theta';                                            
                                                    
        if j~=Mx
            Patch_bdy_temp(Nv/2+1:Nv,j+1) = f(Nv/2+1:Nv,end-i_update_r(j));                            
            Patch_bdy_temp(Nv+1,j+1) = theta(end-i_update_r(j));   
        end
        
        if j~=1
            Patch_bdy_temp(1:Nv/2,j-1) = f(1:Nv/2,i_update_l(j));
            Patch_bdy_temp(Nv+2,j-1) = theta(i_update_l(j));
        end
    end
        
%     iter_residual = sum(max(Patch_bdy_temp-Patch_bdy));
    iter_residual = sum(w_bdy'*(Patch_bdy_temp-Patch_bdy).^2);
    Patch_bdy = Patch_bdy_temp;
                         
end


%%% Output

%%% Partition of Unity 

POU_f = Partition_of_Unity(t,s,dx,Nv);
POU_theta = Partition_of_Unity(t,s,dx,1);

f = zeros(Nv,Nx);
theta = zeros(Nx,1);
for j = 1:Mx
    
    f_temp = f_patch{1,j}; theta_temp = theta_patch{1,j};
   
    f(:,round(t(j)/dx)+1:round(s(j)/dx)+1) = f(:,round(t(j)/dx)+1:round(s(j)/dx)+1) + f_temp.*POU_f(:,round(t(j)/dx)+1:round(s(j)/dx)+1,j);
    theta(round(t(j)/dx)+1:round(s(j)/dx)+1) = theta(round(t(j)/dx)+1:round(s(j)/dx)+1) + theta_temp.*POU_theta(round(t(j)/dx)+1:round(s(j)/dx)+1,j);
    
end

end