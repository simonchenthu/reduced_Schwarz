function [u,q] = Schwarz_elliptic(f,del_f,a,x_nw_o,y_nw_o,x_se_o,y_se_o,...
                                  bdy_w,bdy_e,bdy_s,bdy_n,dx)
% Classical Schwarz iteration for semilinear elliptic equations

%%    Initialization

Mx = size(x_nw_o,2); My = size(y_nw_o,2);
Ny_patch = (y_se_o-y_nw_o)/dx+1; 
Nx_patch = (x_se_o-x_nw_o)/dx+1; 
Nx = length(x_nw_o(1):dx:x_se_o(end)); Ny = length(y_nw_o(1):dx:y_se_o(end));

% restriction indices
j_update_s = (y_se_o(1:end-1)-y_nw_o(2:end))/dx+1; j_update_s = [0,j_update_s];
i_update_e = (x_se_o(1:end-1)-x_nw_o(2:end))/dx; i_update_e = [i_update_e,0];
j_update_n = (y_se_o(1:end-1)-y_nw_o(2:end))/dx; j_update_n = [j_update_n,0];
i_update_w = (x_se_o(1:end-1)-x_nw_o(2:end))/dx+1; i_update_w = [0,i_update_w];




% intial boundary on each patch
Patch_bdy_s = Nx_patch'*ones(1,My); Patch_bdy_s = num2cell(Patch_bdy_s);
Patch_bdy_s = cellfun(@(x)zeros(x,1),Patch_bdy_s,'UniformOutput',false);

Patch_bdy_e = ones(Mx,1)*Ny_patch; Patch_bdy_e = num2cell(Patch_bdy_e);
Patch_bdy_e = cellfun(@(x)zeros(x,1),Patch_bdy_e,'UniformOutput',false);

Patch_bdy_n = Nx_patch'*ones(1,My); Patch_bdy_n = num2cell(Patch_bdy_n);
Patch_bdy_n = cellfun(@(x)zeros(x,1),Patch_bdy_n,'UniformOutput',false);

Patch_bdy_w = ones(Mx,1)*Ny_patch; Patch_bdy_w = num2cell(Patch_bdy_w);
Patch_bdy_w = cellfun(@(x)zeros(x,1),Patch_bdy_w,'UniformOutput',false);

for j = 1:Mx
    Patch_bdy_s{j,1} = bdy_s(x_nw_o(j)/dx+1:x_se_o(j)/dx+1);
    Patch_bdy_n{j,My} = bdy_n(x_nw_o(j)/dx+1:x_se_o(j)/dx+1);
end
for k = 1:My
    Patch_bdy_w{1,k} = bdy_w(y_nw_o(k)/dx+1:y_se_o(k)/dx+1);
    Patch_bdy_e{Mx,k} = bdy_e(y_nw_o(k)/dx+1:y_se_o(k)/dx+1);
end
Patch_bdy_s_new = Patch_bdy_s;
Patch_bdy_e_new = Patch_bdy_e;
Patch_bdy_n_new = Patch_bdy_n;
Patch_bdy_w_new = Patch_bdy_w;



%% Schwarz iteration
res = 1; q = 0;
tol = 1e-5;
u_patch = cell(Mx,My);
while res > tol && q+1<=800 % q = 1:iter
    
    q = q+1;
    
    for k = 1:My
        for j = 1:Mx
            
            x_patch = x_nw_o(j):dx:x_se_o(j); 
            y_patch = y_nw_o(k):dx:y_se_o(k);
            
            u_patch{j,k} = semilinear_elliptic_newton(x_patch,y_patch,dx,f,del_f,a,...
                                         Patch_bdy_w{j,k},Patch_bdy_e{j,k},...
                                         Patch_bdy_s{j,k},Patch_bdy_n{j,k});
            
        end
    end
    
    
    %%%  Update Boundaries
    % update patches on the south (north BC)
    for k = 2:My
        for j = 1:Mx
            Patch_bdy_n_new{j,k-1} = u_patch{j,k}(j_update_s(k),:)';
        end
    end
    
    % update patches on the east (west BC)
    for k = 1:My
        for j = 1:Mx-1
            Patch_bdy_w_new{j+1,k} = u_patch{j,k}(:,end-i_update_e(j));
        end
    end
    
    % update patches on the north (south BC)
    for k = 1:My-1
        for j = 1:Mx
            Patch_bdy_s_new{j,k+1} = u_patch{j,k}(end-j_update_n(k),:)';
        end
    end
    
    % update patches on the west (east BC)
    for k = 1:My
        for j = 2:Mx
            Patch_bdy_e_new{j-1,k} = u_patch{j,k}(:,i_update_w(j));
        end
    end
    
    res_s = cellfun(@(x,y)x-y,Patch_bdy_s_new,Patch_bdy_s,'UniformOutput',false);
    res_e = cellfun(@(x,y)x-y,Patch_bdy_e_new,Patch_bdy_e,'UniformOutput',false);
    res_n = cellfun(@(x,y)x-y,Patch_bdy_n_new,Patch_bdy_n,'UniformOutput',false);
    res_w = cellfun(@(x,y)x-y,Patch_bdy_w_new,Patch_bdy_w,'UniformOutput',false);
    
    res = ones(1,Mx)*sqrt(dx*(cellfun(@(x)x'*x,res_s)+cellfun(@(x)x'*x,res_e)...
             +cellfun(@(x)x'*x,res_n)+cellfun(@(x)x'*x,res_w)))*ones(My,1);
     
    Patch_bdy_s = Patch_bdy_s_new;
    Patch_bdy_e = Patch_bdy_e_new;
    Patch_bdy_n = Patch_bdy_n_new;
    Patch_bdy_w = Patch_bdy_w_new;
    
end



%% Output

% Partition of Unity  
POU = Partition_of_Unity(x_nw_o,y_nw_o,x_se_o,y_se_o,Nx,Ny);

% Global solution
u = zeros(Ny,Nx);
for k = 1:My
    for j = 1:Mx
        
        u(y_nw_o(k)/dx+1:y_se_o(k)/dx+1,x_nw_o(j)/dx+1:x_se_o(j)/dx+1) ...
                  = u(y_nw_o(k)/dx+1:y_se_o(k)/dx+1,x_nw_o(j)/dx+1:x_se_o(j)/dx+1)...
                    + u_patch{j,k}.*POU(y_nw_o(k)/dx+1:y_se_o(k)/dx+1,x_nw_o(j)/dx+1:x_se_o(j)/dx+1,j,k);
        
    end
end


end