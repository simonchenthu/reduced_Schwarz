function [u,q] = reduced_Schwarz_elliptic_l2(x_nw_o,y_nw_o,x_se_o,y_se_o,...
                                          Dic,Ind,knn,dx)

%%    Initialization

[Mx,My] = size(Dic); 
Ny_patch = (y_se_o-y_nw_o)/dx+1; 
Nx_patch = (x_se_o-x_nw_o)/dx+1; 
N_dict = cellfun(@(x)size(x,3),Dic);
Nx = length(x_nw_o(1):dx:x_se_o(end)); Ny = length(y_nw_o(1):dx:y_se_o(end));

% restriction
Rst1 = cell(Mx,My); % S
Rst2 = cell(Mx,My); % E
Rst3 = cell(Mx,My); % N
Rst4 = cell(Mx,My); % W
for k = 2:My
    for j = 1:Mx
        Rst1{j,k} = reshape(Dic{j,k}((y_se_o(k-1)-y_nw_o(k))/dx+1,1:end-1,:),[],N_dict(j,k));
    end
end
for k = 1:My
    for j = 1:Mx-1
        Rst2{j,k} = reshape(Dic{j,k}(2:end-1,end-(x_se_o(j)-x_nw_o(j+1))/dx,:),[],N_dict(j,k));
    end
end
for k = 1:My-1
    for j = 1:Mx
        Rst3{j,k} = reshape(Dic{j,k}(end-(y_se_o(k)-y_nw_o(k+1))/dx,:,:),[],N_dict(j,k));
    end
end
for k = 1:My
    for j = 2:Mx
        Rst4{j,k} = reshape(Dic{j,k}(2:end,(x_se_o(j-1)-x_nw_o(j))/dx+1,:),[],N_dict(j,k));
    end
end




% intial boundary on each patch
% Patch_bdy = cell(Mx,My); 
Patch_bdy = cellfun(@(x)size(x,1),Ind,'UniformOutput',false);
Patch_bdy = cellfun(@(x)zeros(x,1),Patch_bdy,'UniformOutput',false);                   
for j = 1:Mx
    Patch_bdy{j,1}(1:Nx_patch(j)) = Ind{j,1}(1:Nx_patch(j),1);
    Patch_bdy{j,My}(Nx_patch(j)+Ny_patch(My):2*Nx_patch(j)-2+Ny_patch(My)) ...
        = Ind{j,My}(Nx_patch(j)+Ny_patch(My):2*Nx_patch(j)-2+Ny_patch(My),1);
end
for k = 1:My
    Patch_bdy{1,k}(2*Nx_patch(1)-1+Ny_patch(k):end) ...
        = Ind{1,k}(2*Nx_patch(1)-1+Ny_patch(k):end,1);
    Patch_bdy{Mx,k}(Nx_patch(Mx)+1:Nx_patch(Mx)+Ny_patch(k)-1) ...
        = Ind{Mx,k}(Nx_patch(Mx)+1:Nx_patch(Mx)+Ny_patch(k)-1,1);
end
Patch_bdy_new = Patch_bdy;






% index for neighbors and c
I_knn = zeros(knn,Mx,My);
coef = zeros(knn-1,Mx,My);

%% Schwarz iteration
res = 1; q = 0;
tol = 1e-5;
while res > tol && q+1<=800
    
    q = q+1;
    
    for k = 1:My
        for j = 1:Mx
            
            %%%   Search for k-nearest neighbors
            distance = sqrt(dx*sum((Ind{j,k}-Patch_bdy{j,k}).^2));
            [~,I_knn(:,j,k)] = mink(distance,knn);
            
            %%%   Linear Approximation
            Bdy_c = Ind{j,k}(:,I_knn(1,j,k));              % use nearest point as center
            Bdy_n = Ind{j,k}(:,I_knn(2:knn,j,k));          % use the other knn as basis
            
            Bdy_n_centered = Bdy_n-Bdy_c;
            coef(:,j,k) = Bdy_n_centered\(Patch_bdy{j,k}-Bdy_c);
            
        end
    end
    
    
    %%  Update Boundaries
    % update patches on the south (north BC)
    for k = 2:My
        for j = 1:Mx
            
            Patch_bdy_new{j,k-1}(Nx_patch(j)+Ny_patch(k-1):2*Nx_patch(j)-2+Ny_patch(k-1))...
                                        = Rst1{j,k}(:,I_knn(1,j,k))...
                                          +(Rst1{j,k}(:,I_knn(2:knn,j,k))-Rst1{j,k}(:,I_knn(1,j,k)))*coef(:,j,k);
            
        end
    end
    
    % update patches on the east (west BC)
    for k = 1:My
        for j = 1:Mx-1
            Patch_bdy_new{j+1,k}(2*Nx_patch(j+1)-1+Ny_patch(k):end)...
                                        = Rst2{j,k}(:,I_knn(1,j,k))...
                                          +(Rst2{j,k}(:,I_knn(2:knn,j,k))-Rst2{j,k}(:,I_knn(1,j,k)))*coef(:,j,k);
        end
    end
    
    % update patches on the north (south BC)
    for k = 1:My-1
        for j = 1:Mx
            Patch_bdy_new{j,k+1}(1:Nx_patch(j))...
                                        = Rst3{j,k}(:,I_knn(1,j,k))...
                                          +(Rst3{j,k}(:,I_knn(2:knn,j,k))-Rst3{j,k}(:,I_knn(1,j,k)))*coef(:,j,k);
        end
    end
    
    % update patches on the west (east BC)
    for k = 1:My
        for j = 2:Mx
            Patch_bdy_new{j-1,k}(Nx_patch(j-1)+1:Nx_patch(j-1)+Ny_patch(k)-1)...
                                        = Rst4{j,k}(:,I_knn(1,j,k))...
                                          +(Rst4{j,k}(:,I_knn(2:knn,j,k))-Rst4{j,k}(:,I_knn(1,j,k)))*coef(:,j,k);
        end
    end
    
    res = cellfun(@(x,y)x-y,Patch_bdy,Patch_bdy_new,'UniformOutput',false);
    res = ones(1,Mx)*sqrt(dx*cellfun(@(x)x'*x,res))*ones(My,1);

    Patch_bdy = Patch_bdy_new;
    
end





%% Output

% Partition of Unity  
POU = Partition_of_Unity(x_nw_o,y_nw_o,x_se_o,y_se_o,Nx,Ny);

% global solution
u = zeros(Ny,Nx);
for k = 1:My
    for j = 1:Mx
        
        u_center = Dic{j,k}(:,:,I_knn(1,j,k));
        u_center = u_center(:);                              % offline solution at the nearest point
        
        u_neighbor = Dic{j,k}(:,:,I_knn(2:knn,j,k));
        u_neighbor = reshape(u_neighbor,[],knn-1);                      % offline solutions at the other knn
        
        u_n_centered = u_neighbor-u_center*ones(1,knn-1);
        
        u_patch = u_center+u_n_centered*coef(:,j,k); u_patch = reshape(u_patch,Ny_patch(k),Nx_patch(j),[]);
        
        u(y_nw_o(k)/dx+1:y_se_o(k)/dx+1,x_nw_o(j)/dx+1:x_se_o(j)/dx+1) ...
                  = u(y_nw_o(k)/dx+1:y_se_o(k)/dx+1,x_nw_o(j)/dx+1:x_se_o(j)/dx+1)...
                    + u_patch.*POU(y_nw_o(k)/dx+1:y_se_o(k)/dx+1,x_nw_o(j)/dx+1:x_se_o(j)/dx+1,j,k);
        
    end
end


end