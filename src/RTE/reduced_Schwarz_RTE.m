function [f,theta,q,discrep,iter_res] = reduced_Schwarz_RTE(t,s,Dic,Ind,knn,dx,w_bdy) %,v)


[Nv,N_sample,Mx] = size(Ind); Nv = Nv-2; Nx = length(t(1):dx:s(end));
Rst = zeros(Nv/2+1,N_sample,2,Mx); % 
for m = 2:Mx
    
    f_temp = Dic{1,m}; theta_temp = Dic{2,m};
    
    Rst(:,:,1,m) = [reshape(f_temp(1:Nv/2,round((s(m-1)-t(m))/dx)+1,:),Nv/2,N_sample);...
                    theta_temp(round((s(m-1)-t(m))/dx)+1,:)];
    
end

for m = 1:Mx-1
    
    f_temp = Dic{1,m}; theta_temp = Dic{2,m};
    
    Rst(:,:,2,m) = [reshape(f_temp(Nv/2+1:Nv,end-round((s(m)-t(m+1))/dx),:),Nv/2,N_sample);...
                    theta_temp(end-round((s(m)-t(m+1))/dx),:)];
    
end




%%% Iteration

Patch_bdy = ones(Nv+2,Mx); 
Patch_bdy(Nv/2+1:Nv,1) = Ind(Nv/2+1:Nv,1,1); Patch_bdy(Nv+1,1) = Ind(Nv+1,1,1);
Patch_bdy(1:Nv/2,Mx) = Ind(1:Nv/2,1,Mx); Patch_bdy(end,Mx) = Ind(end,1,Mx);
Patch_bdy_temp = Patch_bdy;   % 1:Nv/2...v<0; Nv/2+1:Nv...v>0; Nv+1...x=t; Nv+2...x=s

I_knn = zeros(knn,Mx);
delta = zeros(knn,Mx);
coef = zeros(knn-1,Mx);

q_max = 800;

discrep = zeros(q_max);
iter_res = zeros(q_max);

% iter_res = 1; 
iter_res_new = 1; 
q = 0;

while iter_res_new>=1e-3 && q+1<=q_max % && iter_res_new<r*iter_res && min(rcondn)>eps
    q = q+1; % iter_res = iter_res_new;
    
    for j = 1:Mx
        
        %%%  compute c
        
        
        %%%   Search for k-nearest neighbors   %%%
        distance = sqrt(w_bdy'*(Ind(:,:,j)-Patch_bdy(:,j)*ones(1,N_sample)).^2);
        % distance = max(Ind(:,:,j)-Patch_bdy(:,j)*ones(1,N_sample));
        [delta(:,j),I_knn(:,j)] = mink(distance,knn);
        
        %%%    Local Linear Approximation    %%%
        Bdy_c = Ind(:,I_knn(1,j),j);              % use nearest point as center
        Bdy_n = Ind(:,I_knn(2:knn,j),j);          % use the other knn as basis
        
        Bdy_n_centered = Bdy_n-Bdy_c*ones(1,knn-1);
        
        
        
%         A = [-Bdy_n_centered,-ones(Nv+2,1);...
%               Bdy_n_centered,-ones(Nv+2,1)];
%         b = [-Patch_bdy(:,j)+Bdy_c;Patch_bdy(:,j)-Bdy_c];
%         options = optimoptions('linprog','Display','off');
%         [coef_temp,~,~] = linprog([zeros(1,knn-1),1],A,b,[],[],[],[],options);
%         coef(:,j) = coef_temp(1:end-1);
        b = sqrt(w_bdy).*(Patch_bdy(:,j)-Bdy_c);
        A = sqrt(w_bdy).*Bdy_n_centered;
%         coef(:,j) = A\b;
        coef(:,j) = (A'*A+(1e-4)*eye(knn-1))\(A'*b);
        
        discrep(q) = discrep(q) + sqrt(ones(1,Nv+2)*(A*coef(:,j)-b).^2);
%         coef(:,j) = (Bdy_n_centered'*(w_bdy.*Bdy_n_centered))\(Bdy_n_centered'*(w_bdy.*(Patch_bdy(:,j)-Bdy_c)));
        

    end
    
    
    
%%%     Update Boundary Conditions   
    for j = 1:Mx-1
       
        Patch_bdy_temp(Nv/2+1:Nv,j+1) = Rst(1:Nv/2,I_knn(1,j),2,j)...
            +(Rst(1:Nv/2,I_knn(2:end,j),2,j)-Rst(1:Nv/2,I_knn(1,j),2,j)*ones(1,knn-1))*coef(:,j);
        Patch_bdy_temp(Nv+1,j+1) = Rst(end,I_knn(1,j),2,j)...
            +(Rst(end,I_knn(2:end,j),2,j)-Rst(end,I_knn(1,j),2,j)*ones(1,knn-1))*coef(:,j);
        
    end
    
    for j = 2:Mx
            
        Patch_bdy_temp(1:Nv/2,j-1) = Rst(1:Nv/2,I_knn(1,j),1,j)...
            +(Rst(1:Nv/2,I_knn(2:end,j),1,j)-Rst(1:Nv/2,I_knn(1,j),1,j)*ones(1,knn-1))*coef(:,j);
        Patch_bdy_temp(Nv+2,j-1) = Rst(end,I_knn(1,j),1,j)...
            +(Rst(end,I_knn(2:end,j),1,j)-Rst(end,I_knn(1,j),1,j)*ones(1,knn-1))*coef(:,j);

    end  
        
        
%     iter_residual = sum(max(Patch_bdy_temp-Patch_bdy));
    iter_res_new = sum(w_bdy'*(Patch_bdy_temp-Patch_bdy).^2);
    Patch_bdy = Patch_bdy_temp;
                         
    iter_res(q) = iter_res_new;
end



        


%%% Output
iter_res = iter_res(1:q);
discrep = discrep(1:q);

%%% Partition of Unity 

POU_f = Partition_of_Unity(t,s,dx,Nv);
POU_theta = Partition_of_Unity(t,s,dx,1);

f = zeros(Nv,Nx);
theta = zeros(Nx,1);
for j = 1:Mx
    
    f_temp = Dic{1,j}; theta_temp = Dic{2,j};
    
    f_center = f_temp(:,:,I_knn(1,j)); 
    f_center = f_center(:);                              % offline solution at the nearest point
    theta_center = theta_temp(:,I_knn(1,j));                       
    
    f_neighbor = f_temp(:,:,I_knn(2:knn,j)); 
    f_neighbor = reshape(f_neighbor,[],knn-1);                      % offline solutions at the other knn
    theta_neighbor = theta_temp(:,I_knn(2:knn,j));
    
    f_n_centered = f_neighbor-f_center*ones(1,knn-1);
    theta_n_centered = theta_neighbor-theta_center*ones(1,knn-1);
    
    f_patch = f_center+f_n_centered*coef(:,j); f_patch = reshape(f_patch,Nv,[]);
    theta_patch = theta_center+theta_n_centered*coef(:,j);
    
    
    
    
    % test
%     x_patch = t(j):dx:s(j);
%     [xx_patch,vv_patch] = meshgrid(x_patch,v);
%     
%     figure(2)
%     if j==1
%         clf;
%     end
%     subplot(2,1,1); plot(x_patch,theta_patch); xlim([t(1),s(end)]); ylim([0,4]); hold on;
%     % title(['Local Linear Approximate Solution, Iteration = ',int2str(q)]);hold on;
%     subplot(2,1,2); mesh(xx_patch,vv_patch,f_patch); xlim([t(1),s(end)]); zlim([0,20]); hold on;
%     %title([int2str(knn),' nearest neighbors are used']);
%     %pause(0.05);
%     pause;
    %
    
    
    
    
    f(:,round(t(j)/dx)+1:round(s(j)/dx)+1) = f(:,round(t(j)/dx)+1:round(s(j)/dx)+1) + f_patch.*POU_f(:,round(t(j)/dx)+1:round(s(j)/dx)+1,j);
    theta(round(t(j)/dx)+1:round(s(j)/dx)+1) = theta(round(t(j)/dx)+1:round(s(j)/dx)+1) + theta_patch.*POU_theta(round(t(j)/dx)+1:round(s(j)/dx)+1,j);
    
end


end