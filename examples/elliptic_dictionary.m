% This code generate the dictionaries on each patch for the semilinear
% elliptic equation
clear;
% addpath(genpath('../src/'))
addpath('../src/elliptic/')
%% Domain parameters
Lx = 1.0; Nx = 2^(9); Dx_buffer = 2^(-4); Dx_overlap = 2^(-4); 
dx = Lx/Nx; x0 = 0:dx:Lx; Nx = length(x0); Mx = 4; 
Ly = 1.0; Ny = 2^(9); Dy_buffer = 2^(-4); Dy_overlap = 2^(-4); 
dy = Ly/Ny; y0 = 0:dy:Ly; Ny = length(y0); My = 4; 

plotz = 2;
expmt = 1;

x_nw_o = 0:Lx/Mx:Lx-Lx/Mx; x_nw_o = max(x_nw_o-Dx_overlap,0);
y_nw_o = 0:Ly/My:Ly-Ly/My; y_nw_o = max(y_nw_o-Dy_overlap,0);
x_se_o = Lx/Mx:Lx/Mx:Lx; x_se_o = min(x_se_o+Dx_overlap,Ly);
y_se_o = Ly/My:Ly/My:Ly; y_se_o = min(y_se_o+Dy_overlap,Ly);


x_nw_b = 0:Lx/Mx:Lx-Lx/Mx; x_nw_b = max(x_nw_b-Dx_overlap-Dx_buffer,0);
y_nw_b = 0:Ly/My:Ly-Ly/My; y_nw_b = max(y_nw_b-Dy_overlap-Dy_buffer,0);
x_se_b = Lx/Mx:Lx/Mx:Lx; x_se_b = min(x_se_b+Dx_overlap+Dx_buffer,Ly);
y_se_b = Ly/My:Ly/My:Ly; y_se_b = min(y_se_b+Dy_overlap+Dy_buffer,Ly);

%% Equation parameters
n = 4;
epsilon = 2^(-n);

f = @(u) u.^3;
del_f = @(u) 3*u.^2;   
f_no = 1;

a = @(x,y) 2+sin(2*pi*x).*cos(2*pi*y)...
    +(2+1.8*sin(2*pi*x/epsilon))./(2+1.8*cos(2*pi*y/epsilon))...
    +(2+sin(2*pi*y/epsilon))./(2+1.8*cos(2*pi*x/epsilon));
a_n = 1;


%% Sampling parameters
N_dict = 64;

bdy_D_n = @(x) sin(2*pi*x); bdy_no_n = 1;
bdy_D_s = @(x) -sin(2*pi*x); bdy_no_s = 1;
bdy_D_w = @(y) sin(2*pi*y); bdy_no_w = 1;
bdy_D_e = @(y) -sin(2*pi*y);  bdy_no_e = 1;

radius_n = 20;
          
dim_r = 5;
          
          
          
%% Generate dictionaries

% rng('default');
rng('shuffle');
for k = 1:My
    for j = 1:Mx
        
        t_start = tic;
        
        % x/y range for buffered patches
        x_patch_b = x_nw_b(j):dx:x_se_b(j);
        y_patch_b = y_nw_b(k):dx:y_se_b(k);
        
        Nx_patch_b = length(x_patch_b);
        Ny_patch_b = length(y_patch_b);
        
        %% sample boundary conditions
        if j~=1 && j~=Mx && k~=1 && k~=My
            
            % compute weights W (order: SENW) and linear transform C
            bdy_x = [x_patch_b,x_patch_b(end)*ones(1,Ny_patch_b-1),...
                     fliplr(x_patch_b(1:end-1)),x_patch_b(1)*ones(1,Ny_patch_b-2)]; 

            bdy_y = [y_patch_b(1)*ones(1,Nx_patch_b),y_patch_b(2:end),...
                     y_patch_b(end)*ones(1,Nx_patch_b-1),fliplr(y_patch_b(2:end-1))]; 
            
            disq =(bdy_x-bdy_x').^2+(bdy_y-bdy_y').^2; 
            
            W = 2*dx^2./disq; W(isnan(W)|isinf(W)) = 0;
            W = diag( ones(1,2*(Nx_patch_b-1)+2*(Ny_patch_b-1))*W + dx ) - W;
            
            C = chol(W);
            
            % generate random samples
            rho = rand(1,N_dict); rho = nthroot(rho,dim_r);
            
            % Generate unit vectors
            x = randn(2*(Nx_patch_b-1)+2*(Ny_patch_b-1),N_dict); 
            x_norm = sqrt(ones(1,2*(Nx_patch_b-1)+2*(Ny_patch_b-1))*x.^2);
            x = C\x;
            x_unit = abs(x./x_norm);
        
            % Boundary conditions
            bdy_rand = radius_n*rho.*x_unit;
        
            bdy_s = bdy_rand(1:Nx_patch_b,:);
            bdy_e = bdy_rand(Nx_patch_b:(Nx_patch_b-1)+(Ny_patch_b-1)+1,:);
            bdy_n = flipud(bdy_rand((Nx_patch_b-1)+(Ny_patch_b-1)+1:2*(Nx_patch_b-1)+(Ny_patch_b-1)+1,:));
            bdy_w = [bdy_rand(1,:);...
                     flipud(bdy_rand(2*(Nx_patch_b-1)+(Ny_patch_b-1)+1:end,:))];
            
        elseif j==1 && k~=1 && k~=My  % west
            
            % compute weights
            bdy_x = [x_patch_b(2:end),x_patch_b(end)*ones(1,Ny_patch_b-1),...
                     fliplr(x_patch_b(2:end-1)),x_patch_b(1)*ones(1,Ny_patch_b)]; 

            bdy_y = [y_patch_b(1)*ones(1,Nx_patch_b-1),y_patch_b(2:end),...
                     y_patch_b(end)*ones(1,Nx_patch_b-2),fliplr(y_patch_b)]; 
            
            disq =(bdy_x-bdy_x').^2+(bdy_y-bdy_y').^2; 
            
            W = 2*dx^2./disq; W(isnan(W)|isinf(W)) = 0;
            W = diag( ones(1,2*(Nx_patch_b-1)+2*(Ny_patch_b-1))*W + dx ) - W;
            
            W11 = W(1:end-Ny_patch_b,1:end-Ny_patch_b);
            W12 = W(1:end-Ny_patch_b,end-Ny_patch_b+1:end);
            W22 = W(end-Ny_patch_b+1:end,end-Ny_patch_b+1:end);
            
            bdy_D_patch = fliplr(bdy_D_w(y_patch_b))';
            r = sqrt(radius_n^2 - bdy_D_patch'*(W22-W12'*(W11\W12))*bdy_D_patch);
            
            C = chol(W11);
            
            % generate random samples
            rho = rand(1,N_dict); rho = nthroot(rho,dim_r);
            
            % Generate unit vectors
            x = randn(2*(Nx_patch_b-1)+(Ny_patch_b-1)-1,N_dict); 
            x_norm = sqrt(ones(1,2*(Nx_patch_b-1)+(Ny_patch_b-1)-1)*x.^2);
            x = C\x;
            x_unit = abs(x./x_norm);
            
            % Boundary conditions
            bdy_rand = r*rho.*x_unit;
            
            bdy_s = [bdy_D_w(y_patch_b(1))*ones(1,N_dict);bdy_rand(1:Nx_patch_b-1,:)];
            bdy_e = bdy_rand(Nx_patch_b-1:(Nx_patch_b-1)+(Ny_patch_b-1),:);
            bdy_n = [bdy_D_w(y_patch_b(end))*ones(1,N_dict);...
                     flipud(bdy_rand((Nx_patch_b-1)+(Ny_patch_b-1):end,:))];
            bdy_w = repmat(bdy_D_w(y_patch_b)',1,N_dict);
                 
        elseif j==Mx && k~=1 && k~=My % east
            
            % compute weights
            bdy_x = [fliplr(x_patch_b(1:end-1)),x_patch_b(1)*ones(1,Ny_patch_b-1),...
                     x_patch_b(2:end-1),x_patch_b(end)*ones(1,Ny_patch_b)]; 

            bdy_y = [y_patch_b(end)*ones(1,Nx_patch_b-1),flipud(y_patch_b(1:end-1)),...
                     y_patch_b(1)*ones(1,Nx_patch_b-2),y_patch_b]; 
            
            disq =(bdy_x-bdy_x').^2+(bdy_y-bdy_y').^2; 
            
            W = 2*dx^2./disq; W(isnan(W)|isinf(W)) = 0;
            W = diag( ones(1,2*(Nx_patch_b-1)+2*(Ny_patch_b-1))*W + dx ) - W;
            
            W11 = W(1:end-Ny_patch_b,1:end-Ny_patch_b);
            W12 = W(1:end-Ny_patch_b,end-Ny_patch_b+1:end);
            W22 = W(end-Ny_patch_b+1:end,end-Ny_patch_b+1:end);
            
            bdy_D_patch = bdy_D_e(y_patch_b)';
            r = sqrt(radius_n^2 - bdy_D_patch'*(W22-W12'*(W11\W12))*bdy_D_patch);
            
            C = chol(W11);
            
            % generate random samples
            rho = rand(1,N_dict); rho = nthroot(rho,dim_r);
            
            % Generate unit vectors
            x = randn(2*(Nx_patch_b-1)+(Ny_patch_b-1)-1,N_dict); 
            x_norm = sqrt(ones(1,2*(Nx_patch_b-1)+(Ny_patch_b-1)-1)*x.^2);
            x = C\x;
            x_unit = abs(x./x_norm);
            
            % Boundary conditions
            bdy_rand = r*rho.*x_unit;
            
            bdy_n = [flipud(bdy_rand(1:Nx_patch_b-1,:));...
                     bdy_D_e(y_patch_b(end))*ones(1,N_dict)];
            bdy_w = flipud(bdy_rand(Nx_patch_b-1:(Nx_patch_b-1)+(Ny_patch_b-1),:));
            bdy_s = [bdy_rand((Nx_patch_b-1)+(Ny_patch_b-1):end,:);...
                     bdy_D_e(y_patch_b(1))*ones(1,N_dict)];
            bdy_e = repmat(bdy_D_e(y_patch_b)',1,N_dict);
            
        elseif k==1 && j~=1 && j~=Mx % south
            
            % compute weights
            bdy_x = [x_patch_b(end)*ones(1,Ny_patch_b-1),fliplr(x_patch_b(1:end-1)),...
                     x_patch_b(1)*ones(1,Ny_patch_b-2),x_patch_b]; 
            
            bdy_y = [y_patch_b(2:end),y_patch_b(end)*ones(1,Nx_patch_b-1),...
                     fliplr(y_patch_b(2:end-1)),y_patch_b(1)*ones(1,Nx_patch_b)]; 
            
            disq =(bdy_x-bdy_x').^2+(bdy_y-bdy_y').^2; 
            
            W = 2*dx^2./disq; W(isnan(W)|isinf(W)) = 0;
            W = diag( ones(1,2*(Nx_patch_b-1)+2*(Ny_patch_b-1))*W + dx ) - W;
            
            W11 = W(1:end-Nx_patch_b,1:end-Nx_patch_b);
            W12 = W(1:end-Nx_patch_b,end-Nx_patch_b+1:end);
            W22 = W(end-Nx_patch_b+1:end,end-Nx_patch_b+1:end);
            
            bdy_D_patch = bdy_D_s(x_patch_b)';
            r = sqrt(radius_n^2 - bdy_D_patch'*(W22-W12'*(W11\W12))*bdy_D_patch);
            
            C = chol(W11);
            
            % generate random samples
            rho = rand(1,N_dict); rho = nthroot(rho,dim_r);
            
            % Generate unit vectors
            x = randn((Nx_patch_b-1)+2*(Ny_patch_b-1)-1,N_dict); 
            x_norm = sqrt(ones(1,(Nx_patch_b-1)+2*(Ny_patch_b-1)-1)*x.^2);
            x = C\x;
            x_unit = abs(x./x_norm);
            
            % Boundary conditions
            bdy_rand = r*rho.*x_unit;
        
            bdy_e = [bdy_D_s(x_patch_b(end))*ones(1,N_dict);bdy_rand(1:Ny_patch_b-1,:)];
            bdy_n = flipud(bdy_rand(Ny_patch_b-1:(Ny_patch_b-1)+(Nx_patch_b-1),:));
            bdy_w = [bdy_D_s(x_patch_b(1))*ones(1,N_dict);...
                     flipud(bdy_rand((Ny_patch_b-1)+(Nx_patch_b-1):end,:))];
            bdy_s = repmat(bdy_D_s(x_patch_b)',1,N_dict);
            
        elseif k==My && j~=1 && j~=Mx % north
            
            % compute weights
            bdy_x = [x_patch_b(1)*ones(1,Ny_patch_b-1),x_patch_b(2:end),...
                     x_patch_b(end)*ones(1,Ny_patch_b-2),fliplr(x_patch_b)]; 

            bdy_y = [fliplr(y_patch_b(2:end)),y_patch_b(1)*ones(1,Nx_patch_b-1),...
                     y_patch_b(2:end-1),y_patch_b(end)*ones(1,Nx_patch_b)]; 
            
            disq =(bdy_x-bdy_x').^2+(bdy_y-bdy_y').^2; 
            
            W = 2*dx^2./disq; W(isnan(W)|isinf(W)) = 0;
            W = diag( ones(1,2*(Nx_patch_b-1)+2*(Ny_patch_b-1))*W + dx ) - W;
            
            W11 = W(1:end-Nx_patch_b,1:end-Nx_patch_b);
            W12 = W(1:end-Nx_patch_b,end-Nx_patch_b+1:end);
            W22 = W(end-Nx_patch_b+1:end,end-Nx_patch_b+1:end);
            
            bdy_D_patch = fliplr(bdy_D_n(x_patch_b))';
            r = sqrt(radius_n^2 - bdy_D_patch'*(W22-W12'*(W11\W12))*bdy_D_patch);
            
            C = chol(W11);
            
            % generate random samples
            rho = rand(1,N_dict); rho = nthroot(rho,dim_r);
            
            % Generate unit vectors
            x = randn((Nx_patch_b-1)+2*(Ny_patch_b-1)-1,N_dict); 
            x_norm = sqrt(ones(1,(Nx_patch_b-1)+2*(Ny_patch_b-1)-1)*x.^2);
            x = C\x;
            x_unit = abs(x./x_norm);
            
            % Boundary conditions
            bdy_rand = r*rho.*x_unit;
            
            bdy_w = [flipud(bdy_rand(1:Ny_patch_b-1,:));...
                     bdy_D_n(x_patch_b(1))*ones(1,N_dict)];
            bdy_s = bdy_rand(Ny_patch_b-1:(Ny_patch_b-1)+(Nx_patch_b-1),:);
            bdy_e = [bdy_rand((Ny_patch_b-1)+(Nx_patch_b-1):end,:);...
                     bdy_D_n(x_patch_b(end))*ones(1,N_dict)];
            bdy_n = repmat(bdy_D_n(x_patch_b)',1,N_dict);
            
        elseif j==1 && k==1 % south-west
            
            % compute weights SENW
            bdy_x = [x_patch_b(end)*ones(1,Ny_patch_b-1),...
                     fliplr(x_patch_b(2:end-1)),...
                     x_patch_b(1)*ones(1,Ny_patch_b),x_patch_b(2:end)]; 

            bdy_y = [y_patch_b(2:end),...
                     y_patch_b(end)*ones(1,Nx_patch_b-2),...
                     fliplr(y_patch_b),y_patch_b(1)*ones(1,Nx_patch_b-1)]; 
            
            disq =(bdy_x-bdy_x').^2+(bdy_y-bdy_y').^2; 
            
            W = 2*dx^2./disq; W(isnan(W)|isinf(W)) = 0;
            W = diag( ones(1,2*(Nx_patch_b-1)+2*(Ny_patch_b-1))*W + dx ) - W;
            
            W11 = W(1:end-Nx_patch_b-Ny_patch_b+1,1:end-Nx_patch_b-Ny_patch_b+1);
            W12 = W(1:end-Nx_patch_b-Ny_patch_b+1,end-Nx_patch_b-Ny_patch_b+2:end);
            W22 = W(end-Nx_patch_b-Ny_patch_b+2:end,end-Nx_patch_b-Ny_patch_b+2:end);
            
            bdy_D_patch = [fliplr(bdy_D_w(y_patch_b))';bdy_D_s(x_patch_b(2:end))'];
            r = sqrt(radius_n^2 - bdy_D_patch'*(W22-W12'*(W11\W12))*bdy_D_patch);
            
            C = chol(W11);
            
            % generate random samples
            rho = rand(1,N_dict); rho = nthroot(rho,dim_r);
            
            % Generate unit vectors
            x = randn((Nx_patch_b-1)+(Ny_patch_b-1)-1,N_dict); 
            x_norm = sqrt(ones(1,(Nx_patch_b-1)+(Ny_patch_b-1)-1)*x.^2);
            x = C\x;
            x_unit = abs(x./x_norm);
        
            % Boundary conditions
            bdy_rand = r*rho.*x_unit;
        
            bdy_e = [bdy_D_s(x_patch_b(end))*ones(1,N_dict);...
                     bdy_rand(1:Ny_patch_b-1,:)];
            bdy_n = [bdy_D_w(y_patch_b(end))*ones(1,N_dict);...
                     flipud(bdy_rand(Ny_patch_b-1:end,:))];
            bdy_w = repmat(bdy_D_w(y_patch_b'),1,N_dict);
            bdy_s = repmat(bdy_D_s(x_patch_b'),1,N_dict);
            
        elseif j==Mx && k==1 % south-east
            
            % compute weights SENW
            bdy_x = [fliplr(x_patch_b(2:end)),...
                     x_patch_b(1)*ones(1,Ny_patch_b-2),...
                     x_patch_b,x_patch_b(end)*ones(1,Ny_patch_b-1)]; 

            bdy_y = [y_patch_b(end)*ones(1,Nx_patch_b-1),...
                     fliplr(y_patch_b(2:end-1)),...
                     y_patch_b(1)*ones(1,Nx_patch_b),y_patch_b(2:end)]; 
            
            disq =(bdy_x-bdy_x').^2+(bdy_y-bdy_y').^2; 
            
            W = 2*dx^2./disq; W(isnan(W)|isinf(W)) = 0;
            W = diag( ones(1,2*(Nx_patch_b-1)+2*(Ny_patch_b-1))*W + dx ) - W;
            
            W11 = W(1:end-Nx_patch_b-Ny_patch_b+1,1:end-Nx_patch_b-Ny_patch_b+1);
            W12 = W(1:end-Nx_patch_b-Ny_patch_b+1,end-Nx_patch_b-Ny_patch_b+2:end);
            W22 = W(end-Nx_patch_b-Ny_patch_b+2:end,end-Nx_patch_b-Ny_patch_b+2:end);
            
            bdy_D_patch = [bdy_D_s(x_patch_b)';bdy_D_e(y_patch_b(2:end))'];
            r = sqrt(radius_n^2 - bdy_D_patch'*(W22-W12'*(W11\W12))*bdy_D_patch);
            
            C = chol(W11);
            
            % generate random samples
            rho = rand(1,N_dict); rho = nthroot(rho,dim_r);
            
            % Generate unit vectors
            x = randn((Nx_patch_b-1)+(Ny_patch_b-1)-1,N_dict); 
            x_norm = sqrt(ones(1,(Nx_patch_b-1)+(Ny_patch_b-1)-1)*x.^2);
            x = C\x;
            x_unit = abs(x./x_norm);
        
            % Boundary conditions
            bdy_rand = r*rho.*x_unit;
        
            bdy_n = [flipud(bdy_rand(1:Nx_patch_b-1,:));...
                     bdy_D_e(y_patch_b(end))*ones(1,N_dict)];
            bdy_w = [bdy_D_s(x_patch_b(1))*ones(1,N_dict);...
                     flipud(bdy_rand(Nx_patch_b-1:end,:))];
            bdy_s = repmat(bdy_D_s(x_patch_b'),1,N_dict);
            bdy_e = repmat(bdy_D_e(y_patch_b'),1,N_dict);
            
        elseif j==Mx && k==My % north-east
            
            % compute weights SENW
            bdy_x = [x_patch_b(1)*ones(1,Ny_patch_b-1),...
                     x_patch_b(2:end-1),...
                     x_patch_b(end)*ones(1,Ny_patch_b),fliplr(x_patch_b(1:end-1))]; 

            bdy_y = [flipud(y_patch_b(2:end)),...
                     y_patch_b(1)*ones(1,Nx_patch_b-2),...
                     y_patch_b,y_patch_b(end)*ones(1,Nx_patch_b-1)]; 
            
            disq =(bdy_x-bdy_x').^2+(bdy_y-bdy_y').^2; 
            
            W = 2*dx^2./disq; W(isnan(W)|isinf(W)) = 0;
            W = diag( ones(1,2*(Nx_patch_b-1)+2*(Ny_patch_b-1))*W + dx ) - W;
            
            W11 = W(1:end-Nx_patch_b-Ny_patch_b+1,1:end-Nx_patch_b-Ny_patch_b+1);
            W12 = W(1:end-Nx_patch_b-Ny_patch_b+1,end-Nx_patch_b-Ny_patch_b+2:end);
            W22 = W(end-Nx_patch_b-Ny_patch_b+2:end,end-Nx_patch_b-Ny_patch_b+2:end);
            
            bdy_D_patch = [bdy_D_e(y_patch_b)';fliplr(bdy_D_n(x_patch_b(1:end-1)))'];
            r = sqrt(radius_n^2 - bdy_D_patch'*(W22-W12'*(W11\W12))*bdy_D_patch);
            
            C = chol(W11);
            
            % generate random samples
            rho = rand(1,N_dict); rho = nthroot(rho,dim_r);
            
            % Generate unit vectors
            x = randn((Nx_patch_b-1)+(Ny_patch_b-1)-1,N_dict); 
            x_norm = sqrt(ones(1,(Nx_patch_b-1)+(Ny_patch_b-1)-1)*x.^2);
            x = C\x;
            x_unit = abs(x./x_norm);
        
            % Boundary conditions
            bdy_rand = r*rho.*x_unit;
            
            bdy_w = [flipud(bdy_rand(1:Ny_patch_b-1,:));...
                     bdy_D_n(x_patch_b(1))*ones(1,N_dict)];
            bdy_s = [bdy_rand(Ny_patch_b-1:end,:);...
                     bdy_D_e(y_patch_b(1))*ones(1,N_dict)];
            bdy_e = repmat(bdy_D_e(y_patch_b'),1,N_dict);
            bdy_n = repmat(bdy_D_n(x_patch_b'),1,N_dict);
            
        elseif j==1 && k==My % north-west
            
            % compute weights SENW
            bdy_x = [x_patch_b(2:end),...
                     x_patch_b(end)*ones(1,Ny_patch_b-2),...
                     fliplr(x_patch_b),x_patch_b(1)*ones(1,Ny_patch_b-1)]; 

            bdy_y = [y_patch_b(1)*ones(1,Nx_patch_b-1),...
                     y_patch_b(2:end-1),...
                     y_patch_b(end)*ones(1,Nx_patch_b),fliplr(y_patch_b(1:end-1))]; 
            
            disq =(bdy_x-bdy_x').^2+(bdy_y-bdy_y').^2; 
            
            W = 2*dx^2./disq; W(isnan(W)|isinf(W)) = 0;
            W = diag( ones(1,2*(Nx_patch_b-1)+2*(Ny_patch_b-1))*W + dx ) - W;
            
            W11 = W(1:end-Nx_patch_b-Ny_patch_b+1,1:end-Nx_patch_b-Ny_patch_b+1);
            W12 = W(1:end-Nx_patch_b-Ny_patch_b+1,end-Nx_patch_b-Ny_patch_b+2:end);
            W22 = W(end-Nx_patch_b-Ny_patch_b+2:end,end-Nx_patch_b-Ny_patch_b+2:end);
            
            bdy_D_patch = [bdy_D_n(x_patch_b)';fliplr(bdy_D_w(y_patch_b(1:end-1)))'];
            r = sqrt(radius_n^2 - bdy_D_patch'*(W22-W12'*(W11\W12))*bdy_D_patch);
            
            C = chol(W11);
            
            % generate random samples
            rho = rand(1,N_dict); rho = nthroot(rho,dim_r);
            
            % Generate unit vectors
            x = randn((Nx_patch_b-1)+(Ny_patch_b-1)-1,N_dict); 
            x_norm = sqrt(ones(1,(Nx_patch_b-1)+(Ny_patch_b-1)-1)*x.^2);
            x = C\x;
            x_unit = abs(x./x_norm);
        
            % Boundary conditions
            bdy_rand = r*rho.*x_unit;
        
            bdy_s = [bdy_D_w(y_patch_b(1))*ones(1,N_dict);...
                     bdy_rand(1:Nx_patch_b-1,:)];
            bdy_e = [bdy_rand(Nx_patch_b-1:end,:);...
                     bdy_D_n(x_patch_b(end))*ones(1,N_dict)];
            bdy_n = repmat(bdy_D_n(x_patch_b'),1,N_dict);
            bdy_w = repmat(bdy_D_w(y_patch_b'),1,N_dict);
            
        end
        
        
        %% Standard Solver  
        
        x_patch_o = x_nw_o(j):dx:x_se_o(j);
        y_patch_o = y_nw_o(k):dx:y_se_o(k);
        
        Nx_patch_o = length(x_patch_o);
        Ny_patch_o = length(y_patch_o);
        
        u = zeros(Ny_patch_o,Nx_patch_o,N_dict);
        phi = zeros(2*(Ny_patch_o-1)+2*(Nx_patch_o-1),N_dict);
        
        for i=1:N_dict
            
            u_temp = semilinear_elliptic_newton(x_patch_b,y_patch_b,dx,f,del_f,a,...
                                        bdy_w(:,i),bdy_e(:,i),bdy_s(:,i),bdy_n(:,i));
            
            ix1 = (x_nw_o(j)-x_nw_b(j))/dx;
            ix2 = (x_se_b(j)-x_se_o(j))/dx;
            iy1 = (y_nw_o(k)-y_nw_b(k))/dx;
            iy2 = (y_se_b(k)-y_se_o(k))/dx;
                                    
            u(:,:,i) = u_temp(iy1+1:end-iy2,ix1+1:end-ix2);
            
            phi(:,i) = [u_temp(iy1+1,ix1+1:end-ix2)';...
                        u_temp(iy1+2:end-iy2,end-ix2);...
                        u_temp(end-iy2,ix1+1:end-ix2-1)';...
                        u_temp(iy1+2:end-iy2-1,ix1+1)];
                            
        end
        
        
        %% compute weights W (order: South/East/North/West)
        bdy_x = [x_patch_o,x_patch_o(end)*ones(1,Ny_patch_o-1),...
                 x_patch_o(1:end-1),x_patch_o(1)*ones(1,Ny_patch_o-2)];
        
        bdy_y = [y_patch_o(1)*ones(1,Nx_patch_o),y_patch_o(2:end),...
                 y_patch_o(end)*ones(1,Nx_patch_o-1),y_patch_o(2:end-1)];
        
        disq =(bdy_x-bdy_x').^2+(bdy_y-bdy_y').^2;
        
        W = 2*dx^2./disq; W(isnan(W)|isinf(W)) = 0;
        W = diag( ones(1,2*(Nx_patch_o-1)+2*(Ny_patch_o-1))*W + dx ) - W;
        
        %% Save
        t_dic = toc(t_start)
        
        if j~=1 && j~=Mx && k~=1 && k~=My
            
            save(fullfile('data_elliptic',['G_H12_fno',int2str(f_no),...
                                  '_eps',int2str(n),...
                                  '_Mx',int2str(Mx),'_My',int2str(My),...
                                  '_(',int2str(j),',',int2str(k),')',...
                                  '_Nx',int2str(Nx),'_Ny',int2str(Ny),...
                                  '_r',num2str(radius_n,2),...
                                  '_Nsample',int2str(N_dict),...
                                  '_Nd',int2str(dim_r),...
                                  '_dxb',num2str(Dx_buffer),'_dxo',num2str(Dx_overlap),...
                                  '_expmt',int2str(expmt),'.mat']),'u','phi','W','t_dic');      
        else
            
            save(fullfile('data_elliptic',['G_H12_fno',int2str(f_no),...
                                  '_bdynoWESN',int2str(bdy_no_w),int2str(bdy_no_e),int2str(bdy_no_s),int2str(bdy_no_n),...
                                  '_eps',int2str(n),...
                                  '_Mx',int2str(Mx),'_My',int2str(My),...
                                  '_(',int2str(j),',',int2str(k),')',...
                                  '_Nx',int2str(Nx),'_Ny',int2str(Ny),...
                                  '_r',num2str(radius_n,2),...
                                  '_Nsample',int2str(N_dict),...
                                  '_Nd',int2str(dim_r),...
                                  '_dxb',num2str(Dx_buffer),'_dxo',num2str(Dx_overlap),...
                                  '_expmt',int2str(expmt),'.mat']),'u','phi','W','t_dic');
                              
        end
        
        
    end
end





