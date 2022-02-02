% This code generate the dictionaries for interior patches
clear;
% addpath(genpath('../src/'))
addpath('../src/RTE/')
%% Domain decomposition
L = 3.0; Nx = 3*2^(11); dx = L/Nx; x0 = 0:dx:L; Nx = length(x0);
t = [0,0.25,0.75,1.25,1.75,2.25,2.75]; s = [0.25,0.75,1.25,1.75,2.25,2.75,3]; 
dx_overlap = 2^(-3); dx_buffer = 2^(-2);

t_b = [t(1),t(2:end)-dx_overlap-dx_buffer];
s_b = [s(1:end-1)+dx_overlap+dx_buffer,s(end)];

t = [t(1),t(2:end)-dx_overlap];
s = [s(1:end-1)+dx_overlap,s(end)];

Mx = length(t);

%% Equation and solver
n = 4;
epsilon=2^(-n);

Nv = 2^7; [v0,w0] = legendre_quad(Nv-1); v0 = v0'; w0 = w0';
w_bdy = [w0';1;1];

% sigma = @(x) 2+1.5*sin(2*pi*x); sigma_n = 1;
sigma = @(x) ones(size(x)); sigma_n = 0;

%% Sampling parameters
n_d = 2;          % dimension used in distribution of random sampling

N_sample = 64;

radius = 25;

expmt = 2;
    

%% Generate dictionaries
rng('shuffle');

if sigma_n == 0
    j_list = 2;
else
    j_list = 2:Mx-1;
end

tic
for j = j_list
    
    t_start = tic;
    %%     Generate Boundary Conditions     %%%
    % Generate radial distribution
    rho = rand(1,N_sample); rho = nthroot(rho,n_d);
    
    % Generate unit vectors
    x = randn(Nv+2,N_sample); x = x./sqrt(2*w_bdy);
    x_norm = sqrt(w_bdy'*x.^2);
    x_unit = abs(x./x_norm);
    
    % Boundary conditions
    boundary_buff = radius*rho.*x_unit;
    
    
    %%%%%%   assemble boundary conditions    %%%%%%
    
    f_boundary_buff = boundary_buff(1:Nv,:);
    theta_boundary_buff= boundary_buff(Nv+1:end,:);
   
    %% Generate local solutions
    x_patch_b = t_b(j):dx:s_b(j);
    sigma_x = sigma(x_patch_b)';
    
    %%% Full solutions
    Nx_patch = length(t(j):dx:s(j));
    
    f_dic = zeros(Nv,Nx_patch,N_sample);
    theta_dic = zeros(Nx_patch,N_sample);
    
    phi_dic = zeros(Nv+2,N_sample);
    
    for i=1:N_sample
        
        theta_bdy = theta_boundary_buff(:,i); %[left;right]
        f_bdy = f_boundary_buff(:,i); %[v<0;v>0]
        [f_temp,theta_temp] = nonlinearRTE_FD_aa_monotone(f_bdy,theta_bdy,epsilon,...
                                                 sigma_x,x_patch_b,w0,v0,dx);
        
        f_dic(:,:,i) = f_temp(:,dx_buffer/dx+1:end-dx_buffer/dx);
        theta_dic(:,i) = theta_temp(dx_buffer/dx+1:end-dx_buffer/dx);
        
        
        phi_dic(:,i) = [f_temp(1:Nv/2,end-dx_buffer/dx);
                       f_temp(Nv/2+1:end,dx_buffer/dx+1);
                       theta_temp(dx_buffer/dx+1);
                       theta_temp(end-dx_buffer/dx)];
        
        
    end
    
    
    t_dic = toc(t_start);
    if sigma_n == 0
        for k = 2:Mx-1
            save(fullfile('data_RTE',['G_eps',int2str(n),'_sigma',...
                int2str(sigma_n),'_Lx',num2str(L,2),'_dx',num2str(dx),'_Nv',...
                int2str(Nv),'_Mx',int2str(Mx),'(',int2str(k),')','_r',num2str(radius,2),...
                '_Nsample',int2str(N_sample),'_dx_overlap',num2str(dx_overlap),...
                '_dx_buffer',num2str(dx_buffer),'_expmt',int2str(expmt),'.mat']),...
                'f_dic','theta_dic','phi_dic','t_dic');
        end
    else
        save(fullfile('data_RTE',['G_eps',int2str(n),'_sigma',int2str(sigma_n),...
            '_Lx',num2str(L,2),'_dx',num2str(dx),'_Nv',int2str(Nv),'_Mx',int2str(Mx),...
            '(',int2str(j),')','_r',num2str(radius,2),'_Nsample',int2str(N_sample),...
            '_dx_overlap',num2str(dx_overlap),'_dx_buffer',num2str(dx_buffer),...
            '_expmt',int2str(expmt),'.mat']),'f_dic','theta_dic','phi_dic','t_dic');
    end

end

toc



