% This code generate the dictionaries for boundary patches
clear;
% addpath(genpath('../src/'))
addpath('../src/RTE/')
%% Domain decomposition
L = 3.0; Nx = 3*2^(11); dx = L/Nx; x0 = 0:dx:L; Nx = length(x0);
t = [0,0.25,0.75,1.25,1.75,2.25,2.75]; s = [0.25,0.75,1.25,1.75,2.25,2.75,3]; 
dx_overlap = 2^(-3); dx_buffer = 2^(-2);

t_b = [t(1),t(2:end)-dx_overlap-dx_buffer];         % left end points of buffered patches
s_b = [s(1:end-1)+dx_overlap+dx_buffer,s(end)];     % right end points of buffered patches

t = [t(1),t(2:end)-dx_overlap];          % left end points of patches
s = [s(1:end-1)+dx_overlap,s(end)];      % right end points of patches

Mx = length(t);

%% Equation and solver
n = 4;
epsilon=2^(-n);

Nv = 2^7; [v0,w0] = legendre_quad(Nv-1); v0 = v0'; w0 = w0'; 
w_l = [w0(Nv/2+1:Nv)';1];
w_r = [w0(1:Nv/2)';1];

% sigma = @(x) 2+1.5*sin(2*pi*x); sigma_n = 1;
sigma = @(x) ones(size(x)); sigma_n = 0;

%%% v<0, x = L %%%
bdy_vn = 2+sin(2*pi*v0(1:Nv/2)'); bdy_xL = 3; bdy_r_no = 1;    
% bdy_vn = 2-(v0(1:Nv/2).^2)'; bdy_xL = 3; bdy_r_no = 2;
% bdy_vn = 3-(v0(1:Nv/2).^2)'; bdy_xL = 1; bdy_r_no = 3;

%%% v>0, x=0 %%%
bdy_vp = 3+sin(2*pi*v0(Nv/2+1:end)'); bdy_x0 = 2; bdy_l_no = 1;  
% bdy_vp = 3-(v0(Nv/2+1:end).^2)'; bdy_x0 = 1; bdy_l_no = 2;
% bdy_vp = 1.5-(v0(Nv/2+1:end).^2)'; bdy_x0 = 2.5; bdy_l_no = 3;
                            

%% Sampling parameters
n_d = 2;          % used in distribution of random sampling

N_sample = 64;

radius = 25;

expmt = 2;
    

%%     Generate Boundary Conditions     %%%
rng('shuffle');

tic;

radius_l = sqrt(radius^2-w_l'*([bdy_vp;bdy_x0].^2));
radius_r = sqrt(radius^2-w_r'*([bdy_vn;bdy_xL].^2));

% Generate radial distribution
rho_l = rand(1,N_sample); rho_l = nthroot(rho_l,n_d);
rho_r = rand(1,N_sample); rho_r = nthroot(rho_r,n_d);

% Generate unit vectors
x_l = randn(Nv/2+1,N_sample); x_l = x_l./sqrt(2*w_r);
x_l_norm = sqrt(w_r'*x_l.^2);
x_l_unit = abs(x_l./x_l_norm);
    
x_r = randn(Nv/2+1,N_sample); x_r = x_r./sqrt(2*w_l);
x_r_norm = sqrt(w_l'*x_r.^2);
x_r_unit = abs(x_r./x_r_norm);

% Boundary conditions
boundary_rand_l = radius_l*rho_l.*x_l_unit;
boundary_rand_r = radius_r*rho_r.*x_r_unit;


%%%%%%   assemble boundary conditions    %%%%%%

f_boundary_buff_l = [boundary_rand_l(1:Nv/2,:);bdy_vp*ones(1,N_sample)];
theta_boundary_buff_l= [bdy_x0*ones(1,N_sample);boundary_rand_l(end,:)];
                  
f_boundary_buff_r = [bdy_vn*ones(1,N_sample);boundary_rand_r(1:Nv/2,:)];
theta_boundary_buff_r= [boundary_rand_r(end,:);bdy_xL*ones(1,N_sample)];


x_patch_b_l = t_b(1):dx:s_b(1);
x_patch_b_r = t_b(end):dx:s_b(end);

sigma_x_l = sigma(x_patch_b_l)';
sigma_x_r = sigma(x_patch_b_r)';


Nx_patch_l = length(t(1):dx:s(1));
Nx_patch_r = length(t(end):dx:s(end));


%% solve dictionary on left patch
time_l = tic;

f_dic = zeros(Nv,Nx_patch_l,N_sample);
theta_dic = zeros(Nx_patch_l,N_sample);
for i=1:N_sample
    
    theta_bdy = theta_boundary_buff_l(:,i); %[left;right]
    f_bdy = f_boundary_buff_l(:,i); %[v<0;v>0]
    [f_temp,theta_temp] = nonlinearRTE_FD_aa_monotone(f_bdy,theta_bdy,epsilon,...
                                               sigma_x_l,x_patch_b_l,w0,v0,dx);
    
    f_dic(:,:,i) = f_temp(:,1:end-dx_buffer/dx);
    theta_dic(:,i) = theta_temp(1:end-dx_buffer/dx);
    
    phi_dic(:,i) = [f_temp(1:Nv/2,end-dx_buffer/dx);
                       bdy_vp;
                       bdy_x0;
                       theta_temp(end-dx_buffer/dx)];
    
end
t_dic = toc(time_l)

save(fullfile('data_RTE',['G_bdy',int2str(bdy_l_no),'_eps',int2str(n),...
    '_sigma',int2str(sigma_n),'_Lx',num2str(L,2),'_dx',num2str(dx),'_Nv',int2str(Nv),...
    '_Mx',int2str(Mx),'(',int2str(1),')','_r',num2str(radius,2),'_Nsample',int2str(N_sample),...
    '_dx_overlap',num2str(dx_overlap),'_dx_buffer',num2str(dx_buffer),'_expmt',int2str(expmt),...
    '.mat']),'f_dic','theta_dic','phi_dic','t_dic');
    


%% solve dictionary on right patch
time_r = tic;

f_dic = zeros(Nv,Nx_patch_r,N_sample);
theta_dic = zeros(Nx_patch_r,N_sample);
for i=1:N_sample
    
    theta_bdy = theta_boundary_buff_r(:,i); %[left;right]
    f_bdy = f_boundary_buff_r(:,i); %[v<0;v>0]
    [f_temp,theta_temp] = nonlinearRTE_FD_aa_monotone(f_bdy,theta_bdy,epsilon,...
                                               sigma_x_l,x_patch_b_r,w0,v0,dx); 
    
    f_dic(:,:,i) = f_temp(:,dx_buffer/dx+1:end);
    theta_dic(:,i) = theta_temp(dx_buffer/dx+1:end);
    
    phi_dic(:,i) = [bdy_vn;
                       f_temp(Nv/2+1:end,dx_buffer/dx+1);
                       theta_temp(dx_buffer/dx+1);
                       bdy_xL];
    
end

t_dic = toc(time_r)

save(fullfile('data_RTE',['G_bdy',int2str(bdy_r_no),'_eps',int2str(n),...
    '_sigma',int2str(sigma_n),'_Lx',num2str(L,2),'_dx',num2str(dx),'_Nv',int2str(Nv),...
    '_Mx',int2str(Mx),'(',int2str(Mx),')','_r',num2str(radius,2),'_Nsample',int2str(N_sample),...
    '_dx_overlap',num2str(dx_overlap),'_dx_buffer',num2str(dx_buffer),'_expmt',int2str(expmt),...
    '.mat']),'f_dic','theta_dic','phi_dic','t_dic');

toc


