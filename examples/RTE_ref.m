% This code generate the reference solution
clear;
% addpath(genpath('../src/'))
addpath('../src/RTE/')
%% Domain parameters
L = 3.0; Nx = 3*2^(14); dx = L/Nx; x0 = 0:dx:L; Nx = length(x0);

Nv = 2^10; [v0,w0] = legendre_quad(Nv-1); v0 = v0'; w0 = w0';
[xx,vv] = meshgrid(x0,v0); 
ploty = 4; plotz = 40;

%% Equation parameters
n = 4;
epsilon=2^(-n);

% v<0, x = L
bdy_vn = 2+sin(2*pi*v0(1:Nv/2)'); bdy_xL = 3; bdy_r_no = 1;    
% bdy_vn = 3-(v0(1:Nv/2).^2)'; bdy_xL = 1; bdy_r_no = 3;

% v>0, x=0
bdy_vp = 3+sin(2*pi*v0(Nv/2+1:end)'); bdy_x0 = 2; bdy_l_no = 1;  
% bdy_vp = 1.5-(v0(Nv/2+1:end).^2)'; bdy_x0 = 2.5; bdy_l_no = 3;

sigma = @(x) ones(size(x)); sigma_n = 0;
% sigma = @(x) 2+1.5*sin(2*pi*x); sigma_n = 1;


%% Compute ref solution %%%

theta_bdy = [bdy_x0;bdy_xL]; %[left;right]
f_bdy = [bdy_vn;bdy_vp]; %[v<0;v>0]

tic
sigma_x = sigma(x0)';
[f_ref,theta_ref] = nonlinearRTE_FD_aa_monotone(f_bdy,theta_bdy,epsilon,...
                                                        sigma_x,x0,w0,v0,dx); 
    
t_ref = toc

figure(2)

subplot(2,1,1); 
plot(x0,theta_ref); 
ylim([0,ploty]); 
title('Reference Solution');

subplot(2,1,2); 
mesh(xx,vv,f_ref); 
zlim([0,plotz]);


save(fullfile('data_RTE',['G_ref_eps',int2str(n),'_bdy',int2str(bdy_l_no),int2str(bdy_r_no),...
    '_sigma',int2str(sigma_n),'_Lx',num2str(L,2),'_dx',num2str(dx),'_Nv',int2str(Nv),...
    '.mat']),'f_ref','theta_ref','t_ref');

