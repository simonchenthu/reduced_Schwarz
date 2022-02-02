% This code computes the reference solution to the elliptic equation using
% finite difference solver and Newton iteration
clear;
% addpath(genpath('../src/'))
addpath('../src/elliptic/')
%% Domain parameters
Lx = 1.0; Nx = 2^(12); dx = Lx/Nx; x0 = 0:dx:Lx; Nx = length(x0); 
Ly = 1.0; Ny = 2^(12); dy = Ly/Ny; y0 = 0:dy:Ly; Ny = length(y0);
plotz = 2; 
expmt = 1;
[xx,yy] = meshgrid(x0,y0);

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

bdy_D_n = @(x) sin(2*pi*x); bdy_no_n = 1;
bdy_D_s = @(x) -sin(2*pi*x); bdy_no_s = 1;
bdy_D_w = @(y) sin(2*pi*y); bdy_no_w = 1;
bdy_D_e = @(y) -sin(2*pi*y);  bdy_no_e = 1;

%% Finite difference solver
bdy_w = bdy_D_w(y0)';
bdy_e = bdy_D_e(y0)';
bdy_n = bdy_D_n(x0)';
bdy_s = bdy_D_s(x0)';

tic

u_ref = semilinear_elliptic_newton(x0,y0,dx,f,del_f,a,...
                                   bdy_w,bdy_e,bdy_s,bdy_n);
                         
t_ref = toc

%% Save
save(fullfile('data_elliptic',['u_ref_fno',int2str(f_no),'_eps',int2str(n),...
    '_bdy',int2str(bdy_no_n),int2str(bdy_no_s),int2str(bdy_no_w),int2str(bdy_no_e),...
    '_a',int2str(a_n),'_dx',num2str(dx),'.mat']),'u_ref','t_ref');
                             
%% Plot
figure(1)
mesh(xx,yy,u_ref); 
xlim([0,Lx]); ylim([0,Ly]); zlim([-plotz,plotz]);


