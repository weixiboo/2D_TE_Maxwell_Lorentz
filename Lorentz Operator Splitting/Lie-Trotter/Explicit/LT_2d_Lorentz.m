
function  [delx,dely,delt,X_dual,Y_dual,H_z_new] = LT_2d_Lorentz(k_mac)
             

%[X_dual_Y_main,Y_main_X_dual,E_x_new]
%[X_main_Y_dual,Y_dual_X_main,E_y_new]
%[X_dual,Y_dual,H_z_new]

final_T =  0.5;

mu0=1;       % Free space permeability
eps0 = 1;
c = 1;
[eps_inf,nu,omega_0,omega_p,sigma] = lorentz_parameters(1);
c_inf = 1/sqrt(eps0*eps_inf*mu0);

ND=20; delx=1/ND;     % Avoid dispersion   
Nx= round(1/delx);
Ny = Nx;
delx = 1/Nx;
dely = delx;

delt= (0.99/c)*(1/sqrt((1/delx^2)+(1/dely^2))); 
NT = round((final_T)/delt); 
delt = (final_T)/NT;

delx = delx/k_mac;
dely = delx;
delt = delt/k_mac;
NT = NT*k_mac;
Nx = Nx*k_mac;
Ny = Ny*k_mac;

%X's and T's
x = 0:delx:Nx*delx;
y = 0:dely:Ny*dely;
x_dual = x(2:end) - delx/2;
y_dual = y(2:end) - dely/2;

[X_dual,Y_dual] = meshgrid(x_dual,y_dual);
[X_dual_Y_main,Y_main_X_dual] = meshgrid(x_dual,y);
[X_main_Y_dual,Y_dual_X_main] = meshgrid(x,y_dual);

%Initial Data
% E_x_old = zeros(size(Y_main_X_dual));
% E_y_old = zeros(size(Y_dual_X_main));
% H_z_old = zeros(size(X_dual));
% 
% J_x_old = zeros(size(Y_main_X_dual));
% J_y_old = zeros(size(Y_dual_X_main));
% P_x_old = zeros(size(Y_main_X_dual));
% P_y_old = zeros(size(Y_dual_X_main));

%PEC Initial Data
wave_k = sqrt(2*pi^2);
theta = give_theta(wave_k,nu);
E_x_old =  -theta.*cos(pi*X_dual_Y_main).*sin(pi*Y_main_X_dual);
E_y_old =  theta.*sin(pi*X_main_Y_dual).*cos(pi*Y_dual_X_main);
H_z_old = ((wave_k^2)./pi).*cos(pi*X_dual).*cos(pi*Y_dual);

J_x_old =  -(theta^2+wave_k^2).*cos(pi*X_dual_Y_main).*sin(pi*Y_main_X_dual);
J_y_old = (theta^2+wave_k^2).*sin(pi*X_main_Y_dual).*cos(pi*Y_dual_X_main);

P_x_old = -(-theta - wave_k^2/theta).*cos(pi*X_dual_Y_main).*sin(pi*Y_main_X_dual);
P_y_old = (-theta - wave_k^2/theta).*sin(pi*X_main_Y_dual).*cos(pi*Y_dual_X_main);


%KF scheme parameters
LT_J1 = 1/delt + nu + (delt.*omega_0.^2)./4 + delt.*omega_p.^2./(4.*eps_inf);
LT_J2 = 1/delt - nu - (delt.*omega_0.^2)./4 - delt.*omega_p.^2./(4.*eps_inf);
LT_J = LT_J2./LT_J1;


%     vidfile = VideoWriter('hmm_vid.avi');
%     vidfile.FrameRate = 60;
%     open(vidfile);

for i = 1:NT
    if(mod(i,round(NT/100)) == 0)
        disp(round(i/round(NT/100)))

        % figure(1)
        % imagesc(x_dual,y_dual,H_z_old,[-1 1])
        % pause(0.01)
    end

    %Step A: 
    H_z_new = H_z_old - (delt/2)./mu0.*(diff(E_y_old,1,2)./delx - diff(E_x_old,1,1)./dely);
        
    E_x_new = E_x_old + ((delt)./(dely*eps0.*eps_inf)).*matrix_padding(diff(H_z_new,1,1),1);
    E_y_new = E_y_old - ((delt)./(delx*eps0.*eps_inf)).*matrix_padding(diff(H_z_new,1,2),2);

    H_z_new = H_z_new - (delt/2)./mu0.*(diff(E_y_new,1,2)./delx - diff(E_x_new,1,1)./dely);
    


    %Step B: Solve for new E and J's
    J_x_new = LT_J.*J_x_old - (omega_0.^2./LT_J1).*P_x_old + eps0*omega_p.^2./LT_J1.*E_x_new;
    J_y_new = LT_J.*J_y_old - (omega_0.^2./LT_J1).*P_y_old + eps0*omega_p.^2./LT_J1.*E_y_new;

    P_x_new = P_x_old + (delt/2)*(J_x_new+J_x_old);
    P_y_new = P_y_old + (delt/2)*(J_y_new+J_y_old);

    E_x_new = E_x_new - delt./(2*eps0*eps_inf).*(J_x_new + J_x_old);
    E_y_new = E_y_new - delt./(2*eps0*eps_inf).*(J_y_new + J_y_old);


    %Update variables
    E_x_old = E_x_new;
    E_y_old = E_y_new;
    H_z_old = H_z_new;

    J_x_old = J_x_new;
    J_y_old = J_y_new;

    P_x_old = P_x_new;
    P_y_old = P_y_new;




%         fig = figure('visible','off');
%         plot(x_dual,E_new)
%         xlabel('x')
%         ylabel('Electric Field')
%         title('Numerical Solution with HMM')
%         axis([0 Nx*delx -0.2 1])
%         pause(1/1000)
%         writeVideo(vidfile,getframe(gcf));

end

%          close(vidfile)

end

function out = maxwell_step_matrix(eps,mu,delw,delt,N)
    c_inf = 1/sqrt(eps*mu);
    eta = c_inf*delt/delw;

    out = eye(N);
    out = out + diag(ones(N-1,1).*(-eta^2/4),-1) ...
        + diag(ones(N-1,1).*(-eta^2/4),1) ...
        + diag(ones(N,1).*(eta^2/2));
    
    
    out(1,:) = 0;
    out(1,1) = 1;
    out(end,:) = 0;
    out(end,end) = 1;

end


function out = matrix_padding(in,dim)
    
    if dim == 1
        pad_len = size(in,2);
        out = [zeros(1,pad_len); in; zeros(1,pad_len)];        
    elseif dim == 2
        pad_len = size(in,1);
        out = [zeros(pad_len,1) in zeros(pad_len,1)];   
    end

end

function [eps_inf,nu,omega_0,omega_p,sigma] = lorentz_parameters(x)
    eps_inf = 1;
    eps_s = 2;
    tau = 0.4;
    omega_0 = 1;
    sigma = 0;

    nu = 1/(2*tau);
    omega_p = omega_0*sqrt(eps_s - eps_inf);
    
end

function theta = give_theta(k,nu)
    f = @(x) x.^4 - 2*nu*x.^3 + (2+k^2).*x.^2 - 2*nu*k^2.*x + k^2;
    theta = fzero(f,[0,1]);
end
