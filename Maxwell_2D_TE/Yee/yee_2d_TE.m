
function  [delx,dely,delt,X_dual,Y_dual,H_z_new] = yee_2d_TE(k_mac)
             

mu0=1;   % Free space parameter
eps0= 1;
c = 1;


ND=20; delx=1/ND;     % Avoid dispersion   
Nx= round(1/delx);
Ny = Nx;
delx = 1/Nx;
dely = delx;

Final_T = 0.25;
delt= (0.9/c)*(1/sqrt((1/delx^2)+(1/dely^2))); 
NT = round((Final_T)/delt); 
delt = (Final_T)/NT;

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

%Initialize array
E_x_new = zeros(size(Y_main_X_dual));
E_y_new = zeros(size(Y_dual_X_main));

%Initial Data
E_x_old = zeros(size(Y_main_X_dual));
E_y_old = zeros(size(Y_dual_X_main));
H_z_old = zeros(size(X_dual));

%PEC Initial Data
% E_x_old = (1/sqrt(2))*cos(sqrt(2)*pi*(-0.5*delt)).*cos(pi*X_dual_Y_main).*sin(pi*Y_main_X_dual);
% E_y_old = -(1/sqrt(2))*cos(sqrt(2)*pi*(-0.5*delt)).*sin(pi*X_main_Y_dual).*cos(pi*Y_dual_X_main);
% H_z_old = zeros(size(X_dual));

%Material parameters
a = 1e-2;
mat_vol_frac = 0.0;
func = @(x) vol_frac(x,a) - mat_vol_frac;

if mat_vol_frac ~= 1
    r = fzero(func,[0,a/sqrt(2)]);
else
    r = a/sqrt(2);
end

eps_inf_Ex = eps0.*epsilon_inf(X_dual_Y_main(2:end-1,:),Y_main_X_dual(2:end-1,:),r,a,delx,dely);
eps_inf_Ey = eps0.*epsilon_inf(X_main_Y_dual(:,2:end-1),Y_dual_X_main(:,2:end-1),r,a,delx,dely);

%Source
source = b((delt/2):delt:(Final_T-delt/2),4);
src_loc = round(size(X_dual,1)*4/10);

    % vidfile = VideoWriter('hmm_vid.avi');
    % vidfile.FrameRate = 60;
    % open(vidfile);

for i = 1:NT
    if(mod(i,round(NT/100)) == 0)
        disp(round(i/round(NT/100)))

        %Uncomment for video
        % figure(1)
        % im = imagesc(x_dual,y_dual,H_z_old,[-5 5]);
        % im.AlphaData = 0.8;
        % hold on
        % im2 = imagesc(x_dual,y,eps_inf_Ex);
        % im2.AlphaData = 0.25;
        % hold off
        % pause(0.01)

        % max(H_z_old,[],'all')
    end


    %Yee's scheme at its essence
    E_x_new(2:end-1,:) = E_x_old(2:end-1,:) + (delt/dely)./eps_inf_Ex.*diff(H_z_old,1,1);
    E_y_new(:,2:end-1) = E_y_old(:,2:end-1) - (delt/delx)./eps_inf_Ey.*diff(H_z_old,1,2);
    H_z_new = H_z_old - delt./mu0.*(diff(E_y_new,1,2)./delx - diff(E_x_new,1,1)./dely);

    %Add source pulse
    H_z_new(src_loc:src_loc+1,src_loc:src_loc+1) = ...
        H_z_new(src_loc:src_loc+1,src_loc:src_loc+1)...
        + (delt)./(delx*dely)./mu0.*source(i)/2;
        

    %Update variables
    E_x_old = E_x_new;
    E_y_old = E_y_new;
    H_z_old = H_z_new;


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

function out = epsilon_inf(X,Y,r,a,delx,dely)
    %Material Parameters
    %Assumes periodic and circular inclusions
    if(any(size(X)~=size(Y)))
        error("Size of X and Y are not the same")
    end

    out = zeros(size(X));
    
    for i = [-1,1]
        for j = [-1,1]
            
            D = sqrt((mod(X+i*delx/4,a)-a/2).^2 + (mod(Y+j*dely/4,a)-a/2).^2);
            %D = max(abs(mod(X+i*delx/4,a)-a/2),abs(mod(Y+j*dely/4,a)-a/2));

            out(D<=r) = out(D<=r) + 1;
            out(D>r) = out(D>r) + 1;


            % out(X+i*delx/4<=(0.5)) = out(X+i*delx/4<=(0.5)) + 2.7;
            % out(X+i*delx/4>(0.5)) = out(X+i*delx/4>(0.5)) + 1.03;

        end
    end

    out = out/4;

    % out(X<=(0.5)) = 2.7;
    % out(X>(0.5)) = 1.03;
    % out(X==0.5) = mean([2.7,1.03]);
    
end

function out = vol_frac(r,a)
    
    out = length(r);

    in_ind = r<=a/2;
    out_ind = r>a/2;

    out(in_ind) = (pi*r(in_ind).^2)./(a^2);

    b = sqrt(r(out_ind).^2-(a.^2./4));
    out(out_ind) = (pi*r(out_ind).^2 - (4*r(out_ind).^2).*atan(2*b./a) + 2.*a.*b)./(a^2) ;

end

function out = b(t,f_bw)
    %Black-Harris Window Pluse
    an = [0.353222222,-0.488,0.145,-0.010222222];
    out = zeros(size(t));
    ind = t > 0 & t < 1.55/f_bw;

    for i = 1:length(an)
        %out(ind) = out(ind) + an(i).*cos(2*pi*(i-1)*(f_bw/1.55).*t(ind));
        out(ind) = out(ind) - (i-1).*an(i).*sin(2*pi*(i-1)*(f_bw/1.55).*t(ind));
    end

end