
mu0=1;   % Free space permeability
eps0= 1;
c = 1;


ND=20; delx=1/ND;     % Avoid dispersion   
Nx= round(1/delx);
Ny = Nx;
delx = 1/Nx;
dely = delx;

Final_T = 0.75;
delt= (0.99/c)*(1/sqrt((1/delx^2)+(1/dely^2))); 
NT = round((Final_T)/delt); 
delt = (Final_T)/NT;

%Source
my_source = b((delt/2):delt:(Final_T-delt/2),6);%.*sin(2*pi*1.*((delt/2):delt:(Final_T-delt/2)));

y = fft(my_source);
fs = 1/delt;
f = (0:length(y)-1)*fs/length(y);
plot(f(1:round(length(f)/2)),abs(y(1:round(length(f)/2))))
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Magnitude')



function out = b(t,f_bw)

    an = [0.353222222,-0.488,0.145,-0.010222222];
    out = zeros(size(t));
    ind = t > 0 & t < 1.55/f_bw;

    for i = 1:length(an)
        %out(ind) = out(ind) + an(i).*cos(2*pi*(i-1)*(f_bw/1.55).*t(ind));
        out(ind) = out(ind) - (i-1).*an(i).*sin(2*pi*(i-1)*(f_bw/1.55).*t(ind));
    end

end