x = -1:0.01:1;
y = x;
[X,Y] = meshgrid(x,y);
surf(X,Y,real(soln(X,Y,1)))




function phi = soln(X,Y,t)
    
    rho = 2*pi*(0.6e9);
    r0 = 0.25;

    
    % X = X - 1.75;
    % Y = Y - 1.75;

    [R,Tht] = cart2pol(X,Y);

    Tht = zeros(size(Tht));

    phi = 0;

    for n = -100:1:100
        phi = phi - (1i)^abs(n).*besselj(abs(n),rho*r0).*exp(1i*n*Tht).*besselh(n,rho.*R)./besselh(n,rho.*r0);
    end
end



