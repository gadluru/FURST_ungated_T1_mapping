function [I2,x,y] = Reg_GS_tv(I1,I0,mask,step,noi)

%--------------------------------------------------------------------------
%   [I2,x,y] = Reg_GS_tv(I1,I0,mask,step,noi)
%--------------------------------------------------------------------------
%   Function to perform image-to-image deformable registration
%--------------------------------------------------------------------------
%   Inputs:      
%       - I1: image to be registered [nx,ny]
%               - nx: number of measurements in spatial x-dimension
%               - ny: number of measurements in spatial y-dimension
%       - I0: reference image for registration [nx,ny]
%       - mask: registration mask for region of interest
%               - nx: number of measurements in spatial x-dimension
%               - ny: number of measurements in spatial y-dimension
%       - step: step size used for each registration update [scalar]
%       - noi: number of iterations used for registration [scalar]
%--------------------------------------------------------------------------
%   Outputs:
%       - I2: rigid registered images [nx,ny]
%          - nx: spatial x-dimension
%          - ny: spatial y-dimension
%       - x: x-deformations used to produce I2 [nx,ny]
%          - nx: spatial x-dimension
%          - ny: spatial y-dimension
%       - y: y-deformations used to produce I2 [nx,ny]
%          - nx: spatial x-dimension
%          - ny: spatial y-dimension
%--------------------------------------------------------------------------

[sx,sy] = size(I0);
Maxiter = noi;
a = 30;
s = 3;

kx = cos(2*pi*(0:sx/(sx-1):sx));
ky = cos(2*pi*(0:sy/(sy-1):sy));
W = 2*(kx+ky.'-2);
W = single((1-a*W).^-s);
 
[y,x] = meshgrid(1:sx,1:sy);

for iter=1:Maxiter
    
    I2 = interp2(I1,y,x);
    
    ddx = 0.5*(I2(3:end,:) - I2(1:end-2,:));
    ddy = 0.5*(I2(:,3:end) - I2(:,1:end-2));
    ddx = cat(1,I2(2,:) - I2(1,:),ddx,I2(end,:) - I2(end-1,:));
    ddy = cat(2,I2(:,2) - I2(:,1),ddy,I2(:,end) - I2(:,end-1));
    
    dI = I2-I0;
    dI = dI - medfilt2(dI,[3,3]);

    J = jacobian(y,x);

    dx = W.*fft2(ddx.*dI.*J);
    dx = -real(ifft2(dx));

    dy = W.*fft2(ddy.*dI.*J);
    dy = -real(ifft2(dy));

    d = sqrt(dx.^2+dy.^2);

    md = max(d(:));
 
    dx = dx/md;
    dy = dy/md;
    
    E(iter) = sum(vec(d .* mask));
    if iter>1 && E(iter)>E(iter-1)
        step = step/2;
    end

    x = x + step*dx;
    y = y + step*dy;
    
    x(x<1) = 1;
    x(x>sx) = sx;
    y(y<1) = 1;
    y(y>sy) = sy;

end

end

function J = jacobian(x,y)
    [fxx,fxy] = gradient(x);
    [fyx,fyy] = gradient(y);
    J = fxx .* fyy - fxy .* fyx;
end
