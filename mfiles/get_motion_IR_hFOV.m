function [registered,idx] = get_motion_IR_hFOV(Image,MBI,para,step,noi)

%--------------------------------------------------------------------------
%   [registered,idx] = get_motion_IR_hFOV(Image,MBI,para,step,noi)
%--------------------------------------------------------------------------
%   Function to perform final deformation registration based on model-based
%   images
%--------------------------------------------------------------------------
%   Inputs:      
%       - Image: image series to be registered [nx,ny,nt,1,nsl]
%          - nx: number of measurements in spatial x-dimension
%          - ny: number of measurements in spatial y-dimension
%          - nt: number of time frames
%          - nsl: number of SMS slices (multi-band = 3)
%       - MBI: reference model-based images to be registered [nx,ny,nt,1,nsl]
%          - nx: number of measurements in spatial x-dimension
%          - ny: number of measurements in spatial y-dimension
%          - nt: number of time frames
%          - nsl: number of SMS slices (multi-band = 3)
%       - para: reconstruction parameters [structure]
%       - step: step size used for each registration iteration
%       - noi: number of iterations for deformable registration [scalar]
%--------------------------------------------------------------------------
%   Outputs:
%       - registered: final registered images [nx,ny,nt,1,nsl]
%          - nx: spatial x-dimension
%          - ny: spatial y-dimension
%          - nt: number of time frames
%          - nsl: number of slices (SMS, multiband=3)
%       - idx: deformable registration indices for final registration [nx,ny,nt,1,nsl]
%          - nx: spatial x-dimension
%          - ny: spatial y-dimension
%          - nt: number of time frames
%          - nsl: number of slices (SMS, multiband=3)
%--------------------------------------------------------------------------

[nx,ny,nt,nset,nsl] = size(Image);

registered = abs(Image); Image = crop_half_FOV(abs(Image)); MBI = crop_half_FOV(abs(MBI));

Image = reshape(Image,[nx/2,ny/2,nt*nset*nsl]);
parfor i=1:nt*nset*nsl
    Image(:,:,i) = medfilt2(Image(:,:,i));
end
Image = reshape(Image,[nx/2,ny/2,nt,nset,nsl]);

[motion_mask,~] = get_cardiac_mask(Image,para); motion_mask = any(any(motion_mask,4),5);

[Y,X,T,N,S] = ndgrid(1:nx,1:ny,1:nt,1:nset,1:nsl);

x = X; y = Y;

x = reshape(x,[nx,ny,nt*nset*nsl]);
y = reshape(y,[nx,ny,nt*nset*nsl]);

dx = (nx - nx/2)/2;
dy = (ny - ny/2)/2;

center_x = dx + 1:dx + nx/2;
center_y = dy + 1:dy + ny/2;

registered = reshape(registered,[nx,ny,nt*nset*nsl]); Image = reshape(Image,[nx/2,ny/2,nt*nset*nsl]); MBI = reshape(MBI,[nx/2,ny/2,nt*nset*nsl]);

parfor i=1:nt*nset*nsl
    [registered(center_x,center_y,i),y(center_x,center_y,i),x(center_x,center_y,i)] = Reg_GS_tv(Image(:,:,i),MBI(:,:,i),motion_mask,step,noi);
end

registered = reshape(registered,[nx,ny,nt,nset,nsl]);

x = reshape(x,[nx,ny,nt,nset,nsl]);
y = reshape(y,[nx,ny,nt,nset,nsl]);

x(center_x,center_y,:,:,:) = x(center_x,center_y,:,:,:) + dx;
y(center_x,center_y,:,:,:) = y(center_x,center_y,:,:,:) + dy;

idx = int32(sub2ind([nx,ny,nt,nset,nsl],round(y),round(x),T,N,S));

