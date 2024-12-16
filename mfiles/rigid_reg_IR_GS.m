function [registered,shifts] = rigid_reg_IR_GS(Image,MBI,para)

%--------------------------------------------------------------------------
%   [registered,shifts] = rigid_reg_IR_GS(Image,MBI,para)
%--------------------------------------------------------------------------
%   Function to perform rigid registration between T1-weighted images and
%   their corresponding model-based images
%--------------------------------------------------------------------------
%   Inputs:      
%       - Image: image series to be rigid registered [nx,ny,nt,1,nsl]
%          - nx: number of measurements in spatial x-dimension
%          - ny: number of measurements in spatial y-dimension
%          - nt: number of time frames
%          - nsl: number of SMS slices (multi-band = 3)
%       - Image: reference model-based images for registration [nx,ny,nt,1,nsl]
%          - nx: number of measurements in spatial x-dimension
%          - ny: number of measurements in spatial y-dimension
%          - nt: number of time frames
%          - nsl: number of SMS slices (multi-band = 3)
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - registered: rigid registered images [nx,ny,nt,1,nsl]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of time frames
%           - nsl: number of slices (SMS, multiband=3)
%       - shifts: calculated x- and y- rigid shifts [2,noi,nt,1,nsl]
%           - noi: number of iterations
%           - nt: number of time frames
%           - nsl: number of slices (SMS, multiband=3)
%--------------------------------------------------------------------------

Image = crop_half_FOV(abs(Image));
MBI = crop_half_FOV(abs(MBI));

[nx,ny,nt,nset,nsl] = size(Image);

Image = reshape(Image,[nx,ny,nt*nset*nsl]);
parfor i=1:nt*nset*nsl
    Image(:,:,i) = medfilt2(Image(:,:,i));
end
Image = reshape(Image,[nx,ny,nt,nset,nsl]);

[motion_mask,~] = get_cardiac_mask(Image,para); motion_mask = any(any(motion_mask,4),5);

Image = reshape(Image,[nx,ny,nt*nset*nsl]); MBI = reshape(MBI,[nx,ny,nt*nset*nsl]); registered = Image;

noi = 4; step = 100; shifts = zeros(2,noi,nt*nset*nsl);
for i=1:noi
    parfor j=1:nt*nset*nsl
        [registered(:,:,j),shifts(:,i,j)] = Reg_GS_rigid(Image(:,:,j),MBI(:,:,j),step,10,motion_mask);
    end
    Image = registered; step = step/2;
end

shifts = single(reshape(shifts,[2,noi,nt,nset,nsl]));

end

