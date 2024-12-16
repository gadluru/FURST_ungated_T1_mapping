function [motion_mask,cardiac_mask] = get_cardiac_mask(Image,para)

%--------------------------------------------------------------------------
%   [motion_mask,cardiac_mask] = get_cardiac_mask(Image,para)
%--------------------------------------------------------------------------
%   Function to get motion mask and cardiac mask
%--------------------------------------------------------------------------
%   Inputs:      
%       - Image: full image series [nx,ny,nt,1,nsl]
%               - nx: number of measurements in spatial x-dimension
%               - ny: number of measurements in spatial y-dimension
%               - nt: number of time frames
%               - nsl: number of SMS slices (multi-band = 3)
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Output:
%       - motion_mask: mask for motion compensation [nx,ny,1,1,nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nt: number of time frames
%           - nsl: number of SMS slices (multi-band = 3)
%       - cardiac_mask: mask for self-gating [nx,ny,1,1,nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nt: number of time frames
%           - nsl: number of SMS slices (multi-band = 3)
%--------------------------------------------------------------------------


[nx,ny,nt,nset,nsl] = size(Image);

Image = reshape(Image,[nx,ny,nt,nset*nsl]);

image_for_std = Image(:,:,1:para.nof_block,:);
std_map = std(image_for_std,0,3);
std_map = imgaussfilt(std_map,5);
filter = fspecial('gaussian',nx,nx/10);
std_map = std_map.*filter;
mask = zeros(nx,nx,1,nset*nsl);
for i=1:nset*nsl
    [x,y] = find(std_map(:,:,1,i)==max(max(std_map(:,:,1,i))));
    mask(x,y,1,i) = 1;
    mask(:,:,1,i) = bwdist(mask(:,:,1,i));
end
mask = reshape(mask,[nx,ny,1,nset,nsl]);

cardiac_mask = mask < 25;
motion_mask = mask < 35;
