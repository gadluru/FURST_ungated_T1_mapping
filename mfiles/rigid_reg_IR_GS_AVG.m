function [registered,shifts] = rigid_reg_IR_GS_AVG(Image,para)

%--------------------------------------------------------------------------
%   [registered,para] = rigid_reg_IR_GS_AVG(Image,para)
%--------------------------------------------------------------------------
%   Function to perform rigid registration between images acquired at
%   similar inversions times betwen different inversion groups
%--------------------------------------------------------------------------
%   Inputs:      
%       - Image: image series to be rigid registered [nx,ny,nt,1,nsl]
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

[nx,ny,nt,nset,nsl] = size(Image);

Image = reshape(Image,[nx,ny,nt*nset*nsl]);
parfor i=1:nt*nset*nsl
    Image(:,:,i) = medfilt2(Image(:,:,i));
end
Image = reshape(Image,[nx,ny,nt,nset,nsl]);

[motion_mask,~] = get_cardiac_mask(Image,para); motion_mask = any(any(motion_mask,4),5);

Image = mean(reshape(Image,[nx,ny,para.nof_block,nt/para.nof_block,nset,nsl]),3);

[nx,ny,nb,ni,nset,nsl] = size(Image);

Image = reshape(permute(Image,[1,2,4,3,5,6]),[nx,ny,ni,nb*nset*nsl]); registered = Image;

noi = 4; shifts = zeros(2,noi,ni,nb*nset*nsl);
for i=1:nb*nset*nsl
    step = 100;
    for j=1:noi
        for k=ni-1:-1:1
            [registered(:,:,k,i),shifts(:,j,k,i)] = Reg_GS_rigid(Image(:,:,k,i),mean(registered(:,:,k+1:end,i),3),step,10,motion_mask);
        end
        Image = registered; step = step/2;
    end
end

shifts = permute(reshape(shifts,[2,noi,ni,nb,nset,nsl]),[1,2,4,3,5,6]);

shifts = single(reshape(repmat(mean(mean(shifts,5),6),[1,1,nt/ni,1,nset,nsl]),[2,noi,nt,nset,nsl]));

end

