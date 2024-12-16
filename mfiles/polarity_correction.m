function [Image,MBI] = polarity_correction(Image,MBI,para)

%--------------------------------------------------------------------------
%   [Image,MBI] = polarity_correction(Image,MBI,para)
%--------------------------------------------------------------------------
%   Function to performs signal polarity correction based on the
%   model-based images generated from the signal dictionary
%--------------------------------------------------------------------------
%   Inputs:      
%       - Image: full image series [nx,ny,nt,1,nsl]
%               - nx: number of measurements in spatial x-dimension
%               - ny: number of measurements in spatial y-dimension
%               - nt: number of time frames
%               - nsl: number of SMS slices (multi-band = 3)
%       - MBI: reference model-based images [nx,ny,nt,1,nsl]
%               - nx: number of measurements in spatial x-dimension
%               - ny: number of measurements in spatial y-dimension
%               - nt: number of time frames
%               - nsl: number of SMS slices (multi-band = 3)
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - Image: polarity corrected image series [nx,ny,nt,1,nsl]
%               - nx: number of measurements in spatial x-dimension
%               - ny: number of measurements in spatial y-dimension
%               - nt: number of time frames
%               - nsl: number of SMS slices (multi-band = 3)
%       - MBI: polarity corrected model-based images [nx,ny,nt,1,nsl]
%               - nx: number of measurements in spatial x-dimension
%               - ny: number of measurements in spatial y-dimension
%               - nt: number of time frames
%               - nsl: number of SMS slices (multi-band = 3)
%--------------------------------------------------------------------------


[nx,ny,nt,nset,nsl] = size(Image);

im1 = zeros(size(MBI));
im2 = zeros(size(Image));
for i=1:para.nof_block:nt

    block1 = MBI(:,:,i:i+para.nof_block-1,:,:);
    block2 = Image(:,:,i:i+para.nof_block-1,:,:);

    [~,idx] = min(block1,[],3);

    for j=1:size(idx,1)
        for k=1:size(idx,2)
            block1(j,k,1:idx(j,k),:,:) = -block1(j,k,1:idx(j,k),:,:);
            block2(j,k,1:idx(j,k),:,:) = -block2(j,k,1:idx(j,k),:,:);
        end
    end

    im1(:,:,i:i+para.nof_block-1,:,:) = block1;
    im2(:,:,i:i+para.nof_block-1,:,:) = block2;
end

MBI = im1;
Image = im2;

