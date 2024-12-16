
function Image = apply_mbr_shifts(Image,para,mode)

%--------------------------------------------------------------------------
%   Image = apply_mbr_shifts(Image,para,mode)
%--------------------------------------------------------------------------
%   Deforms and reforms the input Image based on the pre-calculated
%   deformable registration indices
%--------------------------------------------------------------------------
%   Inputs:      
%       - Image: motion compensated/non-motion compensated images [nx,ny,nt,1,nsl]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of time frames
%           - nsl: number of slices (SMS, multiband=3)
%       - para: reconstruction parameters [structure]
%       - mode: deforms/reforms images ['forward','backward']
%--------------------------------------------------------------------------
%   Outputs:
%       - Image: motion compensated/non-motion compensated images [nx,ny,nt,1,nsl]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of time frames
%           - nsl: number of slices (SMS, multiband=3)
%--------------------------------------------------------------------------

if para.setting.ifGPU
    Image = gather(Image);
end

shifts = para.Motion.MBR.idx;

if contains(mode,'forward')
    for i=1:size(shifts,6)
        Image = Image(shifts(:,:,:,:,:,i));
    end
elseif contains(mode,'backward')
    for i=size(shifts,6):-1:1
        Image(shifts(:,:,:,:,:,i)) = Image;
    end
end

if para.setting.ifGPU
    Image = gpuArray(Image);
end

end