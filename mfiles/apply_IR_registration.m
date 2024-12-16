
function Image = apply_IR_registration(Image,para,mode)

%--------------------------------------------------------------------------
%   Image = apply_IR_registration(Image,para,mode)
%--------------------------------------------------------------------------
%   Deforms and reforms the input Image based on pre-calculated rigid
%   shifts and deformable registration indices
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

if contains(mode,'forward')
    Image = apply_rigid_shifts(Image,para,mode);
    Image = apply_mbr_shifts(Image,para,mode);
elseif contains(mode,'backward')
    Image = apply_mbr_shifts(Image,para,mode);
    Image = apply_rigid_shifts(Image,para,mode);
end

end
