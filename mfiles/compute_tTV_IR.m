function update_ttv = compute_tTV_IR(Image,weight,epsilon,para)

%--------------------------------------------------------------------------
%   update_ttv = compute_tTV_IR(Image,weight,epsilon,para)
%--------------------------------------------------------------------------
%   computes the temporal total variation update term for time frames within
%   the same inversion group.
%--------------------------------------------------------------------------
%   Inputs:
%       - image: image series to be reconstructed [nx,ny,nt,1,nsl]
%       - weight: temporal total variation regularization weight [scalar]
%       - epsilon: small value to prevent singularity [scalar]
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - update_ttv: temporal total vaiaration update term [nx,ny,nt,1,nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nt: number of time frames
%           - nsl: number of SMS slices (multi-band = 3)
%--------------------------------------------------------------------------

if weight ~= 0

    [sx,sy,nt,nset,nsl] = size(Image);

    Image = reshape(Image,[sx,sy,para.nof_block,nt/para.nof_block,nset,nsl]);

    temp_a = diff(Image,1,3);
    temp_b = temp_a./(sqrt((abs(temp_a).^2) + epsilon)); clear temp_a
    temp_c = diff(temp_b,1,3);
    update_ttv = cat(3,temp_b(:,:,1,:,:,:),temp_c,-temp_b(:,:,end,:,:,:)); clear temp_b temp_c
    update_ttv = weight .* reshape(update_ttv,[sx,sy,nt,nset,nsl]);
    
else
    update_ttv = 0;
end

end