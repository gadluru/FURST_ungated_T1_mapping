function tTV_update = compute_tTV_IR_reorder(Image,Data,weight,epsilon,para)

%--------------------------------------------------------------------------
%   tTV_update = compute_tTV_IR_reorder(Image,Data,weight,epsilon,para)
%--------------------------------------------------------------------------
%   computes the reordered temporal total variation update term 
%--------------------------------------------------------------------------
%   Inputs:
%       - image: image series to be reconstructed [nx,ny,nt,1,nsl]
%       - Data: Reconstruction variables [structure]
%           - kSpace: CAIPI modulated Cartesian k-space data [nx,ny,nt,nc,1,1,nsl]
%           - mask: mask of Cartesian k-space samples [nx,ny,nt,1,1,1,nsl]
%           - SMS: CAIPI phase modulation matrix [1,1,1,1,nsl,1,nsl]
%           - first_est: preliminary STCR [nx,ny,nt,1,nsl]
%           - sens_map: coil sensitivity maps [nx,ny,1,nc,nsl]
%           - ramp_filter: radial sampling filter [nx,ny]
%           - sms_filter: SMS acquisition filter [nx,ny,nt]
%           - basis: basis vectors for subspace constraint [nt,20]
%           - reorder: contains reordering indices for reordered temporal
%                      and spatial total variation constraints
%               - nx: number of measurements in spatial x-dimension
%               - ny: number of measurements in spatial y-dimension
%               - nt: number of time frames
%               - nc: number of PCA coils
%               - nsl: number of SMS slices (multi-band = 3)
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

if weight~=0

    real_update = compute_tTV(real(Image),Data,epsilon,'real');

    imag_update = compute_tTV(imag(Image),Data,epsilon,'imag');

    tTV_update = weight .* (real_update + 1i.*imag_update);

else
    tTV_update = 0;
end

end

function tTV_update = compute_tTV(Image,Data,epsilon,mode)

    if contains(mode,'real')
        img = Image(Data.reorder.temporal_real);
    else
        img = Image(Data.reorder.temporal_imag);
    end

    temp_a = diff(img,1,3); clear img
    temp_b = temp_a./sqrt(abs(temp_a).^2 + epsilon); clear temp_a
    temp_c = diff(temp_b,1,3);
    tTV_update = cat(3,temp_b(:,:,1,:,:,:),temp_c,-temp_b(:,:,end,:,:,:)); clear temp_b temp_c

    if contains(mode,'real')
        tTV_update(Data.reorder.temporal_real) = tTV_update;
    else
        tTV_update(Data.reorder.temporal_imag) = tTV_update;
    end

end
