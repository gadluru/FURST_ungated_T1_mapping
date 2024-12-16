function [Cost,totalCost] = Cost_STCR_MBR_IR_bins_reorder(Image, fNorm, mNorm, Data, para)

%--------------------------------------------------------------------------
%   [Cost,totalCost] = Cost_STCR_MBR_IR_bins_reorder(Image, fNorm, mNorm, Data, para)
%--------------------------------------------------------------------------
%   computes the cost of the subspace-constrained model-based STCR reconstruction problem 
%--------------------------------------------------------------------------
%   Inputs:
%       - image: image series to be reconstructed [nx,ny,nt,1,nsl]
%       - fNorm: fidelity norm used for conjugate gradient step sizes [nk, nc]
%               - nk: number of acquired k-space samples
%               - nc: number of PCA coils
%       - mNorm: model-based images used for model-based norm [nx,ny,nt,1,nsl]
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
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - Cost: variable containing cost of each regularization term 
%               for each iteration [structure]
%       - totalCost: total cost of the current step of the STCR
%               reconstruction
%--------------------------------------------------------------------------

Cost = para.Cost; N = numel(Image);

fNorm = sum((abs(vec(fNorm)).^2)./prod(para.Recon.kSpace_size)/64);

mNorm = para.Recon.weight_tf .* sum(abs(vec(Image - mNorm)).^2);

img_real = real(Image); img_imag = imag(Image);

if para.Recon.weight_tTV ~= 0

    tNorm1 = diff(img_real(Data.reorder.temporal_real),1,3);
    tNorm2 = diff(img_imag(Data.reorder.temporal_imag),1,3);

    tNorm = para.Recon.weight_tTV .* sum(abs(vec(tNorm1 + 1i.*tNorm2))); clear tNorm1 tNorm2
else
    tNorm = 0;
end

if para.Recon.weight_sTV ~= 0

    sx_real_norm = diff(img_real(Data.reorder.spatial_x_real),1,2); sx_real_norm(:,end+1,:,:,:) = 0;
    sy_real_norm = diff(img_real(Data.reorder.spatial_y_real),1,1); sy_real_norm(end+1,:,:,:,:) = 0;

    sx_imag_norm = diff(img_imag(Data.reorder.spatial_x_imag),1,2); sx_imag_norm(:,end+1,:,:,:) = 0;
    sy_imag_norm = diff(img_imag(Data.reorder.spatial_y_imag),1,1); sy_imag_norm(end+1,:,:,:,:) = 0;

    sx_norm = sx_real_norm + 1i.*sx_imag_norm; sy_norm = sy_real_norm + 1i.*sy_imag_norm; clear sx_real_norm sy_real_norm sx_imag_norm sy_imag_norm

    sNorm = para.Recon.weight_sTV .* sum(sqrt(vec(abs(sx_norm)).^2 + vec(abs(sy_norm)).^2));
else
    sNorm = 0;
end

if para.Recon.weight_iTV ~= 0

    InversionTimes = para.self_gating.InversionTimes; time_mask = unique(InversionTimes);
    inversion_bins = para.self_gating.inversion_bins; inversion_mask = unique(inversion_bins);

    inNorm = 0;
    for i=1:length(inversion_mask)
        for j=1:length(time_mask)
            
            if sum(inversion_bins == inversion_mask(i) & InversionTimes == time_mask(j)) ~= 0
                inNorm = inNorm + para.Recon.weight_iTV .* sum(abs(vec(diff(Image(:,:,inversion_bins == inversion_mask(i) & InversionTimes == time_mask(j),:,:),1,3))));
            end
        end
    end
else
    inNorm = 0;
end

inNorm = inNorm ./ N;
sNorm = sNorm ./ N;
tNorm = tNorm ./ N;
fNorm = fNorm ./ N;
mNorm = mNorm ./ N;

totalCost = inNorm + sNorm + tNorm + fNorm + mNorm;

if isempty(Cost.fidelityNorm)==1
    Cost.fidelityNorm = gather(fNorm);
    Cost.temporalNorm = gather(tNorm);
    Cost.inversionNorm = gather(inNorm);
    Cost.modelNorm = gather(mNorm);
    Cost.spatialNorm = gather(sNorm);
    Cost.totalCost = gather(totalCost);
else    
    Cost.fidelityNorm(end+1) = gather(fNorm);
    Cost.temporalNorm(end+1) = gather(tNorm);
    Cost.inversionNorm(end+1) = gather(inNorm);
    Cost.modelNorm(end+1) = gather(mNorm);
    Cost.spatialNorm(end+1) = gather(sNorm);
    Cost.totalCost(end+1) = gather(totalCost);
end

end
