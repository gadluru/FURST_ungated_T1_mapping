function [Cost,totalCost] = Cost_STCR_IR(Image, fNorm, para)

%--------------------------------------------------------------------------
%   [Cost,totalCost] = Cost_STCR_IR(Image, fNorm, para)
%--------------------------------------------------------------------------
%   computes the cost of the STCR reconstruction problem 
%--------------------------------------------------------------------------
%   Inputs:
%       - image: image series to be reconstructed [nx,ny,nt,1,nsl]
%       - fNorm: fidelity norm used conjugate gradient step sizes [nk, nc]
%               - nk: number of acquired k-space samples
%               - nc: number of PCA coils
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - Cost: variable containing cost of each regularization term 
%               for each iteration [structure]
%       - totalCost: total cost of the current step of the STCR
%               reconstruction
%--------------------------------------------------------------------------

[nx,ny,nt,nset,nsl] = size(Image);

Cost = para.Cost; N = numel(Image);

fNorm = sum((abs(vec(fNorm)).^2)./prod(para.Recon.kSpace_size)/64);

if para.Recon.weight_tTV ~= 0
    tNorm = para.Recon.weight_tTV .* sum(abs(vec(diff(reshape(Image,[nx,ny,para.nof_block,nt/para.nof_block,nset,nsl]),1,3))));
else
    tNorm = 0;
end

if para.Recon.weight_sTV ~= 0
    sx_norm = diff(Image,1,2); sx_norm(:,end+1,:,:,:) = 0;
    sy_norm = diff(Image,1,1); sy_norm(end+1,:,:,:,:) = 0;
    
    sNorm = para.Recon.weight_sTV .* sum(sqrt(vec(abs(sx_norm)).^2 + vec(abs(sy_norm)).^2)); clear sx_norm sy_norm
else
    sNorm = 0;
end

sNorm = sNorm ./ N;
tNorm = tNorm ./ N;
fNorm = fNorm ./ N;

totalCost = sNorm + tNorm + fNorm;

if isempty(Cost.fidelityNorm)
    Cost.fidelityNorm = gather(fNorm);
    Cost.temporalNorm = gather(tNorm);
    Cost.spatialNorm = gather(sNorm);
    Cost.totalCost = gather(totalCost);
else    
    Cost.fidelityNorm(end+1) = gather(fNorm);
    Cost.temporalNorm(end+1) = gather(tNorm);
    Cost.spatialNorm(end+1) = gather(sNorm);
    Cost.totalCost(end+1) = gather(totalCost);
end
end
