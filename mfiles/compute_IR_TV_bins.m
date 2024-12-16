function update = compute_IR_TV_bins(Image,weight,epsilon,para)

%--------------------------------------------------------------------------
%   update = compute_IR_TV_bins(Image,weight,epsilon,para)
%--------------------------------------------------------------------------
%   computes the inversion total variation update term which minimizes the
%   difference between images acquired at the same inversion time
%--------------------------------------------------------------------------
%   Inputs:
%       - image: image series to be reconstructed [nx,ny,nt,1,nsl]
%       - weight: temporal total variation regularization weight [scalar]
%       - epsilon: small value to prevent singularity [scalar]
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - update: inversion total vaiaration update term [nx,ny,nt,1,nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nt: number of time frames
%           - nsl: number of SMS slices (multi-band = 3)
%--------------------------------------------------------------------------

if weight ~= 0

    [sx,sy,nt,nset,nsl] = size(Image);

    InversionTimes = para.self_gating.InversionTimes; time_mask = unique(InversionTimes);
    inversion_bins = para.self_gating.inversion_bins; inversion_mask = unique(inversion_bins);

    update  = zeros(size(Image),class(Image));
    for i=1:length(inversion_mask)
        for j=1:length(time_mask)

            if sum(inversion_bins == inversion_mask(i) & InversionTimes == time_mask(j)) ~= 0
                update(:,:,inversion_bins == inversion_mask(i) & InversionTimes == time_mask(j),:,:) = compute_tv(Image(:,:,inversion_bins == inversion_mask(i) & InversionTimes == time_mask(j),:,:),epsilon);
            end
            
        end
    end

    update = weight .* reshape(update,[sx,sy,nt,nset,nsl]);

else
    update = 0;
end

end

function update = compute_tv(Image,epsilon)

if size(Image,3) > 2
    temp_a = diff(Image,1,3);
    temp_b = temp_a./sqrt(temp_a.^2 + epsilon); clear temp_a
    temp_c = diff(temp_b,1,3);
    update = cat(3,temp_b(:,:,1,:,:),temp_c,-temp_b(:,:,end,:,:)); clear temp_b temp_c
else
    update = 0;
end

end
