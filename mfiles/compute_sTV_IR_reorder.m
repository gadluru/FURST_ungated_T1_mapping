function [sTV_update] = compute_sTV_IR_reorder(Image, Data, weight, epsilon, para)

%--------------------------------------------------------------------------
%   [sTV_update] = compute_sTV_IR_reorder(Image, Data, weight, epsilon, para)
%--------------------------------------------------------------------------
%   computes the reordered spatial total variation update term 
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
%       - sTV_update: spatial total vaiaration update term [nx,ny,nt,1,nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nt: number of time frames
%           - nsl: number of SMS slices (multi-band = 3)
%--------------------------------------------------------------------------

if weight~=0

    real_update = compute_sTV(real(Image),Data,epsilon,'real');

    imag_update = compute_sTV(imag(Image),Data,epsilon,'imag');

    sTV_update = weight .* (real_update + 1i.*imag_update);

else
    sTV_update = 0;
end
end

function [sTV_update] = compute_sTV(img,Data,epsilon,mode)

siz = size(img);

if contains(mode,'real')
    img_x = img(Data.reorder.spatial_x_real);
    img_y = img(Data.reorder.spatial_y_real);
else
    img_x = img(Data.reorder.spatial_x_imag);
    img_y = img(Data.reorder.spatial_y_imag);
end

diff_x = diff(img_x,1,2);
diff_y = diff(img_y,1,1);

T1_num = cat(2,diff_x,zeros([siz(1),1,siz(3:end)]));
T2_num = cat(2,zeros([siz(1),1,siz(3:end)]),diff_x); clear diff_x
T3_num = cat(1,diff_y,zeros([1,siz(2:end)]));
T4_num = cat(1,zeros([1,siz(2:end)]),diff_y); clear diff_y

T1_den = sqrt(abs(T1_num).^2 + abs((T3_num+T4_num)/2).^2 + epsilon); clear T4_num
T3_den = sqrt(abs(T3_num).^2 + abs((T1_num+T2_num)/2).^2 + epsilon); clear T2_num

T1 = T1_num./T1_den; clear T1_den T1_num
T3 = T3_num./T3_den; clear T3_den T3_num

T2 = cat(2,zeros([siz(1),1,siz(3:end)]),T1(:,1:end-1,:,:,:,:));
T4 = cat(1,zeros([1,siz(2:end)]),T3(1:end-1,:,:,:,:,:));

update_x = T1 - T2; update_y = T3 - T4;

if contains(mode,'real')
    update_x(Data.reorder.spatial_x_real) = update_x;
    update_y(Data.reorder.spatial_y_real) = update_y;
else
    update_x(Data.reorder.spatial_x_imag) = update_x;
    update_y(Data.reorder.spatial_y_imag) = update_y;
end

sTV_update = update_x + update_y;

end
