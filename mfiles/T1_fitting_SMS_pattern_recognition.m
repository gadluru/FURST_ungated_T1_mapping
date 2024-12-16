function [Image,T1,para] = T1_fitting_SMS_pattern_recognition(Image,para)

%--------------------------------------------------------------------------
%   [Image,T1,para] = T1_fitting_SMS_pattern_recognition(Image,para)
%--------------------------------------------------------------------------
%   Function to generate model-based images and their corresponding T1 maps
%   based off of the input image.
%--------------------------------------------------------------------------
%   Inputs:      
%       - Image: image series to be reconstructed [nx,ny,nt,1,nsl]
%               - nx: number of measurements in spatial x-dimension
%               - ny: number of measurements in spatial y-dimension
%               - nt: number of time frames
%               - nc: number of PCA coils
%               - nsl: number of SMS slices (multi-band = 3)
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - Image: Model-based images [nx,ny,nt,1,nsl]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of time frames
%           - nsl: number of slices (SMS, multiband=3)
%       - T1: T1 maps determined from the input Image [nx,ny,1,1,nsl]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nsl: number of slices (SMS, multiband=3)
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------

[nx,ny,nt,nset,nsl] = size(Image);

if para.setting.ifGPU
    Image = gather(Image);
end

dx = (nx - nx/2)/2;
dy = (ny - ny/2)/2;

center_x = dx + 1:dx + nx/2;
center_y = dy + 1:dy + ny/2;

norm = sqrt(sum(crop_half_FOV(abs(Image)).^2,3));
phase = exp(1i .* crop_half_FOV(angle(Image)));

im_m = crop_half_FOV(abs(Image));
im_m = im_m./sqrt(sum(im_m.^2,3));

Dic = abs(para.fitting.Dic);
Dic = Dic./sqrt(sum(Dic.^2,2)); Dict = Dic;

[im_m,Dic] = bin_averaging(im_m,Dic,para);

%% fit
idx = single(ones(nx,ny,nset,nsl)); idxh = single(zeros(nx/2,ny/2,nset,nsl));
for i=1:nsl
    for j=1:nset
        parfor k=1:ny/2
            im = im_m(:,k,:,j,i);

            im = reshape(im,[nx/2,size(im,3)]);

            d0 = abs(im - Dic(i,:,:));

            d0 = sqrt(squeeze(sum(d0.^2,2)));

            [~,idx_T1] = min(d0,[],2);

            idxh(:,k,j,i) = idx_T1;
        end
    end
end
idx(center_x,center_y,:) = idxh; clear idxh im im_m d0 idx_T1 Dic

[T1,FA,IV] = ind2sub([length(para.T1),length(para.FA),length(para.IV)],idx);

parfor i=1:nsl
    for j=1:nset
        MBI(:,:,j,i) = Dict(i,:,idx(:,:,j,i));
    end
end
clear Dict idx

T1 = para.T1(T1); T1 = reshape(T1,[nx,ny,1,nset,nsl]);
FA = para.FA(FA); FA = reshape(FA,[nx,ny,1,nset,nsl]); para.fitting.FA = gather(crop_half_FOV(FA(:,:,:,para.Recon.order)));
IV = para.IV(IV); IV = reshape(IV,[nx,ny,1,nset,nsl]); para.fitting.IV = gather(crop_half_FOV(IV(:,:,:,para.Recon.order))); clear FA IV

MBI = reshape(MBI,[nt,nx,ny,nset,nsl]);
MBI = permute(MBI,[2,3,1,4,5]);

Image(center_x,center_y,:,:,:) = MBI(center_x,center_y,:,:,:) .* norm .* phase; clear MBI norm phase

Image = reshape(Image,[nx,ny,nt,nset,nsl]);

if para.setting.ifGPU
    Image = gpuArray(Image);
end

end