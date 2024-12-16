function Image = get_MBI(Image,para)

%--------------------------------------------------------------------------
%   Image = get_MBI(Image,para)
%--------------------------------------------------------------------------
%   Function to generate model-based images based off of the input image
%   for self-gating and registraiton.
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
%--------------------------------------------------------------------------

[nx,ny,nt,nset,nsl] = size(Image);

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

parfor i=1:nsl
    for j=1:nset
        MBI(:,:,j,i) = Dict(i,:,idx(:,:,j,i));
    end
end

MBI = reshape(MBI,[nt,nx,ny,nset,nsl]);
MBI = permute(MBI,[2,3,1,4,5]);

Image(center_x,center_y,:,:,:) = MBI(center_x,center_y,:,:,:) .* norm .* phase; clear MBI norm phase

Image = reshape(Image,[nx,ny,nt,nset,nsl]);

end