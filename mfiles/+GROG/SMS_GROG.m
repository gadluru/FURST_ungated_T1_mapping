function kSpace_cart = SMS_GROG(kSpace_all, kx, ky, phase_mod, G_pre_calculated)

%--------------------------------------------------------------------------
%   kSpace_cart = SMS_GROG(kSpace_all, kx, ky, phase_mod, G_pre_calculated)
%--------------------------------------------------------------------------
%   GROG interpolation function for converison of 2D MRI radial k-space
%   data to 2D MRI cartesian k-space data
%--------------------------------------------------------------------------
%   Inputs:      
%       - kSpace_all: Processed kSpace data [nx,nr,nt,nc]
%           - nx: number of measurements along a ray
%           - nr: number of rays per frame
%           - nt: number of frames
%           - nc: number of PCA coils
%       - kx: k-space coordinates kx, normalized to 2x matrix size [nx,nr,nt]
%           - nx: number of measurements along a ray
%           - nr: number of rays per frame
%           - nt: number of frames
%       - ky: k-space coordinates ky, normalized to 2x matrix size [nx,nr,nt]
%           - nx: number of measurements along a ray
%           - nr: number of rays per frame
%           - nt: number of frames
%       - phase_mod: k-space coordinates ky, normalized to 2x matrix size [1,nr,nt]
%           - nr: number of rays per frame
%           - nt: number of frames
%       - G_pre_calculated: GROG gridding operators [structure]
%--------------------------------------------------------------------------
%   Output:
%       - kSpace_cart: CAIPI modulated Cartesian k-space data [nx,ny,nt,nc,1,1,nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nt: number of time frames
%           - nc: number of PCA coils
%           - nsl: number of SMS slices (multi-band = 3)
%--------------------------------------------------------------------------
%   Reference:
%       - Ye Tian, et al. (2019) PLOS ONE 14(2): e0211738.
%       - Nicole Seiberlich, et al. (2008) MRM 59:930-935.
%--------------------------------------------------------------------------
% Copyright: University of Utah Cardiovascular MRI Group
% https://medicine.utah.edu/radiology/radiology-research/research-labs/dibella2/
%--------------------------------------------------------------------------
%   Author:
%       Ye Tian
%       E-mail: phye1988@gmail.com
%--------------------------------------------------------------------------

nSMS = max(phase_mod(:))+1;

[sx,nor,nof,nc,~,NSlice] = size(kSpace_all);

G = cell(1,nSMS);

for i=1:NSlice
    
    kSpace_slice = kSpace_all(:,:,:,:,i);
    
    for j=1:nSMS
        SMS = phase_mod == j-1;

        kSpace_SMS = kSpace_slice(repmat(SMS,[sx 1 1 nc]));
        kSpace_SMS = reshape(kSpace_SMS,[sx,nor/nSMS,nof,nc]);

        kx_SMS = kx(repmat(SMS,[sx 1 1]));
        ky_SMS = ky(repmat(SMS,[sx 1 1]));

        kx_SMS = reshape(kx_SMS,[sx nor/nSMS nof]);
        ky_SMS = reshape(ky_SMS,[sx nor/nSMS nof]);

        if exist('G_pre_calculated')
            G{j} = G_pre_calculated{j};
            kSpace_cart(:,:,:,:,:,:,j) = GROG.GROG_Dictionary_interp_new(kSpace_SMS,G_pre_calculated{j}.Gx,G_pre_calculated{j}.Gy,kx_SMS,ky_SMS);
        else
            [G{j}.Gx, G{j}.Gy] = GROG.get_Gx_Gy(kSpace_SMS, kx_SMS, ky_SMS);
            kSpace_cart(:,:,:,:,:,:,j) = GROG.GROG_Dictionary_interp_new(kSpace_SMS,G{j}.Gx,G{j}.Gy,kx_SMS,ky_SMS);
        end
    end
end

sx_over = size(kSpace_cart,1);
kSpace_cart = reshape(kSpace_cart,[sx_over,sx_over,nof,nc,1,NSlice,nSMS]);
if size(kSpace_cart,1)~=sx
    kSpace_cart([1,end],:,:,:,:,:,:) = [];
    kSpace_cart(:,[1,end],:,:,:,:,:) = [];
end
