function kSpace_cart = GROG_Dictionary_interp_new(kSpace_radial,Gx,Gy,kx,ky)

%--------------------------------------------------------------------------
%   kSpace_cart = GROG_Dictionary_interp_new(kSpace_radial,Gx,Gy,kx,ky)
%--------------------------------------------------------------------------
%   GROG interpolation function for converison of 2D MRI radial k-space
%   data to 2D MRI cartesian k-space data. Pre-calculate a dictionary to 
%   store all the possible shifts for GROG operator from -0.50 to 0.50 at 
%   0.01 interval.
%--------------------------------------------------------------------------
%   Inputs:      
%       - kSpace_radial: Processed kSpace data [nx,nr,nt,nc]
%           - nx: number of measurements along a ray
%           - nr: number of rays per frame in one SMS set
%           - nt: number of frames
%           - nc: number of PCA coils
%       - Gx: GROG gridding operator along x-dimension [nc, nc]
%           - nc: number of PCA coils
%       - Gy: GROG gridding operator along y-dimension [nc, nc]
%           - nc: number of PCA coils
%       - kx: k-space coordinates kx, normalized to 2x matrix size [nx,nr,nt]
%           - nx: number of measurements along a ray
%           - nr: number of rays per frame in one SMS set
%           - nt: number of frames
%       - ky: k-space coordinates ky, normalized to 2x matrix size [nx,nr,nt]
%           - nx: number of measurements along a ray
%           - nr: number of rays per frame in one SMS set
%           - nt: number of frames
%--------------------------------------------------------------------------
%   Output:
%       - kSpace_cart: CAIPI modulated Cartesian k-space data for one SMS slice [nx,ny,nt,nc]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nt: number of time frames
%           - nc: number of PCA coils
%--------------------------------------------------------------------------
%   Reference:
%       - Ye Tian, et al. (2019) PLOS ONE 14(2): e0211738.
%       - Nicole Seiberlich, et al. (2008) MRM 59:930-935.
%--------------------------------------------------------------------------
% Copyright: University of Utah Cardiovascular MRI Group
% https://medicine.utah.edu/radiology/radiology-research/research-labs/dibella2/

[sx,nor,nof,nc] = size(kSpace_radial);

kSpace_radial = reshape(kSpace_radial,[sx*nor*nof,1,nc]);

skx = size(kx,1);

% get shift distance
kx_round = round(kx);
ky_round = round(ky);
dx = round((kx_round - kx)*100)/100;
dy = round((ky_round - ky)*100)/100;
weight = 1 - sqrt(dx.^2 + dy.^2);

% calculate the dictionary for GROG operator
GxDict = single(zeros([nc nc 101]));
GyDict = single(zeros([nc nc 101]));

d = -0.50:0.01:0.50;
for di = 1:101
    GxDict(:,:,di) = Gx^d(di);
    GyDict(:,:,di) = Gy^d(di);
end

% calculate shift operator for all radial points
x_cart = kx_round + skx/2 + 2;
y_cart = ky_round + skx/2 + 2;

xDict = round((dx + 0.5)*100+1);
yDict = round((dy + 0.5)*100+1);

xDict = xDict(:);
yDict = yDict(:);

Gx_shift_all = GxDict(:,:,xDict);
Gx_shift_all = permute(Gx_shift_all,[3 1 2]);
Gy_shift_all = GyDict(:,:,yDict);
Gy_shift_all = permute(Gy_shift_all,[3 4 1 2]);
G_shift = bsxfun(@times,Gx_shift_all,Gy_shift_all);
G_shift = squeeze(sum(G_shift,3));

% calculate shifted radial points
k_target = bsxfun(@times,G_shift,kSpace_radial); clear kSpace_radial G*
k_target = squeeze(sum(k_target,3));
k_target = bsxfun(@times,weight(:),k_target);
k_target = reshape(k_target,[sx nor nof nc]);

indx = sub2ind([(skx+3),(skx+3),nof,1],x_cart,y_cart);

k_target = reshape(k_target,[sx*nor*nof,nc]);
weight = reshape(weight,skx*nor,nof);
indx = reshape(indx,[skx*nor,nof]);

indx = bsxfun(@plus,indx,0:(skx+3)^2:(skx+3)^2*(nof-1));
radial_num = bsxfun(@plus,(1:skx*nor).',0:skx*nor:skx*nor*(nof-1));
rad2cart = sparse(indx(:),radial_num(:),1,(skx+3)^2*nof,skx*nor*nof);
kSpace_cart = single(rad2cart*double(k_target));
weight_cart = single(rad2cart*double(weight(:)));

% apply the weight and cut the end of cart space
weight_cart(weight_cart==0)=1;
kSpace_cart = bsxfun(@rdivide,kSpace_cart,weight_cart);
kSpace_cart = reshape(kSpace_cart,[(skx+3) (skx+3) nof nc]);
kSpace_cart(1,:,:,:,:,:,:) = [];
kSpace_cart(:,1,:,:,:,:,:) = [];
kSpace_cart(end-1:end,:,:,:,:,:,:) = [];
kSpace_cart(:,end-1:end,:,:,:,:,:) = [];