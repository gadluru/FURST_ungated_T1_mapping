function [Gx, Gy] = get_Gx_Gy(kSpace_in, kx, ky)

%--------------------------------------------------------------------------
%   [Gx, Gy] = get_Gx_Gy(kSpace_in, kx, ky)
%--------------------------------------------------------------------------
%   Function that estimate GROG operators, Gx and Gy, for interpolating
%   radial data to Cartesian
%--------------------------------------------------------------------------
%   Inputs:
%       - kSpace_in:radial k-space data [nx, nr, 1, nc]
%           - nx: number of measurements along a ray
%           - nr: number of rays
%           - nc: number of PCA coils
%       - kx: k-space coordinates kx, normalized to 2x matrix size [nx, nr]
%           - nx: number of measurements along a ray
%           - nr: number of rays
%       - ky: k-space coordinates ky, normalized to 2x matrix size [nx, nr]
%           - nx: number of measurements along a ray
%           - nr: number of rays
%--------------------------------------------------------------------------
%   Outputs:
%       - Gx: GROG gridding operator along x-dimension [nc, nc]
%           - nc: number of PCA coils
%       - Gy: GROG gridding operator along y-dimension [nc, nc]
%           - nc: number of PCA coils
%--------------------------------------------------------------------------
% Copyright: University of Utah Cardiovascular MRI Group
% https://medicine.utah.edu/radiology/radiology-research/research-labs/dibella2/
%--------------------------------------------------------------------------
%   Author:
%       Ye Tian
%       E-mail: phye1988@gmail.com
%--------------------------------------------------------------------------

[sx,nor,nof,nc] = size(kSpace_in);
kSpace_in = reshape(kSpace_in,[sx,nor*nof,nc]);
kx = reshape(kx,[sx,nor*nof]);
ky = reshape(ky,[sx,nor*nof]);

idx = sum(kSpace_in(:,:,1),1) == 0;
kSpace_in(:,idx,:) = [];
kx(:,idx) = [];
ky(:,idx) = [];
n_drop = sum(idx);
N = nor*nof - n_drop;

logG = zeros(nc,nc,N);

for i=1:N
    k_back = squeeze(kSpace_in(1:end-1,i,:)).';
    k_forward = squeeze(kSpace_in(2:end,i,:)).';
    G = k_forward * pinv(k_back);
    logG(:,:,i) = logm(G);
end

dkx = kx(2,:) - kx(1,:);
dky = ky(2,:) - ky(1,:);

logGx = zeros(nc);
logGy = zeros(nc);

pinv_dkx_dky = pinv([dkx.' dky.']);

for i=1:nc
    for j=1:nc
        temp = pinv_dkx_dky * squeeze(logG(i,j,:));
        logGx(i,j) = temp(1);
        logGy(i,j) = temp(2);
    end
end

Gx = expm(logGx);
Gy = expm(logGy);
