function [csm] = ismrm_estimate_csm_walsh_optimized_yt(img, smoothing)

%--------------------------------------------------------------------------
%   [csm] = ismrm_estimate_csm_walsh_optimized_yt(img, smoothing)
%--------------------------------------------------------------------------
%   Estimates relative coil sensitivity maps from a set of coil images
%   using the eigenvector method. Code is based on an original implementation 
%   by Peter Kellman, NHLBI, NIH. Code made available for the ISMRM 2013 
%   Sunrise Educational Course. The Walsh method is the same with SOS
%   combination if not smoothing.
%--------------------------------------------------------------------------
%   Inputs:      
%       - img: temporal averaged zero-filled coil images [nx,ny,nc]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nc: number of PCA coils
%       - smoothing: smoothing factor [scalar]
%--------------------------------------------------------------------------
%   Outputs:
%       - sens_map: coil sensitivity maps [nx,ny,nc]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nc: number of PCA coils
%--------------------------------------------------------------------------
%   Reference:
%       - Walsh et al. (Magn Reson Med 2000;43:682-90.)
%--------------------------------------------------------------------------
%   Michael S. Hansen michael.hansen@nih.gov
%--------------------------------------------------------------------------

if nargin < 2
    smoothing = 5;
end

ncoils = size(img, 3);

% normalize by root sum of squares magnitude
mag = sqrt(sum(img .* conj(img),3));
s_raw = bsxfun(@rdivide,img,mag+eps); clear mag;

% compute sample correlation estimates at each pixel location
Rs = permute(conj(s_raw), [1,2,4,3]).*(s_raw);

% apply spatial smoothing to sample correlation estimates (NxN convolution)
if smoothing>1
	h_smooth = ones(smoothing)/(smoothing^2); % uniform smoothing kernel
    for m = 1:ncoils
        for n = 1:ncoils
            Rs(:,:,m,n) = conv2(Rs(:,:,m,n),h_smooth,'same');
        end
    end
end

% compute dominant eigenvectors of sample correlation matrices
[csm,~] = ismrm_eig_power(Rs); % using power method

return 

function [v,d] = ismrm_eig_power(R)

%--------------------------------------------------------------------------
%   [v,d] = ismrm_eig_power(R)
%--------------------------------------------------------------------------
%   Vectorized method for calculating the dominant eigenvector based on
%   power method. Input, R, is an image of sample correlation matrices
%   where: R(y,x,:,:) are sample correlation matrices (ncoil x ncoil) for each pixel
%--------------------------------------------------------------------------
%   Inputs:      
%       - R : sample correlation matrices for each pixel [nx,ny,nc, nc]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nc: number of PCA coils
%           - nc: number of PCA coils
%--------------------------------------------------------------------------
%   Outputs:
%       - v: dominant eigenvector (coil sensitivity maps) [nx,ny,nc]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nc: number of PCA coils
%       - d: dominant maximum eigenvalue
%--------------------------------------------------------------------------
%               ***************************************
%               *  Peter Kellman  (kellman@nih.gov)   *
%               *  Laboratory for Cardiac Energetics  *
%               *  NIH NHLBI                          *
%               ***************************************

[rows,cols,ncoils,~] = size(R);
N_iterations=2;
v=ones(rows,cols,ncoils); % initialize e.v.

d=zeros(rows,cols);
for i=1:N_iterations
	v = squeeze(sum(bsxfun(@times,R,v),3));
    d = sos(v);
    d( d < eps) = eps;
    v = bsxfun(@rdivide,v,d);
end

p1=angle(conj(v(:,:,1)));
v = bsxfun(@times,v,exp(1i*p1));
v = conj(v);

return
