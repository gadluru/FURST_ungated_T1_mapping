function shifted = apply_avg_rigid_shifts(Image,shifts)

%--------------------------------------------------------------------------
%   shifted = apply_avg_rigid_shifts(Image,shifts)
%--------------------------------------------------------------------------
%   Deforms the input Image based on the pre-calculated rigid
%   shifts
%--------------------------------------------------------------------------
%   Inputs:      
%       - Image: non-motion compensated images [nx,ny,nt,1,nsl]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of time frames
%           - nsl: number of slices (SMS, multiband=3)
%       - shifts: calculated x- and y- rigid shifts [2,noi,nt,1,nsl]
%           - noi: number of iterations
%           - nt: number of time frames
%           - nsl: number of slices (SMS, multiband=3)
%--------------------------------------------------------------------------
%   Outputs:
%       - shifted: motion compensated images [nx,ny,nt,1,nsl]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of time frames
%           - nsl: number of slices (SMS, multiband=3)
%--------------------------------------------------------------------------

[nx,ny,nt,nset,nsl] = size(Image);

shifts = reshape(squeeze(sum(shifts,2)),[2,nt*nset*nsl]);

Image = reshape(Image,[nx,ny,nt*nset*nsl]);

PhaseX = reshape(-2i*pi*((1:nx)-floor(nx/2)-1)/nx,nx,1);
PhaseY = reshape(-2i*pi*((1:ny)-floor(ny/2)-1)/ny,1,ny);

Image_fft = fftshift2(fft2(fftshift2(Image)));

for i=1:nt*nset*nsl
    Phase = exp(PhaseX*shifts(1,i) + PhaseY*shifts(2,i));
    Image_fft(:,:,i) = Image_fft(:,:,i).*Phase;
end

shifted = fftshift2(ifft2(fftshift2(Image_fft)));

shifted = reshape(shifted,[nx,ny,nt,nset,nsl]);

end
