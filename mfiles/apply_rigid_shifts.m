function shifted = apply_rigid_shifts(Image,para,mode)

%--------------------------------------------------------------------------
%   shifted = apply_rigid_shifts(Image,para,mode)
%--------------------------------------------------------------------------
%   Deforms and reforms the input Image based on the pre-calculated rigid
%   shifts
%--------------------------------------------------------------------------
%   Inputs:      
%       - Image: motion compensated/non-motion compensated images [nx,ny,nt,1,nsl]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of time frames
%           - nsl: number of slices (SMS, multiband=3)
%       - para: reconstruction parameters [structure]
%       - mode: deforms/reforms images ['forward','backward']
%--------------------------------------------------------------------------
%   Output:
%       - shifted: motion compensated/non-motion compensated images [nx,ny,nt,1,nsl]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of time frames
%           - nsl: number of slices (SMS, multiband=3)
%--------------------------------------------------------------------------

[nx,ny,nt,nset,nsl] = size(Image);

shifts = reshape(squeeze(sum(para.Motion.rigid.shifts,2)),[2,nt*nset*nsl,size(para.Motion.rigid.shifts,6)]);

Image = reshape(Image,[nx,ny,nt*nset*nsl]);

PhaseX = reshape(-2i*pi*((1:nx)-floor(nx/2)-1)/nx,nx,1);
PhaseY = reshape(-2i*pi*((1:ny)-floor(ny/2)-1)/ny,1,ny);

Image_fft = fftshift2(fft2(fftshift2(Image)));

if contains(mode,'backward')
    shifts = -shifts;
end

for i=1:size(shifts,3)
    for j=1:nt*nset*nsl
        Phase = exp(PhaseX*shifts(1,j,i) + PhaseY*shifts(2,j,i));
        Image_fft(:,:,j) = Image_fft(:,:,j).*Phase;
    end
end

shifted = fftshift2(ifft2(fftshift2(Image_fft)));

shifted = reshape(shifted,[nx,ny,nt,nset,nsl]);

end