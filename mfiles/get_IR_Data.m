function [Data,para] = get_IR_Data(Data,para)

%--------------------------------------------------------------------------
%   [Data,para] = get_IR_Data(Data,para)
%--------------------------------------------------------------------------
%   Data preparation prior to iterative reconstruction. Generates the
%   initial estimate, coil sensitivity maps, and filters from the
%   interpolated Cartesian k-space data.
%--------------------------------------------------------------------------
%   Inputs:      
%       - Data_all: Reconstruction variables [structure]
%           - kSpace: CAIPI modulated Cartesian k-space data [nx,ny,nt,nc]
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - Data_all: Reconstruction variables [structure]
%           - kSpace: CAIPI modulated Cartesian k-space data [nx,ny,nt,nc]
%           - mask: mask of Cartesian k-space samples [nx,ny,nt]
%           - first_est: preliminary STCR [nx,ny,nt]
%           - sens_map: coil sensitivity maps [nx,ny,1,nc]
%           - ramp_filter: radial sampling filter [nx,ny]
%               - nx: number of measurements in spatial x-dimension
%               - ny: number of measurements in spatial y-dimension
%               - nt: number of time frames
%               - nc: number of PCA coils
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------

Data.kSpace = fftshift2(Data.kSpace);
Data.mask = logical(abs(Data.kSpace(:,:,:,1,:,:)));

Data.kSpace = ifft2(Data.kSpace);
Data.kSpace = fftshift2(Data.kSpace);
Data.kSpace = fft2(Data.kSpace);

Data.kSpace = Data.kSpace.*Data.mask; kSpace = Data.kSpace;
Data.first_est = ifft2(Data.kSpace);

if ~isfield(Data,'sens_map')
    Data.sens_map = get_sens_map(Data.first_est,'2D');
end

para.Recon.kSpace_size = [size(Data.kSpace,1),size(Data.kSpace,2)];
para.Recon.image_size = [size(Data.kSpace,1),size(Data.kSpace,2)];
para.Recon.sx = size(Data.kSpace,1);

Data.ramp_filter = ramp_IR_filter_for_pre_interp(para);
kSpace = kSpace .* Data.ramp_filter;

Data.first_est = ifft2(kSpace) .* conj(Data.sens_map);
Data.first_est = sum(Data.first_est,4);

end

