function [Data,para] = get_IR_Data_SMS(Data,para)

%--------------------------------------------------------------------------
%   [Data,para] = get_IR_Data_SMS(Data,para)
%--------------------------------------------------------------------------
%   Data preparation prior to iterative reconstruction. Generates the
%   initial estimate, coil sensitivity maps, and filters from the
%   interpolated Cartesian k-space data.
%--------------------------------------------------------------------------
%   Inputs:      
%       - Data_all: Reconstruction variables [structure]
%           - kSpace: CAIPI modulated Cartesian k-space data [nx,ny,nt,nc,1,1,nsl]
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - Data_all: Reconstruction variables [structure]
%           - kSpace: CAIPI modulated Cartesian k-space data [nx,ny,nt,nc,1,1,nsl]
%           - mask: mask of Cartesian k-space samples [nx,ny,nt,1,1,1,nsl]
%           - SMS: CAIPI phase modulation matrix [1,1,1,1,nsl,1,nsl]
%           - first_est: preliminary STCR [nx,ny,nt,1,nsl]
%           - sens_map: coil sensitivity maps [nx,ny,1,nc,nsl]
%           - ramp_filter: radial sampling filter [nx,ny]
%           - sms_filter: SMS acquisition filter [nx,ny,nt]
%               - nx: number of measurements in spatial x-dimension
%               - ny: number of measurements in spatial y-dimension
%               - nt: number of time frames
%               - nc: number of PCA coils
%               - nsl: number of SMS slices (multi-band = 3)
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------

Data.kSpace = fftshift2(Data.kSpace);
Data.mask = logical(abs(Data.kSpace(:,:,:,1,:,:,:)));

Data.kSpace = ifft2(Data.kSpace);
Data.kSpace = fftshift2(Data.kSpace);
Data.kSpace = fft2(Data.kSpace);
Data.kSpace = Data.kSpace.*Data.mask;

Data.SMS = exp(1i*(0:para.Recon.nSMS-1)*2*pi/para.Recon.nSMS.*(0:para.Recon.nSMS-1).');
Data.SMS = permute(Data.SMS,[3,4,5,6,1,7,2]);

kSpace_sms = sum(Data.kSpace.*conj(Data.SMS),7);
Data.first_est = ifft2(kSpace_sms);

if ~isfield(Data,'sens_map')
    Data.sens_map = get_sens_map(Data.first_est,'SMS');
end

para.Recon.kSpace_size = [size(Data.kSpace,1),size(Data.kSpace,2)];
para.Recon.image_size = [size(Data.kSpace,1),size(Data.kSpace,2)];
para.Recon.sx = size(Data.kSpace,1);

Data.ramp_filter = ramp_IR_filter_for_pre_interp(para);
Data.sms_filter = sms_filter_for_pre_interp(Data.mask);

kSpace_sms = kSpace_sms .* Data.ramp_filter;
kSpace_sms = kSpace_sms .* Data.sms_filter;

Data.first_est = ifft2(kSpace_sms).*conj(Data.sens_map);
Data.first_est = sum(Data.first_est,4);

end
