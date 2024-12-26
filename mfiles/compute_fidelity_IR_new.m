function [fidelity_update, fidelity_norm] = compute_fidelity_IR_new(image, Data, para)

%--------------------------------------------------------------------------
%   [fidelity_update, fidelity_norm] = compute_fidelity_yt_new(image, Data, para)
%--------------------------------------------------------------------------
%   Compute fidelity update term of the MRI reconstruction problem and
%   ensures consistency to acquired k-space data.
%--------------------------------------------------------------------------
%   Inputs:
%       - image: image series to be reconstructed [nx,ny,nt,1,nsl]
%       - Data: Reconstruction variables [structure]
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
%   Outputs:
%       - fidelity_update: data fidelity update term [nx,ny,nt,1,nsl]
%               - nx: number of measurements in spatial x-dimension
%               - ny: number of measurements in spatial y-dimension
%               - nt: number of time frames
%               - nsl: number of SMS slices (multi-band = 3)
%       - fidelity_norm: fidelity norm used conjugate gradient step sizes [nk, nc]
%               - nk: number of acquired k-space samples
%               - nc: number of PCA coils
%--------------------------------------------------------------------------
%   Reference:
%       - Acquisition and reconstruction of undersampled radial data 
%         for myocardial perfusion MRI. JMRI, 2009, 29(2):466-473.
%--------------------------------------------------------------------------
%   Author:
%       Ye Tian
%       E-mail: phye1988@gmail.com
%--------------------------------------------------------------------------

switch para.Recon.type
    
    case 'separate SMS'
        fidelity_update_all = single(zeros(size(image)));
        fidelity_norm = single(zeros(size(Data.kSpace)));
        for i=1:para.Recon.no_comp
            fidelity_update = image.*Data.sens_map(:,:,:,i,:,:);
            fidelity_update = fft2(fidelity_update);
            fidelity_update = sum(fidelity_update.*Data.SMS,5);
            fidelity_update_temp = single(zeros(size(fidelity_update)));
            fidelity_update = fidelity_update(Data.mask);
            fidelity_update = Data.kSpace(:,i) - fidelity_update;
            fidelity_norm(:,i) = gather(fidelity_update);
            fidelity_update_temp(Data.mask) = fidelity_update;
            fidelity_update = gpuArray(fidelity_update_temp); clear fidelity_update_temp
            fidelity_update = sum(fidelity_update.*conj(Data.SMS),7);
            
            if isfield(Data,'ramp_filter')
                fidelity_update = bsxfun(@times,fidelity_update,Data.ramp_filter);
            end
            
            if isfield(Data,'sms_filter')
                fidelity_update = bsxfun(@times,fidelity_update,Data.sms_filter);
            end
            
            fidelity_update = ifft2(fidelity_update);
            fidelity_update_all = fidelity_update_all + bsxfun(@times,fidelity_update,conj(Data.sens_map(:,:,:,i,:,:)));
        end
        fidelity_update = fidelity_update_all; clear fidelity_update_all
        return
        
    case '2D'
        
        fidelity_update_all = single(zeros(size(image)));
        fidelity_norm = single(zeros(size(Data.kSpace)));
        for i=1:para.Recon.no_comp
            fidelity_update = image.*Data.sens_map(:,:,:,i,:,:);
            fidelity_update = fft2(fidelity_update);
            fidelity_update_temp = single(zeros(size(fidelity_update)));
            fidelity_update = fidelity_update(Data.mask);
            fidelity_update = Data.kSpace(:,i) - fidelity_update;
            fidelity_norm(:,i) = gather(fidelity_update);
            fidelity_update_temp(Data.mask) = fidelity_update;
            fidelity_update = gpuArray(fidelity_update_temp); clear fidelity_update_temp

            if isfield(Data,'ramp_filter')
                fidelity_update = bsxfun(@times,fidelity_update,Data.ramp_filter);
            end
            
            fidelity_update = ifft2(fidelity_update);
            fidelity_update_all = fidelity_update_all + bsxfun(@times,fidelity_update,conj(Data.sens_map(:,:,:,i,:,:)));
        end
        fidelity_update = fidelity_update_all; clear fidelity_update_all
        return
  
end
