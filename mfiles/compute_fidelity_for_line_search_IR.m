function fidelity_norm = compute_fidelity_for_line_search_IR(image,Data,para)
%--------------------------------------------------------------------------
%   [fidelity_norm] = compute_fidelity_yt_new(image, Data, para)
%--------------------------------------------------------------------------
%   Compute fidelity norm of a MRI reconstruction problem
%--------------------------------------------------------------------------
%   Inputs:      
%       - image             [sx, sy, nof, ...]
%       - Data              [structure]
%       - para              [structure]
%
%       - image             image
%       - Data              see 'help STCR_conjugate_gradient.m'
%       - para              see 'help STCR_conjugate_gradient.m'
%--------------------------------------------------------------------------
%   Output:
%       - fidelity_norm     [scalar]
%
%       - fidelity_norm     || Am - d ||_2^2
%--------------------------------------------------------------------------
%   A standard fidelity term it solves is:
%
%   || Am - d ||_2^2
%
%   see 'help STCR_conjugate_gradient.m' for more information.
%--------------------------------------------------------------------------
%   Reference:
%       [1]     Acquisition and reconstruction of undersampled radial data 
%               for myocardial perfusion MRI. JMRI, 2009, 29(2):466-473.
%--------------------------------------------------------------------------
%   Author:
%       Ye Tian
%       E-mail: phye1988@gmail.com
%--------------------------------------------------------------------------

switch para.Recon.type  
    case 'separate SMS'
        fidelity_norm = single(zeros(size(Data.kSpace)));
        for i=1:para.Recon.no_comp
            fidelity_update = image.*Data.sens_map(:,:,:,i,:,:);
            fidelity_update = fft2(fidelity_update);
            fidelity_update = sum(fidelity_update.*Data.SMS,5);
            fidelity_update = fidelity_update(Data.mask);
            fidelity_update = Data.kSpace(:,i) - fidelity_update;
            fidelity_norm(:,i) = gather(fidelity_update);
        end
        return
    case '2D'
        fidelity_norm = single(zeros(size(Data.kSpace)));
        for i=1:para.Recon.no_comp
            fidelity_update = image.*Data.sens_map(:,:,:,i,:,:);
            fidelity_update = fft2(fidelity_update);
            fidelity_update = fidelity_update(Data.mask);
            fidelity_update = Data.kSpace(:,i) - fidelity_update;
            fidelity_norm(:,i) = gather(fidelity_update);
        end
        return
end
        
end