function step_size = line_search_MBR_IR_bins_reorder(old, update, update_M, Data, para)

%--------------------------------------------------------------------------
%   step_size = line_search_MBR_IR_bins_reorder(old, update, update_M, Data, para)
%--------------------------------------------------------------------------
%   This function trys to find a suitable step size to perform a CG update.
%   The function starts with a step size adopted from last iteration, and
%   multiply it by 1.3 (magic number). If the step size yeilds a cost that
%   is larger than the previous cost, it shrinks the step size by 0.8
%   (magic number again). If it yeilds a cost that is smaller than the
%   previous cost, it will increase the step size by 1.3 until it no longer
%   yeild a smaller cost. The maximum number of trys is 10.
%--------------------------------------------------------------------------
%   Inputs:
%       - old: image series from previous iteration [nx,ny,nt,1,nsl]
%       - update: conjugate gradient update [nx,ny,nt,1,nsl]
%       - update_M: model-based images used for model-based norm [nx,ny,nt,1,nsl]
%       - Data: Reconstruction variables [structure]
%           - kSpace: CAIPI modulated Cartesian k-space data [nx,ny,nt,nc,1,1,nsl]
%           - mask: mask of Cartesian k-space samples [nx,ny,nt,1,1,1,nsl]
%           - SMS: CAIPI phase modulation matrix [1,1,1,1,nsl,1,nsl]
%           - first_est: preliminary STCR [nx,ny,nt,1,nsl]
%           - sens_map: coil sensitivity maps [nx,ny,1,nc,nsl]
%           - ramp_filter: radial sampling filter [nx,ny]
%           - sms_filter: SMS acquisition filter [nx,ny,nt]
%           - basis: basis vectors for subspace constraint [nt,20]
%           - reorder: contains reordering indices for reordered temporal
%                      and spatial total variation constraints
%               - nx: number of measurements in spatial x-dimension
%               - ny: number of measurements in spatial y-dimension
%               - nt: number of time frames
%               - nc: number of PCA coils
%               - nsl: number of SMS slices (multi-band = 3)
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - step_size: conjugate gradient step size to compute updated image
%                    series [scalar]
%--------------------------------------------------------------------------

step_start = para.Recon.step_size(end)*1.3;%magic number
tau = 0.8; max_try = 10; step_size = step_start; flag = 0;

cost_old = para.Cost.totalCost(end);

for i=1:max_try
    
    new = old + step_size*update;
    
    fidelity_new = compute_fidelity_for_line_search_IR(apply_IR_registration(new,para,'backward'), Data, para);

    [~,cost_new] = Cost_STCR_MBR_IR_bins_reorder(new, fidelity_new, update_M, Data, para); clear fidelity_new
    
    if cost_new > cost_old && flag == 0
        step_size = step_size*tau;
    elseif cost_new < cost_old 
        step_size = step_size*1.3;
        cost_old = cost_new;
        flag = 1;
    elseif cost_new > cost_old && flag == 1
        step_size = step_size/1.3;
        return
    end
    
end
end