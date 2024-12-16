function [Image,T1,para] = MBR_T1_FA_Bins(Data,para)

%--------------------------------------------------------------------------
%   [Image,T1,para] = MBR_T1_FA_Bins(Data,para)
%--------------------------------------------------------------------------
%   Final subspace-constrained model-based STCR reconstruction using 
%   conjugate gradient iterations to optimize the cost function 
%   from Euler-Lagrange derivative.
%--------------------------------------------------------------------------
%   Inputs:      
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
%       - Image: Reconstructed Images [matrix stored as structure]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nt: number of time frames
%           - nsl: number of slices (SMS, multiband=3)
%       - T1: Reconstructed T1 maps [matrix stored as structure]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nsl: number of slices (SMS, multiband=3)
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------

disp('Performing Subspace Constrained Model-based STCR reconstruction...');
disp('Showing progress...')

ifGPU          = para.setting.ifGPU;
weight_tTV     = para.Recon.weight_tTV;
weight_iTV     = para.Recon.weight_iTV;
weight_sTV     = para.Recon.weight_sTV;
weight_tf      = para.Recon.weight_tf;
epsilon        = para.Recon.epsilon;
para.Recon.step_size = para.Recon.step_size(1);

if isfield(Data,'first_est')
    new_img_x = Data.first_est;
end

for n=1:para.Recon.no_comp
    k = Data.kSpace(:,:,:,n,:,:,:);
    kSpace(:,n) = k(Data.mask);
end
Data.kSpace = kSpace; clear kSpace k

if ifGPU
    Data.kSpace   = gpuArray(Data.kSpace);
    new_img_x     = gpuArray(new_img_x);
    Data.sens_map = gpuArray(Data.sens_map);
    epsilon       = gpuArray(epsilon);
    if isfield(Data,'ramp_filter')
        Data.ramp_filter = gpuArray(Data.ramp_filter);
    end
    if isfield(Data,'sms_filter')
        Data.ramp_filter = gpuArray(Data.ramp_filter);
    end    
end

para.Cost = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'modelNorm',[],'inversionNorm',[],'totalCost',[]);

fidelity  = @(im) compute_fidelity_IR_new(im,Data,para);
temporal  = @(im) compute_tTV_IR_reorder(im,Data,weight_tTV,epsilon,para);
inversion = @(im) compute_IR_TV_bins(im,weight_iTV,epsilon,para);
spatial   = @(im) compute_sTV_IR_reorder(im,Data,weight_sTV,epsilon);
subspace  = @(im) compute_subspace_constraint(im,Data);
T1fit     = @(im,para) T1_fitting_SMS_pattern_recognition(im,para);

para.Recon.noi = 25;
for iter_no = 1:para.Recon.noi
    
    if mod(iter_no,5) == 1
        t1 = tic;
    end
    
    %% fidelity and regularization terms

    [update_term,fidelity_norm] = fidelity(new_img_x);

    new_img_x = apply_IR_registration(new_img_x,para,'forward'); update_term = apply_IR_registration(update_term,para,'forward');

    [MBI,~,~] = T1fit(new_img_x,para);
    update_term = update_term + weight_tf * (MBI - new_img_x);

    update_term = update_term + inversion(new_img_x);
    update_term = update_term + temporal(new_img_x);
    update_term = update_term + spatial(new_img_x);

%% conjugate gradient
    if iter_no > 1
        beta = update_term(:)'*update_term(:)/(update_term_old(:)'*update_term_old(:)+epsilon);
        update_term = update_term + beta*update_term_old;
    end
    update_term_old = update_term; clear update_term

%% line search  
    
    [para.Cost,~] = Cost_STCR_MBR_IR_bins_reorder(new_img_x,fidelity_norm,MBI,Data,para);
    step_size = line_search_MBR_IR_bins_reorder(new_img_x,update_term_old,MBI,Data,para);

    para.step_size(iter_no) = step_size;

    new_img_x = new_img_x + step_size * update_term_old;

    new_img_x = subspace(new_img_x);

    new_img_x = apply_IR_registration(new_img_x,para,'backward');

%% stop criteria
    if para.Recon.break && iter_no > 1
        if step_size<1e-4
            break
        end
    end

    if mod(iter_no,5) == 0
        fprintf(['Iteration = ' num2str(iter_no) '...']);
        para.Recon.time(iter_no) = toc(t1);
        toc(t1)
    end
end
Image = apply_IR_registration(new_img_x,para,'forward'); 

[~,T1,para] = T1fit(Image,para);

Image = gather(crop_half_FOV(abs(Image(:,:,:,para.Recon.order))));
T1 = gather(crop_half_FOV(abs(T1(:,:,:,para.Recon.order))));

para.Recon.time_total = sum(para.Recon.time);
figure, plotCost(para.Cost); drawnow
fprintf(['Iterative reconstruction running time is ' num2str(para.Recon.time_total) 's' '\n'])

end
