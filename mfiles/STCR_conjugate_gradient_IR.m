function [Image,para] = STCR_conjugate_gradient_IR(Data,para)

%--------------------------------------------------------------------------
%   [Image,para] = STCR_conjugate_gradient_IR(Data,para)
%--------------------------------------------------------------------------
%   Preliminary STCR reconstruction using conjugate gradient iterations to
%   optimize the cost function from Euler-Lagrange derivative.
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
%               - nx: number of measurements in spatial x-dimension
%               - ny: number of measurements in spatial y-dimension
%               - nt: number of time frames
%               - nc: number of PCA coils
%               - nsl: number of SMS slices (multi-band = 3)
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - Image: preliminary STCR reconstruction [nx,ny,nt,1,nsl]
%          - nx: number of measurements in spatial x-dimension
%          - ny: number of measurements in spatial y-dimension
%          - nt: number of time frames
%          - nsl: number of SMS slices (multi-band = 3)
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------

disp('Performing iterative STCR reconstruction...');
disp('Showing progress...')

ifGPU          = para.setting.ifGPU;
weight_tTV     = para.Recon.weight_tTV;
weight_sTV     = para.Recon.weight_sTV;
epsilon        = para.Recon.epsilon;
para.Recon.step_size = para.Recon.step_size(1);

if isfield(Data,'first_est')
    new_img_x = Data.first_est;
end

for j=1:para.Recon.no_comp
    k = Data.kSpace(:,:,:,j,:,:,:);
    kSpace(:,j) = k(Data.mask);
end
Data.kSpace = kSpace; clear kSpace k

if ifGPU
    Data.kSpace    = gpuArray(Data.kSpace);
    new_img_x      = gpuArray(new_img_x);
    Data.first_est = gpuArray(Data.first_est);
    Data.sens_map  = gpuArray(Data.sens_map);
    epsilon        = gpuArray(epsilon);
    if isfield(Data,'ramp_filter')
        Data.ramp_filter = gpuArray(Data.ramp_filter);
    end
    if isfield(Data,'sms_filter')
        Data.sms_filter = gpuArray(Data.sms_filter);
    end    
end

para.Cost = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'totalCost',[]);

fidelity  = @(im) compute_fidelity_IR_new(im,Data,para);
temporal  = @(im) compute_tTV_IR(im,weight_tTV,epsilon,para);
spatial   = @(im) compute_sTV_IR(im,weight_sTV,epsilon);

for iter_no = 1:para.Recon.noi
   
    if mod(iter_no,10) == 1
        t1 = tic;
    end
    
    [update_term,fidelity_norm] = fidelity(new_img_x);
    
%% fidelity term/temporal/spatial TV

    update_term = update_term + temporal(new_img_x);
    update_term = update_term + spatial(new_img_x);

%% conjugate gradient
    if iter_no > 1
        beta = update_term(:)'*update_term(:)/(update_term_old(:)'*update_term_old(:)+epsilon);
        update_term = update_term + beta*update_term_old;
    end
    update_term_old = update_term; clear update_term
    
%% line search
    [para.Cost,~] = Cost_STCR_IR(new_img_x, fidelity_norm, para);
    step_size = line_search_IR(new_img_x, update_term_old, Data, para);

    para.Recon.step_size(iter_no) = step_size;

    new_img_x = new_img_x + step_size * update_term_old;

%% stop criteria
    if para.Recon.break && iter_no > 1
        if step_size<1e-4
            break
        end
    end
    
    if mod(iter_no,10) == 0
        fprintf(['Iteration = ' num2str(iter_no) '...']);
        para.Recon.time(iter_no) = toc(t1);
        toc(t1)
    end
end

Image = gather(new_img_x);
para.Recon.time_total = sum(para.Recon.time);
figure, plotCost(para.Cost); drawnow
fprintf(['Iterative reconstruction running time is ' num2str(para.Recon.time_total) 's' '\n'])

end
