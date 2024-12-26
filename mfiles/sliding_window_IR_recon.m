function [Data_all,para_all] = sliding_window_IR_recon(kSpace_all,para)

%--------------------------------------------------------------------------
%   [Data_all,para_all] = sliding_window_IR_recon(kSpace_all,para)
%--------------------------------------------------------------------------
%   Function to perform preliminary STCR reconstruction prior to the final
%   subspace-constrained model-based STCR reconstruction
%--------------------------------------------------------------------------
%   Inputs:      
%       - kSpace_all: Processed kSpace data [nx, nr, nc]
%           - nx: number of measurements along a ray
%           - nr: number of acquired rays
%           - nc: number of PCA coils
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

%% Pre-calculate GROG gridding operators for interpolating radial data to Cartesian
for i=1:para.Recon.nset

    set = para.kSpace_info.set == i-1;
    
    kSpace_set = kSpace_all(:,set,:);
    theta_set = para.kSpace_info.angle_mod(set);
    phase_set = para.kSpace_info.phase_mod(set);
    
    for j=1:para.Recon.nSMS
        phase_SMS = phase_set == j-1;
        
        kSpace_SMS = kSpace_set(:,phase_SMS,:); kSpace_SMS = permute(kSpace_SMS,[1,2,4,3]);
        theta_SMS = theta_set(phase_SMS);
        
        [kx, ky] = get_k_coor(para.Recon.sx,theta_SMS,para.Recon.sx/2+1);
        
        if isfield(para,'trajectory_correction')
            correction = para.trajectory_correction.*permute([cos(theta_SMS);sin(theta_SMS)],[4,1,2,3]);
            correction = squeeze(sum(correction,2));
            kx = kx - correction(1,:,:);
            ky = ky - correction(2,:,:);
        end

        [para.G{i,j}.Gx, para.G{i,j}.Gy] = GROG.get_Gx_Gy(kSpace_SMS, kx, ky);
    end
end

%% Data binning and radial to Cartesian interpolation

% Bins radial data based on the specified number rays per frame and the
% number of rays between frames (only bins rays from the same inversion
% group). Binned radial rays are then interpolated onto a Cartesian grid
if para.kSpace_info.TimeBetweenImages == 0
    patch_shift = para.kSpace_info.RadialViews*para.kSpace_info.ImagesInFirstGroup;
    patch_size = para.kSpace_info.RadialViews*para.kSpace_info.ImagesInFirstGroup;
else
    patch_shift = para.kSpace_info.RadialViews;
    patch_size = para.kSpace_info.RadialViews;
end

for i=1:para.Recon.nset
    count = 0;

    set = para.kSpace_info.set == i-1;

    kSpace_set = kSpace_all(:,set,:);
    theta_set = para.kSpace_info.angle_mod(set);
    phase_set = para.kSpace_info.phase_mod(set);

    t_begin_all = 1:patch_shift:size(kSpace_set,2)-patch_size+1;

    if t_begin_all(end) + patch_size-1 < size(kSpace_set,2)
        t_begin_all(end+1) = size(kSpace_set,2)-patch_size + 1;
    end

    for j=t_begin_all
        count = count + 1;

        kSpace_block = kSpace_set(:,j:j+patch_size-1,:);
        theta_block = theta_set(j:j+patch_size-1);
        phase_block = phase_set(j:j+patch_size-1);
        
        nor_total = size(kSpace_block,2);
        
        nof = floor((nor_total-(para.Recon.nor-para.Recon.nor_sl))/para.Recon.nor_sl);
        nor_total = nof*para.Recon.nor_sl+(para.Recon.nor-para.Recon.nor_sl);
        
        kSpace_block(:,nor_total+1:end,:) = [];
        theta_block(nor_total+1:end) = [];
        phase_block(nor_total+1:end) = [];
        
        kSpace = zeros(para.Recon.sx,para.Recon.nor,nof,para.Recon.no_comp);
        theta = zeros(1,para.Recon.nor,nof);
        phase = zeros(1,para.Recon.nor,nof);

        for k=1:nof
            ray_idx = (k-1)*para.Recon.nor_sl+1:(k-1)*para.Recon.nor_sl+para.Recon.nor;

            kSpace(:,:,k,:) = kSpace_block(:,ray_idx,:);
            theta(1,:,k) = theta_block(ray_idx);
            phase(1,:,k) = phase_block(ray_idx);
        end

        [kx,ky] = get_k_coor(para.Recon.sx,theta,para.Recon.sx/2+1);
        
        if isfield(para,'trajectory_correction')
            correction = para.trajectory_correction.*permute([cos(theta);sin(theta)],[4,1,2,3]);
            correction = squeeze(sum(correction,2));
            kx = kx - correction(1,:,:);
            ky = ky - correction(2,:,:);
        end
        
        Data_all{i}.kSpace{count} = GROG.SMS_GROG(kSpace,kx,ky,phase,para.G(i,:));
        
    end
    
    Data_all{i}.kSpace = cat(3,Data_all{i}.kSpace{:}); Data_all{i}.kSpace = orientate_image(Data_all{i}.kSpace,para.kSpace_info.orintation);

    if para.Recon.nSMS == 1
        [Data_all{i},para] = get_IR_Data(Data_all{i},para);
    else
        [Data_all{i},para] = get_IR_Data_SMS(Data_all{i},para);
    end
    
end
clearvars -except Data_all para

%% Preliminary STCR reconstruction

% preliminary STCR reconstruction is performed in two data blocks to reduce
% memory requirements
siz = size(Data_all{1}.first_est); siz(4) = para.Recon.nset; Image = zeros(siz,'single');

patch_shift = floor(length(para.kSpace_info.InversionTimes)/2)*para.nof_block;
patch_size = floor(length(para.kSpace_info.InversionTimes)/2)*para.nof_block;

t_begin_all = 1:patch_shift:size(Image,3)-patch_size+1;

if t_begin_all(end) + patch_size-1 < size(Image,3)
    t_begin_all(end+1) = size(Image,3)-patch_size + 1;
end

% Preliminary STCR reconstruction replaces the initial zero-filled image
% estimate prior to the final subpsace-constrained model-based STCR
% reconstruction
for i=1:para.Recon.nset
    tracker = zeros(1,1,size(Data_all{i}.first_est,3));
    
    para.Recon.scale = max(abs(Data_all{i}.first_est(:)));
    para.Recon.weight_tTV = para.Recon.scale*para.weight_tTV;
    para.Recon.weight_sTV = para.Recon.scale*para.weight_sTV;
    for n=t_begin_all
        
        tracker(n:n+patch_size-1) = tracker(n:n+patch_size-1)+1;
        
        Data.kSpace = Data_all{i}.kSpace(:,:,n:n+patch_size-1,:,:,:,:);
        Data.mask = Data_all{i}.mask(:,:,n:n+patch_size-1,:,:,:,:);
        Data.first_est = Data_all{i}.first_est(:,:,n:n+patch_size-1,:,:);
        Data.sens_map = Data_all{i}.sens_map;
        if isfield(Data_all{i},'sms_filter')
            Data.sms_filter = Data_all{i}.sms_filter(:,:,n:n+patch_size-1);
        end
        if isfield(Data_all{i},'ramp_filter')
            Data.ramp_filter = Data_all{i}.ramp_filter;
        end
        if isfield(Data_all{i},'SMS')
            Data.SMS = Data_all{i}.SMS;
        end

        Image(:,:,n:n+patch_size-1,:,:) = Image(:,:,n:n+patch_size-1,:,:) + STCR_conjugate_gradient_IR(Data,para); clear Data
    end
    Data_all{i}.first_est = Image./tracker; para_all{i} = para;
end
