function [Data_all,para_all] = prepare_Data(Data,para)

%--------------------------------------------------------------------------
%   [Data_all,para_all] = prepare_Data(Data,para)
%--------------------------------------------------------------------------
%   Data preparation determined by all-data reconstruction or binned
%   near-systole/near-diastole reconstruction. Data preparation includes self-gating
%   (if binned reconstruction), pre-calculating rigid and deformable
%   motion compensation, pre-calculating reordering indices, pre-calculating 
%   inversion time and group indices, and pre-calculating basis functions
%   for subspace constraint.
%--------------------------------------------------------------------------
%   Inputs:      
%       - Data: Reconstruction variables [structure]
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
%   Output:
%       - Data_all: Reconstruction variables [structure]
%           - kSpace: CAIPI modulated Cartesian k-space data [nx,ny,nt,nc]
%           - mask: mask of Cartesian k-space samples [nx,ny,nt]
%           - first_est: preliminary STCR [nx,ny,nt]
%           - sens_map: coil sensitivity maps [nx,ny,1,nc]
%           - ramp_filter: radial sampling filter [nx,ny]
%           - basis: basis vectors for subspace constraint [nt,20]
%           - reorder: contains reordering indices for reordered temporal
%                      and spatial total variation constraints
%               - nx: number of measurements in spatial x-dimension
%               - ny: number of measurements in spatial y-dimension
%               - nt: number of time frames
%               - nc: number of PCA coils
%       - para_all: reconstruction parameters [structure]
%--------------------------------------------------------------------------

fprintf(['\n' 'Data Preparation for Final Reconstruction...' '\n'])

% Binned near-systole/near-diastole data preparation
t1 = tic;
for i=1:size(Data,1)
    if para{i}.self_gating.flag

        % rigid registration of corresponding images across inversion groups 
        [~,shifts] = rigid_reg_IR_GS_AVG(Data{i}.first_est,para{i});
        
        % calculates inversion group indices
        para{i}.self_gating.inversion_bins = mod(1:size(Data{i}.first_est,3),para{i}.nof_block*para{i}.inversion_number);
        para{i}.self_gating.inversion_bins(para{i}.self_gating.inversion_bins == 0) = para{i}.nof_block*para{i}.inversion_number;
 
        MBI = get_MBI(apply_avg_rigid_shifts(Data{i}.first_est,shifts),para{i});

        % self gating to get near-systole/near-diastole images
        para{i} = self_gate_T1_IR(apply_avg_rigid_shifts(Data{i}.first_est,shifts),MBI,para{i});

        % data binning for near-systole/near-diastole
        Data_all = repmat(Data,[1,size(para{i}.self_gating.cardiac_bins,1)]); para_all = repmat(para,[1,size(para{i}.self_gating.cardiac_bins,1)]);
        for j=1:size(para{i}.self_gating.cardiac_bins,1)

            Data_all{i,j}.kSpace = Data{i}.kSpace(:,:,para{i}.self_gating.cardiac_bins(j,:),:,:,:,:);
            Data_all{i,j}.mask = Data{i}.mask(:,:,para{i}.self_gating.cardiac_bins(j,:),:,:,:,:);
            Data_all{i,j}.first_est = Data{i}.first_est(:,:,para{i}.self_gating.cardiac_bins(j,:),:,:);

            para_all{i,j}.fitting.Dic = para{i}.fitting.Dic(:,para{i}.self_gating.cardiac_bins(j,:),:);

            para_all{i,j}.self_gating.inversion_bins = para{i}.self_gating.inversion_bins(:,para{i}.self_gating.cardiac_bins(j,:));
            para_all{i,j}.self_gating.InversionTimes = para{i}.self_gating.InversionTimes(:,para{i}.self_gating.cardiac_bins(j,:));

            % final rigid and deformable registration of near-systole/near-diastole images
            [~,para_all{i,j}] = MBR_IR(Data_all{i,j}.first_est,para_all{i,j},shifts(:,:,para{i}.self_gating.cardiac_bins(j,:),:,:),1);

            % get near-systole/near-diastole reordering 
            Data_all{i,j}.reorder = get_reordering_matrices(Data_all{i,j}.first_est,para_all{i,j});

            if isfield(Data{i},'sms_filter')
                Data_all{i,j}.sms_filter = Data{i}.sms_filter(:,:,para{i}.self_gating.cardiac_bins(j,:));
            end

            % generates basis functions for subspace constraint
            basis = permute(para{i}.fitting.Dic(:,para{i}.self_gating.cardiac_bins(j,:),:),[2,3,1]);
            basis = squeeze(basis(:,:));
            [u,~,~] = svd(basis,'econ');
            Data_all{i,j}.basis = u(:,1:20);
        end
    else
        % All-data data preparation
        Data_all = Data; para_all = para; 
        
        % calculates inversion group indices
        para_all{i}.self_gating.inversion_bins = mod(1:size(Data_all{i}.first_est,3),para_all{i}.nof_block*para{i}.inversion_number);
        para_all{i}.self_gating.inversion_bins(para_all{i}.self_gating.inversion_bins == 0) = para_all{i}.nof_block*para{i}.inversion_number;

        % rigid and deformable all-data registration
        [~,shifts] = rigid_reg_IR_GS_AVG(Data_all{i}.first_est,para_all{i}); [~,para_all{i}] = MBR_IR(Data_all{i}.first_est,para_all{i},shifts,1);

        % get all-data reordering  indices
        Data_all{i}.reorder = get_reordering_matrices(Data_all{i}.first_est,para_all{i});

        % generates basis functions for subspace constraint
        basis = permute(para_all{i}.fitting.Dic,[2,3,1]);
        basis = squeeze(basis(:,:));
        [u,~,~] = svd(basis,'econ');
        Data_all{i}.basis = u(:,1:20);
    end
end

toc(t1); fprintf('\n');
