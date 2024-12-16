function Recon_UCAIR_ungatedIR_self_gating(kSpace_all,para)

%--------------------------------------------------------------------------
%   Recon_UCAIR_ungated_IR_self_gating(kSpace_all,para)
%--------------------------------------------------------------------------
%   Reconstruction for Free-breathing Ungated Radial SMS Cardiac T1 Mapping
%--------------------------------------------------------------------------
%   Inputs:      
%       - kSpace_all: Processed kSpace data [nx, nr, nc]
%           - nx: number of measurements along a ray
%           - nr: number of acquired rays
%           - nc: number of PCA coils
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

%% Data Initialization and Dictionary Estimation

[para.Recon.sx,~,para.Recon.no_comp] = size(kSpace_all);

para.Recon.nset = max(para.kSpace_info.set(:))+1; % number of sets of SMS slices
para.Recon.nSMS = max(para.kSpace_info.phase_mod(:))+1; % number of SMS slices

para.Recon.order = vec((1:para.Recon.nset)'+([1,para.Recon.nSMS:-1:2]-1)*para.Recon.nset);
[~,para.Recon.order_back] = sort(para.Recon.order);

para.self_gating.InversionTimes = vec(repmat(para.kSpace_info.InversionTimes',[para.nof_block,1]))';

para.fitting.Dic = get_ungatedIR_Full_dic(para); % dictionary estimation

%% RING Trajectory Correction

set = para.kSpace_info.set == 0; set = find(set);

para.trajectory_correction = RING_IR_SMS(kSpace_all(:,set,:),para.kSpace_info.angle_mod(set),para.Recon.nSMS);

%% Preliminary Reconstruction and Data Preparation

[Data,para] = sliding_window_IR_recon(kSpace_all,para); [Data,para] = prepare_Data(Data,para); clearvars -except Data para

%% Subspace Constrained Model-based Reconstruction

for i=1:size(Data,1)
    for j=1:size(Data,2)
        para{i,j}.Recon.weight_tTV = para{i,j}.Recon.scale*para{i,j}.weight_tTV;
        para{i,j}.Recon.weight_iTV = para{i,j}.Recon.scale*para{i,j}.weight_iTV;
        para{i,j}.Recon.weight_sTV = para{i,j}.Recon.scale*para{i,j}.weight_sTV;
        [Image{i,j},T1{i,j},para{i,j}] = MBR_T1_FA_Bins(Data{i,j},para{i,j});
    end
end

for i=1:size(T1,1)
    for j=1:size(T1,2)
        figure,imagesc(T1{i,j}(:,:),[0,2500]),axis off;axis equal;colormap hot;colorbar
    end
end

save(fullfile(para{1,1}.dir.save_recon_img_dir,para{1,1}.dir.save_recon_img_name),'Image','T1','para','-v7.3')

end

