
%% Free-breathing Ungated Radial Simultaneous Multi-Slice Cardiac T1 Mapping (FURST)
% This code performs the FURST reconstruction from the JMRI publication,
% Le JV, Mendes JK, Sideris K, et al. J Magn Reson Imaging 2024. This code
% was developed and tested on a RockLinux 8.6 operating system, with an AMD
% EPYC Milan 7543 32 2.8GHz 258 MB cache, 512 GB RAM and Nvidia A100 gpus.
% Systems with ~64 GB RAM may encounter issues with insufficient memory 
% that require code changes to reduce memory requirements during certain 
% function calls.
%--------------------------------------------------------------------------
%   Author:
%       Johnathan Le
%       E-mail: le.johnv@outlook.com
%--------------------------------------------------------------------------

addpath mfiles/

if ~exist('RawData','dir')
    mkdir('RawData')
end

% example human dataset download
if ~isfile('RawData/example_human_dataset.mat')
    url = 'https://dataverse.harvard.edu/api/access/datafile/10775588';
    websave('RawData/example_human_dataset.mat',url);
end

all_mat = dir('RawData/example_human_dataset.mat');

for i=1:length(all_mat)
    close all
    clc
    
    para.dir.load_kSpace_name = all_mat(i).name;
    para.dir.load_kSpace_dir = [all_mat(i).folder,'/'];
    
    load([para.dir.load_kSpace_dir,para.dir.load_kSpace_name])

    [nx,nr,nc] = size(kSpace);

    para.Recon.interp_method = 'GROG';
    para.Recon.type = 'separate SMS';

    para.Recon.kSpace_center = floor(nx/2)+1;
    para.over_sampling = 1;
    para.core_size = [1,1];

    para.setting.ifGPU = 1; % This code requires ~34 GB of GPU memory, set flag to 0 to run on CPU
    para.self_gating.flag = 0; % flag for all-data/binned reconstruction [0: all-data, 1: binned]

    para.Recon.noi = 50; % Number of iterations for reconstruction 

    para.Recon.nor = 30; % Number of rays per frame, temporal footprint
    para.Recon.nor_sl = 30; % Number of rays between each frame, temporal resolution [all-data: 30, binned: 15]

    para.Recon.weight_tf = 0.05; % Model-based constraint regularization weight
    para.weight_sTV = 0.00001; % spatial TV constraint regularization weight
    para.weight_tTV = 0.008; % temporal TV constraint regularization weight
    para.weight_iTV = 0.004; % inversion TV constraint regularization weight

    % Dictionary parameters
    para.IV = single(1); % inversion efficiency
    para.FA = single(0.6:0.05:1.2); % readout flip angle, B1 variations
    para.T1 = single(1:10:2000); % range of T1

    para = prepare_para(para);
    
    cut = 1:size(kSpace,2);

    % Number of groups with different inversion recovery times
    if kSpace_info.ImagesInSecondGroup == 0
        para.inversion_number = 1;
    elseif kSpace_info.ImagesInThirdGroup == 0
        para.inversion_number = 2;
    else
        para.inversion_number = 3;
    end

    if para.kSpace_info.TimeBetweenImages == 0
        para.nof_block = floor((para.kSpace_info.RadialViews*para.kSpace_info.ImagesInFirstGroup+1-(para.Recon.nor-para.Recon.nor_sl))/para.Recon.nor_sl);
        para.nof_new = size(cut,2)/para.kSpace_info.RadialViews;
    else
        para.nof_block = floor((para.kSpace_info.RadialViews-(para.Recon.nor-para.Recon.nor_sl))/para.Recon.nor_sl)*para.kSpace_info.ImagesInFirstGroup;
        para.nof_new = size(cut,2)/para.kSpace_info.RadialViews;
    end

    Recon_UCAIR_ungatedIR_self_gating(kSpace,para);

end

