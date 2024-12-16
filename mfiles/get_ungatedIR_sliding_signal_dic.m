
function sliding_signal = get_ungatedIR_sliding_signal_dic(SI_set_all,para)

%--------------------------------------------------------------------------
%   sliding_signal = get_ungatedIR_sliding_signal_dic(SI_set_all,para)
%--------------------------------------------------------------------------
%   Bins the dictionary rays modeled from 
%   get_Bloch_3_ungatedIR_dic_with_phase_interleaving_1_SMS_3.m using a
%   sliding window approach based on the chosen parameters for the number
%   of rays per frame and the number of rays between frames.
%--------------------------------------------------------------------------
%   Inputs: 
%       - SI_set_all: Raw Dictionary [1, nT1, nr]
%           - nT1: range of T1 used for dictionary estimation
%           - nr: number of acquired rays
%         (see get_Bloch_3_ungatedIR_dic_with_phase_interleaving_1_SMS_3.m)
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - sliding_signal: Processed Dictionary [nsl, nt, nv]
%           - nsl: number of SMS slices (multiband=3)
%           - nt: number of time frames
%           - nv: number of IV/FA/T1 variations
%--------------------------------------------------------------------------


%% Data Initialization

RadialViews = para.kSpace_info.RadialViews;
short_recovery = para.kSpace_info.TimeBetweenImages*1000;

nor = para.Recon.nor;
nor_sl = para.Recon.nor_sl;
nor_all = para.kSpace_info.RadialViews*para.repetitions*para.ImagesInGroup;

if short_recovery == 0
    block = RadialViews*para.kSpace_info.ImagesInFirstGroup;
else
    block = RadialViews;
end

nof_all = floor((block-(nor-nor_sl))/nor_sl);

sliding_signal = single(zeros(3,nof_all*length(1:block:nor_all),size(SI_set_all,3),size(SI_set_all,4),size(SI_set_all,5),size(SI_set_all,6)));

%% Models SMS reconstuction with sliding window approach

% only rays from the same inversion groups can be binned into a time frame
count = 0;
for i=1:block:nor_all

    Mzb = SI_set_all(:,i:i+block-1,:,:,:,:);

    nor_total = size(Mzb,2);

    nof = floor((nor_total-(nor-nor_sl))/nor_sl);
    nor_total = nof*nor_sl+(nor-nor_sl);

    Mzb(:,nor_total+1:end,:,:,:,:) = [];

    for iframe=1:nof
        count = count + 1;

        ray_idx = (iframe-1)*nor_sl+1:(iframe-1)*nor_sl+nor;

        t = mean(Mzb(:,ray_idx(1:3:end),:,:,:,:) + Mzb(:,ray_idx(2:3:end),:,:,:,:) + Mzb(:,ray_idx(3:3:end),:,:,:,:),2);
        sliding_signal(1,count,:,:,:,:) = t;

        t = mean(Mzb(:,ray_idx(1:3:end),:,:,:,:) + Mzb(:,ray_idx(2:3:end),:,:,:,:)*exp(-1i*4*pi/3) + Mzb(:,ray_idx(3:3:end),:,:,:,:)*exp(-1i*2*pi/3),2);
        sliding_signal(2,count,:,:,:,:) = t;

        t = mean(Mzb(:,ray_idx(1:3:end),:,:,:,:) + Mzb(:,ray_idx(2:3:end),:,:,:,:)*exp(-1i*2*pi/3) + Mzb(:,ray_idx(3:3:end),:,:,:,:)*exp(-1i*4*pi/3),2);
        sliding_signal(3,count,:,:,:,:) = t;
    end
end
end

