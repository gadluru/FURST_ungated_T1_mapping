function Dic = get_ungatedIR_Full_dic(para)

%--------------------------------------------------------------------------
%   Dic = get_ungatedIR_Full_dic(para)
%--------------------------------------------------------------------------
%   Dictionary estimation based on the full FURST sequence, modelling the 
%   radial GRE readout and SMS reconstruction. The dictionary
%   is a function of T1, FA (flip angle), IV (inversion efficiency), 
%   IRT (inversion recovery time), RTs (recovery time within a inversion
%   group), and RTL (recovery time between inversion groups).
%--------------------------------------------------------------------------
%   Inputs:      
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - Dic: Processed Dictionary [nsl, nt, nv]
%           - nsl: number of SMS slices (multiband=3)
%           - nt: number of time frames
%           - nv: number of IV/FA/T1 variations
%         (see get_ungatedIR_sliding_signal_dic.m)
%--------------------------------------------------------------------------

%% Data Initialization

fprintf(['Estimate Dictionary...' '\n'])

t1 = tic;

para.nor_total = para.kSpace_info.RadialViews*para.nof_new;
para.ImagesInGroup = (para.kSpace_info.ImagesInFirstGroup+para.kSpace_info.ImagesInSecondGroup+para.kSpace_info.ImagesInThirdGroup);
para.repetitions = para.nof_new/(para.ImagesInGroup);

flip_angle = para.kSpace_info.FlipAngle;
TR = para.kSpace_info.TimePerLine*1000;
short_recovery = para.kSpace_info.TimeBetweenImages*1000;
long_recovery = para.kSpace_info.RecoveryTimeBetweenGroups*1000;
inversion_number = para.inversion_number;

nset = para.Recon.nset;
nSMS = para.Recon.nSMS;

nfa = length(para.FA); FA = para.FA;
nt1 = length(para.T1); T1 = para.T1;
ninvt = length(para.IV); IV = para.IV;

settings.FA = FA;
settings.T1 = T1;
settings.IV = IV;

nor_total = para.nor_total;
nor_all = para.kSpace_info.RadialViews*para.repetitions*para.ImagesInGroup;

%% Dictionary Search

% code searches mfiles/LUT/T1/ for precalculated dictionary files 
% containing match dictionary parameters instead of recalculating 
% dictionary for each function call

file_name = ['mfiles/LUT/T1/Bloch_',num2str(inversion_number),'_ungatedT1_dic_with_phase_interleaving_',num2str(nset)','_SMS_',num2str(nSMS),'_TR_',num2str(TR),'_FA_',num2str(flip_angle),'_short_delay_',num2str(short_recovery),'_long_delay_',num2str(long_recovery),'_',num2str(nor_total),'.mat'];

currentSet = single(zeros(1,nor_all,length(T1),length(FA),length(IV)));

if ~isempty(dir(file_name))

    load(file_name,'SI_set_all','settings')

    T1 = setdiff(T1,settings.T1);
    FA = setdiff(FA,settings.FA);
    IV = setdiff(IV,settings.IV);

    if isempty(T1) && isempty(FA) && isempty(IV)
        if isempty(T1) T1 = para.T1; end
        if isempty(FA) FA = para.FA; end
        if isempty(IV) IV = para.IV; end

        FA_idx = find(ismember(settings.FA,FA) == 1);
        T1_idx = find(ismember(settings.T1,T1) == 1);
        inversion_idx = find(ismember(settings.IV,IV) == 1);

        Dic = get_ungatedIR_sliding_signal_dic(SI_set_all,para);

        Dic = Dic(para.Recon.order_back,:,T1_idx,FA_idx,inversion_idx);

        Dic = single(Dic(:,:,:));

        return
    else
        if isempty(T1) T1 = para.T1; end
        if isempty(FA) FA = para.FA; end
        if isempty(IV) IV = para.IV; end
    end

    T1 = sort(cat(2,settings.T1,T1)); T1 = unique(T1);
    FA = sort(cat(2,settings.FA,FA)); FA = unique(FA);
    IV = sort(cat(2,settings.IV,IV)); IV = unique(IV);

    FA_idx = find(ismember(FA,settings.FA) == 1);
    T1_idx = find(ismember(T1,settings.T1) == 1);
    inversion_idx = find(ismember(IV,settings.IV) == 1);

    nfa = length(FA); nt1 = length(T1); ninvt = length(IV);

    currentSet = single(zeros(size(SI_set_all,1),nor_all,nt1,nfa,ninvt));

    currentSet(:,:,T1_idx,FA_idx,inversion_idx) = SI_set_all;

end
mask = squeeze(sum(sum(currentSet,1),2));

%% Dictionary Estimation Code

parfor ifa = 1:nfa
    for invt = 1:ninvt
        if ~all(squeeze(mask(:,ifa,invt)))
            currentSet(:,:,:,ifa,invt) = get_Bloch_3_ungatedIR_dic_with_phase_interleaving_1_SMS_3(para,T1,FA(ifa),IV(invt));
        end
    end
end

SI_set_all = currentSet;

% Dictionary saving and binning based on sliding window specifications
if ~isempty(dir(file_name))
    settings.FA = FA;
    settings.T1 = T1;
    settings.IV = IV;

    save(file_name,'SI_set_all','settings','-v7.3')

    FA_idx = find(ismember(settings.FA,para.FA) == 1);
    T1_idx = find(ismember(settings.T1,para.T1) == 1);
    inversion_idx = find(ismember(settings.IV,para.IV) == 1);

    Dic = get_ungatedIR_sliding_signal_dic(SI_set_all,para);

    Dic = Dic(:,:,T1_idx,FA_idx,inversion_idx);
else
    mkdir('mfiles/LUT/T1/')

    save(file_name,'SI_set_all','settings','-v7.3')

    Dic = get_ungatedIR_sliding_signal_dic(SI_set_all,para);
end

Dic = Dic(para.Recon.order_back,:,:,:,:,:);

Dic = single(Dic(:,:,:));

toc(t1); fprintf('\n');

end
