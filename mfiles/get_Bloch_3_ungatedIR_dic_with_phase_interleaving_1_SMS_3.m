function SI_set_all = get_Bloch_3_ungatedIR_dic_with_phase_interleaving_1_SMS_3(para,T1,flip_angle_ratio,inversion_ratio)

%--------------------------------------------------------------------------
%   SI_set_all = get_Bloch_3_ungatedIR_dic_with_phase_interleaving_1_SMS_3(para)
%--------------------------------------------------------------------------
%   Dictionary estimation based on the full FURST sequence, modelling the 
%   radial GRE readout and SMS reconstruction. The dictionary
%   is a function of T1, FA (flip angle), IV (inversion efficiency), 
%   IRT (inversion recovery time), RTs (recovery time within a inversion
%   group), and RTL (recovery time between inversion groups).
%--------------------------------------------------------------------------
%   Inputs:      
%       - para: reconstruction parameters [structure]
%       - T1: range of T1 for estimation [vector, 1:10:2000]
%       - flip_angle_ratio: flip angle B1 variation [scalar, 0.8:0.05:1.2]
%       - inversion_ratio: inversion efficiency variation [scalar, 1]
%--------------------------------------------------------------------------
%   Outputs:
%       - SI_set_all: Dictionary [1, nT1, nr]
%           - nT1: range of T1 used for dictionary estimation
%           - nr: number of acquired rays
%--------------------------------------------------------------------------

%% Data Initialization

RadialViews = para.kSpace_info.RadialViews;
repetitions = para.repetitions;
InversionTime = para.kSpace_info.InversionTimes;
flip_angle = para.kSpace_info.FlipAngle*flip_angle_ratio;

nfg = para.kSpace_info.ImagesInFirstGroup;
nsg = para.kSpace_info.ImagesInSecondGroup;
tsg = para.kSpace_info.ImagesInThirdGroup;

TR = para.kSpace_info.TimePerLine*1000;
short_recovery = para.kSpace_info.TimeBetweenImages*1000;
long_recovery = para.kSpace_info.RecoveryTimeBetweenGroups*1000;

%Slice profile used for dictionary estimation
load('mfiles/SliceProfile/SliceFast.mat','Slice')
slice_profile_excitation = Slice.Profile(1:1000);
slice_phase_excitation = exp(-1i.*Slice.Phase(1:1000));

if InversionTime(1) < 1
    InversionTime = InversionTime*1000;
end
InversionTime = reshape(InversionTime,[length(InversionTime)/repetitions,repetitions]);

flip_angle_slice = flip_angle.*slice_profile_excitation;

Slab_T1 = ones(2200,1) .* T1;

M0 = 1;
Mi = -M0*inversion_ratio;

SI_set_all = single(zeros(1,length(T1),RadialViews*repetitions*para.ImagesInGroup));
Mz = single(zeros(2200,length(T1),1000)); Mz(:,:,1) = Mi * ones(2200,length(T1));

%% Sequence Modeling and Dictionary Estimation

count = 0; % Number of repetitions of the three inversion groups
for i=1:repetitions

    inversion_group = InversionTime(:,i);

    for j=1:nfg
        
        count = count + 1;
        
        % clears all except previous magnetization history to save memory
        idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1; M = Mz(:,:,idx); Mz(:) = 0; Mz(:,:,1) = M; clear M
        
        if j == 1
            % models inversion recovery time at the start of each inversion
            % group
            t = inversion_group(1); idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1;
            Mz(:,:,idx+1:idx+length(t)) = M0 - (M0 - Mz(:,:,idx)) .* exp(-(t)./Slab_T1);
        else
            % models recovery time within each inversion group
            t = short_recovery; idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1;
            Mz(:,:,idx+1:idx+length(t)) = M0 - (M0 - Mz(:,:,idx)) .* exp(-(t)./Slab_T1);
        end
        
        % First Inversion Group Measurement
        set = ones(1,RadialViews);
        phase = vec(repmat([0,1,2],[1,RadialViews]))';
        SI_set = zeros(1,length(T1),RadialViews);
        for k=1:RadialViews
            flip_angle_all = zeros(2200,1);
            
            SMS_slices = [1:1000,401:1400,801:1800];
            
            % Models slice excitation
            flip_angle_all(SMS_slices(1:1000)) = flip_angle_all(SMS_slices(1:1000)) + flip_angle_slice;
            flip_angle_all(SMS_slices(1001:2000)) = flip_angle_all(SMS_slices(1001:2000)) + flip_angle_slice;
            flip_angle_all(SMS_slices(2001:3000)) = flip_angle_all(SMS_slices(2001:3000)) + flip_angle_slice;
            
            idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1;
            Mz(:,:,idx+1) = Mz(:,:,idx).*cos(flip_angle_all/180*pi);
            
            switch phase(k)
                case 0
                    n_phase = [0,0,0];
                case 1
                    n_phase = [0,2,1];
                case 2
                    n_phase = [0,1,2];
            end
            
            % Models SMS phase encoding
            N = ceil(k);
            for iii = [1,2,3]-1+set(k)
                n_SMS_slice = ceil(iii);
                
                slice_location = (1:1000) + 400*(iii-1);
                
                idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1;
                SI_set(set(k),:,N) = SI_set(set(k),:,N) + sum(Mz(slice_location,:,idx-1) .* sin(flip_angle_slice/180*pi) .* slice_phase_excitation) .* exp(1i*n_phase(n_SMS_slice)*2*pi/3);
            end
            
            idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1;
            Mz(:,:,idx+1) = M0 - (M0 - Mz(:,:,idx)) .* exp(-TR./Slab_T1);
            
        end
        idx = find(squeeze(SI_set_all(1,1,:)) == 0); idx = idx(1)-1;
        SI_set_all(:,:,idx+1:idx+size(SI_set,3)) = SI_set;
    end
    
    % Models the recovery time between each inversion group
    t = long_recovery; idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1;
    Mz(:,:,idx+1:idx+length(t)) = M0 - (M0 - Mz(:,:,idx)) .* exp(-(t)./Slab_T1);
    
    % Models inversion pulse with each inversion group
    idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1;
    Mz(:,:,idx+1) = - Mz(:,:,idx)*inversion_ratio;
    
    % Second Inversion Group Measurement
    for j=nfg+1:nfg+nsg
        
        count = count + 1;
        
        idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1; M = Mz(:,:,idx); Mz(:) = 0; Mz(:,:,1) = M; clear M

        if j == nfg+1
            % models inversion recovery time at the start of each inversion
            % group
            t = inversion_group(2); idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1;
            Mz(:,:,idx+1:idx+length(t)) = M0 - (M0 - Mz(:,:,idx)) .* exp(-(t)./Slab_T1);
        else
            % models recovery time within each inversion group
            t = short_recovery; idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1;
            Mz(:,:,idx+1:idx+length(t)) = M0 - (M0 - Mz(:,:,idx)) .* exp(-(t)./Slab_T1);
        end
        
        set = ones(1,RadialViews);
        phase = vec(repmat([0,1,2],[1,RadialViews]))';
        SI_set = zeros(1,length(T1),RadialViews);
        for k=1:RadialViews
            flip_angle_all = zeros(2200,1);
            
            SMS_slices = [1:1000,401:1400,801:1800];

            % Models slice excitation
            flip_angle_all(SMS_slices(1:1000)) = flip_angle_all(SMS_slices(1:1000)) + flip_angle_slice;
            flip_angle_all(SMS_slices(1001:2000)) = flip_angle_all(SMS_slices(1001:2000)) + flip_angle_slice;
            flip_angle_all(SMS_slices(2001:3000)) = flip_angle_all(SMS_slices(2001:3000)) + flip_angle_slice;
            
            idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1;
            Mz(:,:,idx+1) = Mz(:,:,idx).*cos(flip_angle_all/180*pi);
            
            switch phase(k)
                case 0
                    n_phase = [0,0,0];
                case 1
                    n_phase = [0,2,1];
                case 2
                    n_phase = [0,1,2];
            end
            
            % Models SMS phase encoding
            N = ceil(k);
            for iii = [1,2,3]-1+set(k)
                n_SMS_slice = ceil(iii);
                
                slice_location = (1:1000) + 400*(iii-1);
                
                idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1;
                SI_set(set(k),:,N) = SI_set(set(k),:,N) + sum(Mz(slice_location,:,idx-1) .* sin(flip_angle_slice/180*pi) .* slice_phase_excitation) .* exp(1i*n_phase(n_SMS_slice)*2*pi/3);
            end
            
            idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1;
            Mz(:,:,idx+1) = M0 - (M0 - Mz(:,:,idx)) .* exp(-TR./Slab_T1);
            
        end
        idx = find(squeeze(SI_set_all(1,1,:)) == 0); idx = idx(1)-1;
        SI_set_all(:,:,idx+1:idx+size(SI_set,3)) = SI_set;
    end
    
    % Models the recovery time between each inversion group
    t = long_recovery; idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1;
    Mz(:,:,idx+1:idx+length(t)) = M0 - (M0 - Mz(:,:,idx)) .* exp(-(t)./Slab_T1);
    
    % Models inversion pulse with each inversion group
    idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1;
    Mz(:,:,idx+1) = - Mz(:,:,idx)*inversion_ratio;
    
    % Third Inversion Group Measurement
    for j=nfg+nsg+1:nfg+nsg+tsg
        
        count = count + 1;
        
        idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1; M = Mz(:,:,idx); Mz(:) = 0; Mz(:,:,1) = M; clear M

        if j == nfg+nsg+1
            % models inversion recovery time at the start of each inversion
            % group
            t = inversion_group(3); idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1;
            Mz(:,:,idx+1:idx+length(t)) = M0 - (M0 - Mz(:,:,idx)) .* exp(-(t)./Slab_T1);
        else
            % models recovery time within each inversion group
            t = short_recovery; idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1;
            Mz(:,:,idx+1:idx+length(t)) = M0 - (M0 - Mz(:,:,idx)) .* exp(-(t)./Slab_T1);
        end
        
        set = ones(1,RadialViews);
        phase = vec(repmat([0,1,2],[1,RadialViews]))';
        SI_set = zeros(1,length(T1),RadialViews);
        for k=1:RadialViews
            flip_angle_all = zeros(2200,1);
            
            SMS_slices = [1:1000,401:1400,801:1800];
            
            % Models slice excitation
            flip_angle_all(SMS_slices(1:1000)) = flip_angle_all(SMS_slices(1:1000)) + flip_angle_slice;
            flip_angle_all(SMS_slices(1001:2000)) = flip_angle_all(SMS_slices(1001:2000)) + flip_angle_slice;
            flip_angle_all(SMS_slices(2001:3000)) = flip_angle_all(SMS_slices(2001:3000)) + flip_angle_slice;
            
            idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1;
            Mz(:,:,idx+1) = Mz(:,:,idx).*cos(flip_angle_all/180*pi);
            
            switch phase(k)
                case 0
                    n_phase = [0,0,0];
                case 1
                    n_phase = [0,2,1];
                case 2
                    n_phase = [0,1,2];
            end
            
            % Models SMS phase encoding
            N = ceil(k);
            for iii = [1,2,3]-1+set(k)
                n_SMS_slice = ceil(iii);
                
                slice_location = (1:1000) + 400*(iii-1);
                
                idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1;
                SI_set(set(k),:,N) = SI_set(set(k),:,N) + sum(Mz(slice_location,:,idx-1) .* sin(flip_angle_slice/180*pi) .* slice_phase_excitation) .* exp(1i*n_phase(n_SMS_slice)*2*pi/3);
            end
            
            idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1;
            Mz(:,:,idx+1) = M0 - (M0 - Mz(:,:,idx)) .* exp(-TR./Slab_T1);
            
        end
        idx = find(squeeze(SI_set_all(1,1,:)) == 0); idx = idx(1)-1;
        SI_set_all(:,:,idx+1:idx+size(SI_set,3)) = SI_set;
    end
    
    if i ~= repetitions
        % Models the recovery time between each inversion group
        t = long_recovery; idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1;
        Mz(:,:,idx+1:idx+length(t)) = M0 - (M0 - Mz(:,:,idx)) .* exp(-(t)./Slab_T1);
    
        % Models inversion pulse with each inversion group
        idx = find(squeeze(Mz(1,1,:)) == 0); idx = idx(1)-1;
        Mz(:,:,idx+1) = - Mz(:,:,idx)*inversion_ratio;
    end

end

SI_set_all = permute(SI_set_all,[1,3,2]);

