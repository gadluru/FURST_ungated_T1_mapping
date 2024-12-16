
function para = self_gate_T1_IR(Image,MBI,para)

%--------------------------------------------------------------------------
%   para = self_gate_T1_IR(Image,MBI,para)
%--------------------------------------------------------------------------
%   Function to perform self-gating to get near-systole/near-diastole bins
%--------------------------------------------------------------------------
%   Inputs:      
%       - Image: full image series [nx,ny,nt,1,nsl]
%               - nx: number of measurements in spatial x-dimension
%               - ny: number of measurements in spatial y-dimension
%               - nt: number of time frames
%               - nsl: number of SMS slices (multi-band = 3)
%       - MBI: reference model-based images [nx,ny,nt,1,nsl]
%               - nx: number of measurements in spatial x-dimension
%               - ny: number of measurements in spatial y-dimension
%               - nt: number of time frames
%               - nsl: number of SMS slices (multi-band = 3)
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - para: reconstruction parameters [structure]
%           - self_gating.cardiac_bins: near-systole/near-diastole bins [2,nt]
%--------------------------------------------------------------------------

cardiac_signal = get_cardiac_signal_IR(Image,MBI,para);

cardiac_bins = zeros(2,size(Image,3)); window = 1;

dia = local_max(cardiac_signal);
sys = local_max(-cardiac_signal);

if min(sys) < min(dia)
    sys(1) = [];
elseif min(sys) > min(dia)
    dia(1) = [];
end

if max(sys) > max(dia)
    sys(end) = [];
elseif max(sys) < max(dia)
    dia(end) = [];
end

dia_idx = mod(dia,para.nof_block); dia_idx(dia_idx == 0) = para.nof_block;
for i=1:length(dia)
    if dia(i)-window > 0 && dia(i)+window < size(Image,3)
        if dia_idx(i) < 1+window
            cardiac_bins(1,dia(i):dia(i)+window) = 1;
        elseif dia_idx(i) > para.nof_block-window
            cardiac_bins(1,dia(i)-window:dia(i)) = 1;
        else
            cardiac_bins(1,dia(i)-window:dia(i)+window) = 1;
        end
    elseif dia(i) == 1
        cardiac_bins(1,dia(i):dia(i)+window) = 1;
    elseif dia(i) == size(Image,3)
        cardiac_bins(1,dia(i)-window:dia(i)) = 1;
    end
end

sys_idx = mod(sys,para.nof_block); sys_idx(sys_idx == 0) = para.nof_block;
for i=1:length(sys)
    if sys(i)-window > 0 && sys(i)+window < size(Image,3)
        if sys_idx(i) < 1+window
            cardiac_bins(2,sys(i):sys(i)+window) = 1;
        elseif sys_idx(i) > para.nof_block-window
            cardiac_bins(2,sys(i)-window:sys(i)) = 1;
        else
            cardiac_bins(2,sys(i)-window:sys(i)+window) = 1;
        end
    elseif sys(i) == 1
        cardiac_bins(2,sys(i):sys(i)+window) = 1;
    elseif sys(i) == size(Image,3)
        cardiac_bins(2,sys(i)-window:sys(i)) = 1;
    end
end

para.self_gating.cardiac_bins = logical(cardiac_bins);

end

