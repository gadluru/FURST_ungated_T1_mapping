function [Image,Dic] = bin_averaging(Image,Dic,para)

%--------------------------------------------------------------------------
%   [Image,Dic] = bin_averaging(Image,Dic,para)
%--------------------------------------------------------------------------
%   Function to average frames acquired at the same inversion times
%--------------------------------------------------------------------------
%   Inputs:      
%       - Image: image series to be reconstructed [nx,ny,nt,1,nsl]
%               - nx: number of measurements in spatial x-dimension
%               - ny: number of measurements in spatial y-dimension
%               - nt: number of time frames
%               - nc: number of PCA coils
%               - nsl: number of SMS slices (multi-band = 3)
%       - Dic: signal dictionary for estimating T1 [nsl,nt,nv]
%           - nsl: number of SMS slices (multiband=3)
%           - nt: number of time frames
%           - nv: number of IV/FA/T1 variations
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Output:
%       - Image: images averaged across inversion groups [nx,ny,nb,1,nsl]
%           - nx: spatial x-dimension
%           - ny: spatial y-dimension
%           - nb: number of time frames between the acquired inversions
%           - nsl: number of slices (SMS, multiband=3)
%       - Dic: signal dictionary averaged across inversion groups [nsl,nb,nv]
%           - nsl: number of SMS slices (multiband=3)
%           - nb: number of time frames between the acquired inversions
%           - nv: number of IV/FA/T1 variations
%--------------------------------------------------------------------------

InversionTimes = para.self_gating.InversionTimes; time_mask = unique(InversionTimes);
inversion_bins = para.self_gating.inversion_bins; inversion_mask = unique(inversion_bins);

Image_bin = {}; Dic_bin = {}; count = 0;
for i=1:length(inversion_mask)
    for j=1:length(time_mask)

        if sum(inversion_bins == inversion_mask(i) & InversionTimes == time_mask(j)) ~= 0
            count = count + 1;

            Image_bin{count} = mean(Image(:,:,inversion_bins == inversion_mask(i) & InversionTimes == time_mask(j),:,:),3);
            Dic_bin{count} = mean(Dic(:,inversion_bins == inversion_mask(i) & InversionTimes == time_mask(j),:),2);
        end

    end
end

Image = cat(3,Image_bin{:}); Dic = cat(2,Dic_bin{:});