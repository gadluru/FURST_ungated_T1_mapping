function cardiac_signal = get_cardiac_signal_IR(Image,MBI,para)

%--------------------------------------------------------------------------
%   cardiac_signal = get_cardiac_signal_IR(Image,MBI,para)
%--------------------------------------------------------------------------
%   Function to get cardiac signal with bandpass filtering the cardiac
%   motion frequency
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
%       - cardiac_signal: extracted cardiac signal from bandpass filtering [1,nt]
%--------------------------------------------------------------------------

Image = crop_half_FOV(abs(Image));
MBI = crop_half_FOV(abs(MBI));

[~,cardiac_mask] = get_cardiac_mask(Image,para);

dt = para.Recon.nor_sl * para.kSpace_info.TimePerLine*1000 * para.Recon.nset;
Fs = 1000/dt;

[Image,MBI] = polarity_correction(Image,MBI,para);

cardim = (Image - MBI) .* cardiac_mask;

[nx,ny,nt,nset,nsl] = size(cardim);

cardim = permute(cardim,[1,2,4,5,3]);
cardim = reshape(cardim,[nx*ny*nset*nsl,nt]);
cardim = cardim(logical(sum(cardim,2)),:);
cardim = sum(cardim,1);
cardiac_signal = bandpass(cardim,[0.5,2.2],Fs)';

end
