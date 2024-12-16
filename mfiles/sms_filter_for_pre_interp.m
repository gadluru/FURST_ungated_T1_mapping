function filter = sms_filter_for_pre_interp(mask)

%--------------------------------------------------------------------------
%   filter = sms_filter_for_pre_interp(mask)
%--------------------------------------------------------------------------
%   Generates filter to compensate for oversampling of k-space center
%   due to SMS acquisition
%--------------------------------------------------------------------------
%   Inputs:      
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - filter: radial sampling filter [nx,ny,nt]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nt: number of time frames
%--------------------------------------------------------------------------

filter = (sum(mask,7));
filter(filter == 0) = 1;
filter = single(1./filter);

end