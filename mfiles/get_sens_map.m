function sens_map = get_sens_map(im, options)

%--------------------------------------------------------------------------
%   sens_map = get_sens_map(im, options)
%--------------------------------------------------------------------------
%   Generates coil sensitivity maps to calculate data consistency
%   constraint during iterative reconstruction. 
%--------------------------------------------------------------------------
%   Inputs:      
%       - im: zero-filled coil images [nx,ny,nt,nc,1,1,nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nt: number of time frames
%           - nc: number of PCA coils
%           - nsl: number of SMS slices (2D: 1, SMS: 3)
%       - options: reconstruction type flag ['2D','SMS']
%--------------------------------------------------------------------------
%   Outputs:
%       - sens_map: coil sensitivity maps [nx,ny,1,nc,nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nc: number of PCA coils
%           - nsl: number of SMS slices (2D: 1, SMS: 3)
%--------------------------------------------------------------------------

smooth = 10;

if size(im,4) == 1
    sens_map = ones(size(im,1),size(im,2),'single');
    return
end
    
switch options
    case '2D'
        im_for_sens = squeeze(sum(im,3));
        sens_map(:,:,1,:) = ismrm_estimate_csm_walsh_optimized_yt(im_for_sens,smooth);
    case 'SMS'
        [sx,sy,nof,coils,nSMS,ns] = size(im);
        im = reshape(im,[sx,sy,nof,coils,nSMS*ns]);
        im_for_sens = squeeze(sum(im,3));
        sens_map = zeros(sx,sy,1,coils,nSMS*ns);
        for i=1:nSMS*ns
            sens_map(:,:,1,:,i) = ismrm_estimate_csm_walsh_optimized_yt(im_for_sens(:,:,:,i),smooth);
        end
        sens_map = reshape(sens_map,[sx,sy,1,coils,nSMS,ns]);
end

sens_map_scale = max(abs(sens_map(:)));
sens_map = sens_map/sens_map_scale;
sens_map_conj = conj(sens_map);
sens_correct_term = 1./sum(sens_map_conj.*sens_map,4);
sens_correct_term = sqrt(sens_correct_term);
sens_map = single(bsxfun(@times,sens_correct_term,sens_map));

end