function reorder = get_reordering_matrices(Image,para)

%--------------------------------------------------------------------------
%   reorder = get_reordering_matrices(Image,para)
%--------------------------------------------------------------------------
%   Function to get the reordering indices to perform reordered temporal
%   total variation and reordered spatial total variation during the final
%   subspace-constrained model-based STCR reconstruction
%--------------------------------------------------------------------------
%   Inputs:      
%       - Image: image series [nx,ny,nt,1,nsl]
%          - nx: number of measurements in spatial x-dimension
%          - ny: number of measurements in spatial y-dimension
%          - nt: number of time frames
%          - nsl: number of SMS slices (multi-band = 3)
%       - para: reconstruction parameters [structure]
%--------------------------------------------------------------------------
%   Outputs:
%       - reorder: reordering matrices [structure]
%         - temporal_real: real temporal reordering [nx,ny,nt,1,nsl]
%         - temporal_imag: imag temporal reordering [nx,ny,nt,1,nsl]
%         - spatial_x_real: real spatial x reordering [nx,ny,nt,1,nsl]
%         - spatial_y_real: real spatial y reordering [nx,ny,nt,1,nsl]
%         - spatial_y_imag: imag spatial y reordering [nx,ny,nt,1,nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nt: number of time frames
%           - nsl: number of SMS slices (multi-band = 3)
%--------------------------------------------------------------------------

Image = gather(apply_IR_registration(Image,para,'forward'));

[sx,sy,nt,nset,nsl] = size(Image);

[Y,X,T,S,L] = ndgrid(1:sx,1:sy,1:nt,1:nset,1:nsl);

%% temporal reordering

[~,temporal_reorder_real] = sort(real(Image),3);
[~,temporal_reorder_imag] = sort(imag(Image),3);

temporal_reorder_real = sub2ind(size(Image),Y,X,temporal_reorder_real,S,L);
temporal_reorder_imag = sub2ind(size(Image),Y,X,temporal_reorder_imag,S,L);

%% spatial x reordering

[~,spatial_reorder_x_real] = sort(real(Image),2);
[~,spatial_reorder_x_imag] = sort(imag(Image),2);

spatial_reorder_x_real = sub2ind(size(Image),Y,spatial_reorder_x_real,T,S,L);
spatial_reorder_x_imag = sub2ind(size(Image),Y,spatial_reorder_x_imag,T,S,L);

%% spatial y reordering

[~,spatial_reorder_y_real] = sort(real(Image),1);
[~,spatial_reorder_y_imag] = sort(imag(Image),1);

spatial_reorder_y_real = sub2ind(size(Image),spatial_reorder_y_real,X,T,S,L);
spatial_reorder_y_imag = sub2ind(size(Image),spatial_reorder_y_imag,X,T,S,L);

%%

reorder.temporal_real = int32(temporal_reorder_real);
reorder.temporal_imag = int32(temporal_reorder_imag);

reorder.spatial_x_real = int32(spatial_reorder_x_real);
reorder.spatial_x_imag = int32(spatial_reorder_x_imag);

reorder.spatial_y_real = int32(spatial_reorder_y_real);
reorder.spatial_y_imag = int32(spatial_reorder_y_imag);

end

