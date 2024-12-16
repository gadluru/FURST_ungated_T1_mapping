function update = compute_subspace_constraint(Image,Data)

%--------------------------------------------------------------------------
%   update = compute_subspace_constraint(Image,Data)
%--------------------------------------------------------------------------
%   computes the subspace update which restricts the reconstructed image series
%   to the subspace defined by the basis functions, calculated by taking the
%   SVD of the signal dictionary
%--------------------------------------------------------------------------
%   Inputs:
%       - image: image series to be reconstructed [nx,ny,nt,1,nsl]
%       - Data: Reconstruction variables [structure]
%           - kSpace: CAIPI modulated Cartesian k-space data [nx,ny,nt,nc,1,1,nsl]
%           - mask: mask of Cartesian k-space samples [nx,ny,nt,1,1,1,nsl]
%           - SMS: CAIPI phase modulation matrix [1,1,1,1,nsl,1,nsl]
%           - first_est: preliminary STCR [nx,ny,nt,1,nsl]
%           - sens_map: coil sensitivity maps [nx,ny,1,nc,nsl]
%           - ramp_filter: radial sampling filter [nx,ny]
%           - sms_filter: SMS acquisition filter [nx,ny,nt]
%           - basis: basis vectors for subspace constraint [nt,20]
%           - reorder: contains reordering indices for reordered temporal
%                      and spatial total variation constraints
%               - nx: number of measurements in spatial x-dimension
%               - ny: number of measurements in spatial y-dimension
%               - nt: number of time frames
%               - nc: number of PCA coils
%               - nsl: number of SMS slices (multi-band = 3)
%--------------------------------------------------------------------------
%   Outputs:
%       - update: subspace update [nx,ny,nt,1,nsl]
%           - nx: number of measurements in spatial x-dimension
%           - ny: number of measurements in spatial y-dimension
%           - nt: number of time frames
%           - nsl: number of SMS slices (multi-band = 3)
%--------------------------------------------------------------------------

[nx,ny,nt,nset,nsl] = size(Image);

update = permute(Image,[1,2,4,5,3]);
update = reshape(update,[nx*ny*nset*nsl,nt]);
update = Data.basis'*update';
update = Data.basis*update;
update = reshape(update',[nx,ny,nset,nsl,nt]);
update = permute(update,[1,2,5,3,4]);
