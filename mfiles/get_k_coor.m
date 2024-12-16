function [xcoor, ycoor] = get_k_coor(nx,theta,kCenter)

%--------------------------------------------------------------------------
%   [xcoor, ycoor] = get_k_coor(sx,theta,kCenter)
%--------------------------------------------------------------------------
%   Function to return the Cartesian coordinates of radial kSpace data
%--------------------------------------------------------------------------
%   Inputs:      
%       - nx: number of measurements along a ray [scalar]
%       - theta: radial k-space projection angles [vector]
%       - kCenter: center of kSpace [scalar]
%--------------------------------------------------------------------------
%   Outputs:     
%       xcoor: k-space coordinates kx, normalized to 2x matrix size [nx,nr]
%           - nx: number of measurements along a ray
%           - nr: number of rays
%       ycoor: k-space coordinates ky, normalized to 2x matrix size [nx,nr]
%           - nx: number of measurements along a ray
%           - nr: number of rays
%--------------------------------------------------------------------------

xcoor = (1:nx) - kCenter;
ycoor = xcoor;
xcoor = bsxfun(@times,xcoor',cos(theta));
ycoor = bsxfun(@times,ycoor',sin(theta));

end
