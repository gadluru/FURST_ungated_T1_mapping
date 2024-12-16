function [registered,para] = MBR_IR(Image,para,shifts,noi)

%--------------------------------------------------------------------------
%   [registered,para] = MBR_IR(Image,para,shifts,noi)
%--------------------------------------------------------------------------
%   Function to perform final rigid registration followed by deformable
%   registration based on model-based images
%--------------------------------------------------------------------------
%   Inputs:      
%       - Image: image series to be registered [nx,ny,nt,1,nsl]
%          - nx: number of measurements in spatial x-dimension
%          - ny: number of measurements in spatial y-dimension
%          - nt: number of time frames
%          - nsl: number of SMS slices (multi-band = 3)
%       - para: reconstruction parameters [structure]
%       - shifts: preliminary rigid shifts for initial rigid registration [2,noi,nt,1,nsl]
%          - noi: number of iterations
%          - nt: number of time frames
%          - nsl: number of SMS slices (multi-band = 3)
%       - noi: number of iterations for final registration [scalar]
%--------------------------------------------------------------------------
%   Outputs:
%       - registered: final registered images [nx,ny,nt,1,nsl]
%          - nx: spatial x-dimension
%          - ny: spatial y-dimension
%          - nt: number of time frames
%          - nsl: number of slices (SMS, multiband=3)
%       - para: reconstruction parameters [structure]
%          - Motion.rigid.shifts: rigid shifts for final registration [2,noi,nt,1,nsl]
%             - noi: number of iterations
%             - nt: number of time frames
%             - nsl: number of SMS slices (multi-band = 3)
%          - Motion.MBR.idx: deformable registration indices for final registration [nx,ny,nt,1,nsl]
%             - nx: spatial x-dimension
%             - ny: spatial y-dimension
%             - nt: number of time frames
%             - nsl: number of slices (SMS, multiband=3)
%--------------------------------------------------------------------------

temp = apply_avg_rigid_shifts(Image,shifts); MBI = get_MBI(temp,para); registered = Image;
for i=1:noi
    [~,para.Motion.rigid.shifts(:,:,:,:,:,i)] = rigid_reg_IR_GS(registered,MBI,para);
    registered = apply_rigid_shifts(registered,para,'forward');
    [registered,para.Motion.MBR.idx(:,:,:,:,:,i)] = get_motion_IR_hFOV(registered,MBI,para,1,100);
    MBI = get_MBI(registered,para);
end
