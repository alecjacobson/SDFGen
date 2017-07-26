% SDFGen Compute signed distances on a grid
%
% [phi,box] = SDFGen(V,F,dx,padding)
%
% Inputs:
%   V  #V by 3 list of vertex positions
%   F  #F by 3 list of triangle indices
%   dx  spacing between cells
%   padding number of cells beyond bounding box on either side
% Outputs:
%   phi  ny by nx by nz list of signed distance values
%   box  2 by 3 so that box(1,:) is the grid min corner and box(2,:) is the grid
%     max corner
%
% Example
%   dx = 1;
%   [phi,box] = SDFGen(V,F,dx,1);
%   [X,Y,Z] = meshgrid(box(1,1)+dx*(0:size(phi,2)-1),box(1,2)+dx*(0:size(phi,1)-1),box(1,3)+dx*(0:size(phi,3)-1));
%   surf = isosurface(X,Y,Z,phi,0);
%   clf;
%   hold on;
%   tsurf(F,V);
%   tsurf(surf.faces,surf.vertices);
%   hold off
% 
