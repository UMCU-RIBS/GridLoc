function [ coords ] = calcCoords( c, GridSteps, dims , ROI)
% Calculate the coordinates of a grid placed with the point c as
% middlepoint of the space where the grid is going to be projected.
%Also define the GridSteps and GridSize for the density and
% size, respectively, of the mesh that's being calculated. If ROI is an
% input (set to 1), calculate as ROI mesh. Else, calculate as electrode grid.

if nargin < 4 %for elec grids
    
    if length(GridSteps)>1
        %create grid in 2d space (1mm inter ROI point distance)
        gridY = 0:GridSteps(1):(dims(1)-1)*GridSteps(1);
        gridZ = 0:GridSteps(2):(dims(2)-1)*GridSteps(2);
    else
        gridY = 0:GridSteps:(dims(1)-1)*GridSteps;
        gridZ = 0:GridSteps:(dims(2)-1)*GridSteps;
    end
    
else %for establishing ROI
    gridY = 0:GridSteps:dims(1)-GridSteps;
    gridZ = 0:GridSteps:dims(2)-GridSteps;
end


% Create a meshgrid from the z and y coordinates (defined every (z,y))
[meshY,meshZ] = meshgrid(gridY,gridZ);

%Find middlepoint of the grid 
midY = mean(meshY(1,:));
midZ = mean(meshZ(:,1));
mid_grid=[midY midZ];

%establish size of ROI
dimsMesh = size(meshY);
lngth=dimsMesh(1,1)*dimsMesh(1,2);

%Set coords grid
coords(:,1) = reshape(meshY,lngth,1);
coords(:,2) = reshape(meshZ,lngth,1);

%Calculate diff matrix based upon middle of grid
diff_mat_ROI = coords - repmat(mid_grid,[size(coords,1) 1]);

diff_mat = [];
diff_mat(:,1) = zeros(lngth,1);
diff_mat(:,2) = diff_mat_ROI(:,1);
diff_mat(:,3) = diff_mat_ROI(:,2);

%calculate ROI points in fixed x-plane
coords = [];
coords = [diff_mat(:,1)+c(1,1),diff_mat(:,2)+c(1,2),diff_mat(:,3)+c(1,3)];

end

