function [ ROI ] = createGrid( projectedROIpoints, rotation, turns, auxDims, subj_info, hullcortex)
%%This function creates the grids around the projected ROI points and
%rotates them according to the rotation specified. After every rotation,
%the grid is projected onto the hullcortex. It creates the grids based upon
%the gridSpecs set in the previous section. If some of the points of the
%constructed grids cannot be projected onto the brain, it will set their
%[X,Y,Z] to 0.

projectedROI = projectedROIpoints.trielectrodes;
normals = projectedROIpoints.normal; % normalvectors at projected ROI points
numElec = auxDims(1)*auxDims(2); %amount of electrodes
intElec = subj_info.intElec;

%Create ROI struct for all data points
ROI = [];
for ii = 1:length(projectedROIpoints.trielectrodes)
    for turn = 1:turns+1
        ROI(ii).coords(turn).electrodes=zeros(numElec,3);
        ROI(ii).coords(turn).trielectrodes=zeros(numElec,3);
        ROI(ii).coords(turn).normal = zeros(numElec,3);
        ROI(ii).coords(turn).pint = zeros(numElec,3);
    end
end

%set parallel pool to max 4 workers (max 4 cores available on this CURRENT PC).
%This means up to 4 grids can be projected at the same time
delete(gcp('nocreate')) %shutdown parallel pool if already running
parpool('local',4);

parfor roi_punt = 1:length(projectedROI)
    error = 0; %set for error handling
    fprintf('Projecting electrodes of ROI point %d... \n',roi_punt)
    
    %Calculate coordinates of grid for a point in the ROI
    twoDgridElectrodes = calcCoords(projectedROI(roi_punt,:),intElec,auxDims);
    
    %Calculate tangent plane to grid coordinates
    ROI(roi_punt).coords(1).electrodes = calcTangent(hullcortex,projectedROI(roi_punt,:),twoDgridElectrodes,auxDims,numElec,subj_info.hemiVect.hemi);
    
    %calc normals
    normNormals = normals(roi_punt,:)./norm(normals(roi_punt,:));
    
    %normalized matrix for rotations around a point
    for i = 1:numElec
        auxRot = AxelRot(-45,normNormals,projectedROI(roi_punt,:)) * [ROI(roi_punt).coords(1).electrodes(i,:)'; 1];
        rot(i,:) = auxRot(1:3);
    end
    
    %make rotated subj struct with locations
    ROI(roi_punt).coords(1).electrodes = rot;
    
    %Project the rotated electrodes onto the surface of the brain. For
    %explanation of the used neuralAct specifications, type 'help
    %projectElectrodes':
    try
        normdist = 25;
        intersval = 20;
        [ROI(roi_punt).coords(1)]= projectElectrodes(hullcortex, ROI(roi_punt).coords(1), normdist, 0,'fixed',intersval);
    catch
        try
            normdist = 25;
            intersval = intersval + 5;%increase by 5 mm
            [ROI(roi_punt).coords(1)]= projectElectrodes(hullcortex, ROI(roi_punt).coords(1), normdist, 0,'fixed',intersval);
        catch
            try
                intersval = intersval + 5;%increase by 5 mm
                normdist = 25;
                [ROI(roi_punt).coords(1)]= projectElectrodes(hullcortex, ROI(roi_punt).coords(1), normdist, 0,'fixed',intersval);
            catch
                ROI(roi_punt).coords(1).trielectrodes = zeros(numElec,3);
                fprintf('Could not project electrode locations of ROI point %d . \n Setting ALL values of THIS ROI POINT to 0. \n \n',roi_punt)
                error = 1;
            end
        end
    end
    if error == 0 %no need to turn if there are no projected values
        for turn=1:turns
            %preallocating for speed
            rot=zeros(numElec,3);
            
            %normalized matrix for rotations around a point
            for i = 1:numElec
                auxRot = AxelRot(rotation,normNormals,projectedROI(roi_punt,:)) * [ROI(roi_punt).coords(turn).electrodes(i,:)'; 1];
                rot(i,:) = auxRot(1:3);
            end
            %make rotated subj struct with locations
            ROI(roi_punt).coords(turn+1).electrodes = rot;
            
            %project the rotated electrodes onto the surface of the brain:
            try
                normdist = 25;
                intersval = 20;
                [ROI(roi_punt).coords(turn+1)]= projectElectrodes( hullcortex, ROI(roi_punt).coords(turn+1), normdist, 0,'fixed',intersval);
            catch
                try
                    intersval = intersval + 5;%increase by 5 mm
                    [ROI(roi_punt).coords(turn+1)]= projectElectrodes( hullcortex, ROI(roi_punt).coords(turn+1), normdist, 0,'fixed',intersval);
                catch
                    try
                        intersval = intersval + 5;%increase by 5 mm
                        [ROI(roi_punt).coords(turn+1)]= projectElectrodes( hullcortex, ROI(roi_punt).coords(turn+1), normdist, 0,'fixed',intersval);
                    catch
                        ROI(roi_punt).coords(turn).trielectrodes = zeros(numElec,3);
                        fprintf('Could not project electrode locations of ROI point %d , turn %d . \n Setting values to 0. \n \n',roi_punt, turn)
                    end
                end
            end
        end
    end
    disp('Grid positions calculated.')
end

end

