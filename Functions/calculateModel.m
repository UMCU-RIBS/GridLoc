function [ ROI ] = calculateModel( subj_info, ROI, cortex, normAngio)
% The next block is needed for modelling the angiogram to predict activity
% values.
useAngio = exist('normAngio','var');

if useAngio == 1
    %Find the indices of the points within (<=) 7 mm of a projected
    %electrode
    for i = 1:length(ROI)
        for ij = 1:length(ROI(1).coords)
            for ee = 1:length(ROI(1).coords(1).electrodes)
                dist = pdist2(ROI(i).coords(ij).pint(ee,:),cortex.vert,'euclidean');
                idx = find(dist<=2); %2*#mm diameter
                ROI(i).coords(ij).model(ee).idx = idx';
            end
        end
    end
    %finding the values of the angiomaps of the indices (=<2 mm)
    for i = 1:length(ROI)
        for ij = 1:length(ROI(1).coords)
            for ee = 1:length(ROI(1).coords(1).electrodes)
                idx = ROI(i).coords(ij).model(ee).idx;
                Mv = [];
                for iid=1:length(idx)
                    Mv(iid)=normAngio(idx(iid));
                end
                ROI(i).coords(ij).model(ee).Mv = Mv';
            end
        end
    end
    
    %Calculating the mean value of the angio and depth per electrode
    for i = 1:length(ROI)
        for ij = 1:length(ROI(1).coords)
            for ee = 1:length(ROI(1).coords(1).electrodes)
                MvM(ee) = mean(ROI(i).coords(ij).model(ee).Mv);
                McM(ee) = pdist2(ROI(i).coords(ij).pint(ee),ROI(i).coords(ij).trielectrodes(ee),'euclidean');
            end
            maxM=max(McM);
            minM=min(McM);
            McM =(McM-minM)/maxM;
            McM = 1-McM;
            ROI(i).coords(ij).MvM = MvM;
            ROI(i).coords(ij).McM = McM;
        end
    end
    %weighting of angio vs depth
    for i = 1:length(ROI)
        for ij = 1:length(ROI(1).coords)
            ROI(i).coords(ij).weights = (0.5*ROI(i).coords(ij).MvM + 0.5*ROI(i).coords(ij).McM)'; %weighting of the Angio and Curvature means
        end
    end
else
    %Calculating the mean depth per electrode
    for i = 1:length(ROI)
        for ij = 1:length(ROI(1).coords)
            for ee = 1:length(ROI(1).coords(1).electrodes)
                McM(ee) = pdist2(ROI(i).coords(ij).pint(ee),ROI(i).coords(ij).trielectrodes(ee),'euclidean');
            end
            maxM=max(McM);
            minM=min(McM);
            McM =(McM-minM)/maxM;
            McM = 1-McM;
            ROI(i).coords(ij).McM = McM;
        end
    end
        %weighting
    for i = 1:length(ROI)
        for ij = 1:length(ROI(1).coords)
            ROI(i).coords(ij).weights = ROI(i).coords(ij).McM'; %weighting of the Angio and Curvature means
        end
    end
end
end

