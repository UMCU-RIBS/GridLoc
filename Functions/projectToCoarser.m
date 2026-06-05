function [ ROI ] = projectToCoarser( ROI, cortexcoarser, turns)
%This block projects the grids that were projected onto the hullcortex,
%onto the coarser model, whilst using the normals that were used to project
%them onto the hullcortex. This means that they will follow their original
%trajectory, and will allow us to inspect weither or not they fall inside
%sulci or on top of gyrii.

parfor i  = 1:length(ROI)
    fprintf('Calculating projection of electrodes onto coarser model for ROI point %d... \n', i)
    intersval = 25;
    for ij = 1:(turns+1)
        try
            ss=[];
            ss.electrodes = ROI(i).coords(ij).trielectrodes;
            ss.normal     = ROI(i).coords(ij).normal;
            [ss] = projectElectrodes(cortexcoarser,ss,25,1,'fixed',intersval); %if it doesn't find an intersection, increase the number after 'fixed'
            pint = ss.trielectrodes;
            ROI(i).coords(ij).pint = pint;
        catch
            try
                intersval = intersval + 5;%increase by 5 mm
                ss=[];
                ss.electrodes = ROI(i).coords(ij).trielectrodes;
                ss.normal     = ROI(i).coords(ij).normal;
                [ss] = projectElectrodes(cortexcoarser,ss,25,1,'fixed',intersval); %if it doesn't find an intersection, increase the number after 'fixed'
                pint = ss.trielectrodes;
                ROI(i).coords(ij).pint = pint;
            catch
                try
                    intersval = intersval + 5;%increase by 5 mm
                    ss=[];
                    ss.electrodes = ROI(i).coords(ij).trielectrodes;
                    ss.normal     = ROI(i).coords(ij).normal;
                    [ss] = projectElectrodes(cortexcoarser,ss,25,1,'fixed',intersval); %if it doesn't find an intersection, increase the number after 'fixed'
                    pint = ss.trielectrodes;
                    ROI(i).coords(ij).pint = pint;
                catch
                    ROI(i).coords(ij).pint = ROI(i).coords(ij).pint;
                    fprintf('Could not project electrode locations of ROI point %d rotation %d . \n Setting values to 0. \n \n',i,ij)
                end
            end
        end
    end
end
end

