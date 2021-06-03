function [angioMap, normAngio]=calculateAngioMap(subj_info,Tthreshold,voxelDepth,plotAngio)
%Calculation of the T-map of the angiogram. Plots a single value per
%vertice. T-threshold is equal to lower boundary when opening the Angiogram
%in the software package MRIcron. Upper boundary can be adjusted by
%changing the value in this script (see line 25), but was unadjusted during
%our tests.
%N.B.: Requires SPM-12 to run and seemingly had trouble running on OS X Sierra. 
fprintf('Calculating AngioMap...\n')
load (subj_info.neuralAct);

% vol2surf
sfile = subj_info.sfile;
tfile = subj_info.tfile;
if nargin<2
    Tthreshold = subj_info.Tthreshold;
end

% use voxplot_func_gm to calculate the angiogram t-map. 3rd and 4th
% argument set the threshold (grey/white and voxel depth, respectively)
[xyztCortex,t_surf] = voxplot_func_gm(sfile, tfile, cortex, Tthreshold, voxelDepth);
angioMap = ctmr_vox_plot(cortex, xyztCortex,t_surf,1,subj_info.hemiVect.hemi,1);

%Set upper boundary and normalize the angiomap
Z = zscore(angioMap); %zcore angio
T = Z.*(Z>0.1); %change this value to change the upper boundary of the Angio T-map

normAngio = zeros(length(angioMap),1);
for i = 1:length(angioMap)
    normAngio(i) = ((T(i)-min(T))/max(T));
    normAngio(i) = 1 - normAngio(i); %Turn around for model, so 0 is highest value and 1 lowest. 
end

%Set any angiovalues to 0 OR 1
for i = 1:length(normAngio)
    if normAngio(i) ~= 1
       normAngio(i) = 0;
    end
end

% If you want to check what the normalized angioMap looks like, set
% plotAngio to 1
if plotAngio == 1
    
    if ~exist('cortex','var')
        load (subj_info.neuralAct)
    end
    
    FV = [];
    FV.vertices = cortex.vert;
    FV.faces    = cortex.tri;
    
    figure,
    title('Angiomap')
    colormap 'hot'
    caxis([min(normAngio) max(normAngio)]); % color overlay the gaussian curvature
    mesh_h=patch(FV,'FaceVertexCdata',normAngio,'facecolor','interp','EdgeAlpha',0);
    %set some visualization properties
    if subj_info.hemiVect.hemi == 'l'
        view([270,15]);
        light('Position', [-100, 0, 0], 'Style', 'infinite');
    else
        view([90,15]);
        light('Position', [100, 0, 0], 'Style', 'infinite');
    end
    material 'dull';
    lighting gouraud
    colorbar();
    axis off
    axis normal
    daspect([1 1 1])
end
end






