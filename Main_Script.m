
% Gridloc Toolbox V1.3
%  26/04/2021
%
% UMC Utrecht
%  (C) Mariana P Branco  - m.pedrosobranco@umcutrecht.nl
%  (C) Michael Leibbrand - michael.leib7@gmail.com

warning off; %turn off to reduce clutter in command window.
clear all; close all; clc;

%Paths of (support) files/scripts
addpath(genpath('.../GridLoc/'));


%% 1- Create the subject info structure and set the location of the input files.

%set method parameters
subj_info.HFB_mean    = '.../HFB_pattern.mat'; %Set location of file with the caclulated values of mean HFB (i.e., high-frequency band) activity per electrode
subj_info.chansel       = [1:128]; %channel selection HD grid
subj_info.bad_channels  = []; %set bad channels
subj_info.dims          = [8,16]; %amount of electrodes [X,Y]. X-direction is side where the leads are located.
subj_info.intElec       = [3,3]; %inter-electrode distance in mm [X,Y]. 
                                 %check equipment specifications for details.
                                 %square mm ROI coverage, state x*y sizes, minimum is 9 square mm [3,3]
subj_info.ROIsize       = [11,11]; %try to keep ROIsize at/under [11,11] due to processing time (~12 hours for [11,11] 8x4 grid)
subj_info.stepSize      = 1; %set stepsize for the spacing of ROI-points (in mm).

subj_info.hemiVect.hemi = 'r'; %set hemisphere
subj_info.hemiVect.side = 'r'; %side of leads/cables: u (up), d (down), l (left) or r (right)

%enter subject code name
subj_info.subj_name     = 'name';

%Import freesurfer cortex file and convert to triangular surface model
%to use with neuralAct
if ~exist('hullcortex','var')
    %select niftii (.nii) cortex representation to convert into a triangular surface representation
    cortex              = gen_cortex_click_from_FreeSurfer(subj_info.subj_name, 0.5, [15 3], subj_info.hemiVect.hemi);
    cortexcoarser       = coarserModel(cortex, 10000);      %make the model coarser
    hullcortex          = hullModel(cortexcoarser);         %compute the convex hull
    subj_info.neuralAct = ['./Results/' subj_info.subj_name '_neuralAct.mat'];
    
    save(['./Results/' subj_info.subj_name '_neuralAct.mat'],'cortex','cortexcoarser','hullcortex');
end



%% 2 - Process the angiogram files:

subj_info.sfile      = ['./Results/' subj_info.subj_name '_balloon_11_03.nii']; %surface balloon of the brain, to project the vessels
subj_info.tfile      = '.../angiogram.nii'; %Angiogram T-file
subj_info.Tthreshold = 70;    %threshold voor T-statistic for angiogram T-file with vessel information
subj_info.VoxelDepth = 8;      %searchdepth for co-registration of vessels to surface

save(['./Results/' subj_info.subj_name '_subjectInfo.mat'],'subj_info');

%Calculate angioMap - Use this function to use the T-map and surface map
%files of the angiogram to calculate an angiomap which is coregistered to
%the cortex files from the NeuralAct package. T-threshold is equal to the
%lower boundary found for the angiogram via the software package MRIcron.
%Upper boundary can be adjusted from within the function, if necessary.
%N.B.: Save normAngio after running because it might take some time.
[angioMap,normAngio]  = calculateAngioMap(subj_info,subj_info.Tthreshold,subj_info.VoxelDepth,1);

save(['./Results/' subj_info.subj_name '_AngioMap.mat'],'angioMap','normAngio');



%% 3 - Calculate ROI and set Grid properties:

%Load the files previously created:
if ~exist('subj_info','var')
    load(['./Results/' subj_info.subj_name '_subjectInfo.mat']);
    load (subj_info.neuralAct);
    load (['./Results/' subj_info.subj_name '_AngioMap.mat']);
end

%set gridLayout
global dims;          %for gridLayout.m
global inputGrid      %for gridLayout.m
global numElec        %for gridLayout.m

dims      = subj_info.dims;
numElec   = dims(1)*dims(2);

% Define ROI and get tangent plane from ROI
global contOnCommand; %for getROI
contOnCommand = 1;
ROI_no        = 1;
[coords_ROI]  = getROI(cortex,hullcortex,contOnCommand,subj_info);


%Note the channel index, starting in the topleft corner.
%This is useful when the grids aren't numbered in ascending, or descending
%order but randomly. Standard inputGrid is set to 1:64 when clicking 'set layout' whilst leaving the fields empty.
%Meaning, top left corner is electrode 1 and bottom right corner is
%electrode e.g. 128.

%Set layout using interface.
% gridLayout; 
%or input indexes manually here:
inputGrid = 1:128;

%set dimensions of grid according to side of leads for the rest of the script.
%This way it accounts for rectangular grid
auxDims = dims;
if dims(1)~=dims(2)
    if (strcmp(subj_info.hemiVect.side,'l')+strcmp(subj_info.hemiVect.side,'r'))>0
        auxDims(1)=dims(2);
        auxDims(2)=dims(1);
    end
end
clear('ROI_no','contOnCommand','gridSpecs','ans'); %housekeeping

save(['./Results/' subj_info.subj_name '_GridInfo.mat'],'inputGrid','auxDims', 'dims', 'coords_ROI', 'numElec');



%% 4 - Project ROI and grids onto brain and calculate the model:
tic;
disp('Running. Please wait...');

%use neuralact to project ROI-points onto hullcortex
projectedROIpoints = [];
projectedROIpoints.electrodes = coords_ROI.ROI_Tangent;
[ projectedROIpoints ] = projectElectrodes(hullcortex,projectedROIpoints,25,0);

%Set rotation in degrees and derive amount of turns
%set this to 45 to run a quick test of code. Usually set to 1 degree during an actual run
rotation = 1;
turns = round((90/rotation)-1);

%create Grid per ROI point and project on cortex
ROI = createGrid( projectedROIpoints, rotation, turns, auxDims, subj_info, hullcortex);
elapsedTime1 = toc;
tic;

%Project the created grid onto the coarser model
ROI = projectToCoarser( ROI, cortexcoarser, turns);
elapsedTime2 = toc;
tic;

%Calculate model values
useAngio = exist('normAngio','var');
if useAngio == 1
    ROI = calculateModel( subj_info, ROI, cortex, normAngio);
else
    ROI = calculateModel( subj_info, ROI, cortex);
end
elapsedTime3 = toc;

% Save the result:
save(['./Results/' subj_info.subj_name ...
    '_ROI_' num2str(subj_info.ROIsize(1)) 'x' num2str(subj_info.ROIsize(1)) ...
    '_' date '.mat'],'ROI','-v7.3'); %save ROI struct after running



%% 5 - Run model correlation to HFB:

%housekeeping
clear('elapsedTime1','elapsedTime2','elapsedTime3','turns','rotation','projectedROIpoints')

%Load HFB_mean and index according to gridLayout.
load(subj_info.HFB_mean);
HFB_mean(subj_info.bad_channels)  = NaN; %set bad channels to NaN
HFB_mean                          = HFB_mean(subj_info.chansel);
%index it according to gridLayout
[forward_ind, inverse_ind]        = indexFuncLegacy(subj_info);
HFB_mean                          = HFB_mean(forward_ind); 

%Normalize HFB values (set between 0 and 1)
normHFB = rescale(HFB_mean);

%Calculate correlation between normalized HFB values and the predicted
%values for all Grids on all ROI points.
for i = 1:length(ROI)
    for ij = 1:length(ROI(1).coords)
        corrMatrix = corrcoef(ROI(i).coords(ij).weights, normHFB,'rows','complete');
        ROI(i).coords(ij).corr = corrMatrix(1,2);
    end
end

maxima = [];
%find indices of maximum correlation in matrix and set coordinates
for i = 1:length(ROI)
    maxima(i,:)= [ROI(i).coords.corr];
end

[a,b] = find(maxima(:,:) == max(maxima(:))); %find row and collumn maximum correlation
corrPred = ROI(a).coords(b).corr;

%display correlation of predicted
disp(['Best correlation = ' num2str(corrPred)]);



%% 6 - Display the the electrode prediction:

%display the cortex:
figure; 
colormap jet;
subaxis(5,4,[1:4, 5:8, 9:12],'SV',0,'SH',0.5);
ctmr_gauss_plot(cortex,[0 0 0],0); axis off; hold on;

%re-index predicted coordinates to match grid layout: 
coordsPred = ROI(a).coords(b).trielectrodes; 
coordsPred = coordsPred(inverse_ind,:);

%plot coordinates as spheres:
plotSpheres(coordsPred,[0.99,0.07,0.17]);
%or using the electrode label
% label_add(coordsPred);

%set view to correct hemisphere:
switch subj_info.hemiVect.hemi
    case 'l'
        loc_view(-90,0);
    case 'r'
        loc_view(90,0);
end

%save coordinates:
save(['./Results/' subj_info.subj_name '_projected_electrodes_gridloc_V1_' ...
    '_ROI_' num2str(subj_info.ROIsize(1)) 'x' num2str(subj_info.ROIsize(1)) ...
    '_' date '.mat'], 'coordsPred');



%% 7- Display the attenuation patterns:

%plot imagesc to compare normHFB values to predicted values from
%model
useAngio = exist('normAngio','var');
if useAngio == 1
    MvM     = ROI(a).coords(b).MvM;
end
McM     = ROI(a).coords(b).McM;
weights = ROI(a).coords(b).weights;

%Set imagescale plots
HFBSq = reshape(normHFB,auxDims(2),auxDims(1)); %HFB values per Elec
DepthSq = reshape(McM,auxDims(2),auxDims(1));       %Depth Values per Elec
if exist('MvM','var')==1
    AngioSq = reshape(MvM,auxDims(2),auxDims(1));   %Angio values per elec
end
CombinedSq = reshape(weights,auxDims(2),auxDims(1));%Depth+Angio per elec (if angio is not used, depth == combined).
roundCorr = round(corrPred,2); %rounded correlation for output figure

% correct imagescale plots to match view (L or R hemi)
switch subj_info.hemiVect.hemi
    case 'l'
        HFBSq = rot90((HFBSq),2);
        CombinedSq = rot90((CombinedSq),2);
        DepthSq = rot90((DepthSq),2);
        if exist('MvM','var')==1
            AngioSq = rot90((AngioSq),2);
        end
    case 'r'
        HFBSq = flipud(HFBSq);
        CombinedSq = flipud(CombinedSq);
        DepthSq = flipud(DepthSq);
        if exist('MvM','var')==1
            AngioSq = flipud(AngioSq);
        end
end

%set bad channels as background
bad_alpha                 = ones(size(HFBSq));
bad_alpha(isnan(HFBSq)) = 0; %setting transparency for bad channels

%plot HFB pattern
ax1 = subaxis(4,4,13,'Spacing',0);
imagesc(HFBSq, 'AlphaData', bad_alpha), colormap(ax1,parula),
title({'HFB pattern' '(leads = black)'}), axis square, axis off

%Create thick line to indicate lead side
switch subj_info.hemiVect.side
    case 'u'
        line([0.5,dims(1)+0.5],[0.5,0.5], 'Color','k','LineStyle','-','LineWidth',4);
    case 'd'
        line([0.5,dims(1)+0.5],[dims(2)+0.5,dims(2)+0.5], 'Color','k','LineStyle','-','LineWidth',4);
    case 'l'
        line([0.5,0.5],[0.5,dims(1)+0.5], 'Color','k','LineStyle','-','LineWidth',4);
    case 'r'
        line([dims(2)+0.5,dims(2)+0.5],[0.5,dims(1)+0.5], 'Color','k','LineStyle','-','LineWidth',4);
end

%plot combined model
ax2 = subaxis(4,4,14,'Spacing',0);
imagesc(CombinedSq), colormap(ax2,parula),
title({'Combined Model' ['(r = ',num2str(roundCorr),')']}), axis square, axis off

%plot depth model
ax3 = subaxis(4,4,15,'Spacing',0);
imagesc(DepthSq), colormap(ax3,parula),
title('Depth Model'), axis square, axis off

%plot angio model
if useAngio==1
    ax4 = subaxis(4,4,16,'Spacing',0);
    imagesc(AngioSq), colormap(ax4,parula),
    title('Angio Model'), axis square; axis off;
end

