%% Gridloc Toolbox V1.2
%  25/05/2018
%
%  (C) Michael Leibbrand - michael.leib7@gmail.com
%  (C) Mariana P Branco  - m.pedrosobranco@umcutrecht.nl

warning off; %turn off to reduce clutter in command window. 
clear all; close all; clc;

%Paths of (support) files/scripts
addpath(genpath('...\GridLoc Toolbox'));
cd '...\GridLoc Toolbox\' 


%% 1- Create the subject info structure and set the location of the input files.

subj_info.gamma_mean    = 'gamma_mean.mat'; %Set location of file with the caclulated values of mean Gamma (i.e., high-frequency band) activity per electrode
subj_info.chansel       = [1:64]; %channel selection HD grid
subj_info.bad_channels  = [1, 8, 57, 64]; %set bad channels
subj_info.dims          = [8,8]; %amount of electrodes [X,Y]. X-direction is side where the leads are located.
subj_info.intElec       = [3,3]; %inter-electrode distance in mm [X,Y]. Check equipment specifications for details.
%square mm ROI coverage, state x*y sizes, minimum is 9 square mm [3,3]
subj_info.ROIsize       = [3,3]; %try to keep ROIsize at/under [11,11] due to processing time (~12 hours for [11,11] 8x4 grid)
subj_info.stepSize      = [1]; %set stepsize for the spacing of ROI-points (in mm). 

%Viewstruct and other NeuralAct viewBrain options 
%N.B.: Make sure you use the neuralAct package provided with THIS toolbox. It has been edited!
subj_info.hemiVect.hemi = 'l'; %set hemisphere
subj_info.hemiVect.side = 'u'; %side of leads/cables (u, d, l or r)
switch subj_info.hemiVect.hemi
    case 'r'
        subj_info.viewstruct.lightpos = [180,0,0]; 
        subj_info.what2view = {'brain'};
        subj_info.transp    = 1;
        subj_info.colix     = 0;
        subj_info.viewvect  = [90, 0];
    case 'l'
        subj_info.viewstruct.lightpos = [-180,0,0];
        subj_info.what2view = {'brain'};
        subj_info.viewvect  = [270, 0];
        subj_info.transp    = 1;
        subj_info.colix     = 0;
end

%Import freesurfer cortex file and convert to triangular surface model
%to use with neuralAct
if ~exist('hullcortex','var')
    %select niftii (.nii) cortex representation to convert into a triangular surface representation
    cortex = gen_cortex_click_from_FreeSurfer('SetSubjectName', 0.5, [15 3], subj_info.hemiVect.hemi); 
    cortexcoarser = coarserModel(cortex, 10000);      %make the model coarser
    hullcortex = hullModel(cortexcoarser);            %compute the convex hull
    save('neuralAct_SetSubjectName.mat','cortex','cortexcoarser','hullcortex')
    subj_info.neuralAct = 'neuralAct_SetSubjectName.mat';
end



%% 2 - Process the angiogram files:

subj_info.sfile      = 'Subject_surface.img'; %surface balloon of the brain, to project the vessels
subj_info.tfile      = 'Subject_masked_angio.img'; %Angiogram T-file
subj_info.Tthreshold = 311750; %threshold voor T-statistic for angiogram T-file with vessel information
subj_info.VoxelDepth = 8;      %searchdepth for co-registration of vessels to surface 
save('subjectInfo_SetSubjectName.mat','subj_info')

%Calculate angioMap - Use this function to use the T-map and surface map
%files of the angiogram to calculate an angiomap which is coregistered to
%the cortex files from the NeuralAct package. T-threshold is equal to the
%lower boundary found for the angiogram via the software package MRIcron.
%Upper boundary can be adjusted from within the function, if necessary. 
%N.B.: Save normAngio after running because it might take some time.
if ~exist('angioMap','var')
    [angioMap,normAngio]  = calculateAngioMap(subj_info,subj_info.Tthreshold,subj_info.VoxelDepth,1);
    save('AngioMap_SetSubjectName.mat','angioMap','normAngio')
end


%% 3 - Calculate ROI and set Grid properties:

%Load the files previosuly created:
if ~exist('subj_info','var')
    load('subjectInfo_SetSubjectName.mat')
    load (subj_info.neuralAct);
    load ('AngioMap_SetSubjectName.mat');
end

%set gridLayout
global dims;          %for gridLayout.m
global inputGrid      %for gridLayout.m
global numElec        %for gridLayout.m

dims      = subj_info.dims;%[8,8] means it's a squared 8x8=64 electrode grid, [4,8] means it's a 4 by 8=32 electrode grid. 
numElec   = dims(1)*dims(2);
inputGrid = 1:numElec; %set inputgrid to 1:64

% Define ROI and get tangent plane from ROI
global contOnCommand; %for getROI
contOnCommand = 1;
ROI_no        = 1;
[coords_ROI]  = getROI(cortex,hullcortex,contOnCommand,subj_info); 


%Note the channel index, starting in the topleft corner. 
%This is useful when the grids aren't numbered in ascending, or descending
%order but randomly. Standard inputGrid is set to 1:64 when clicking 'set layout' whilst leaving the fields empty.
%Meaning, top left corner is electrode 1 and bottom right corner is
%electrode 64. 

% gridLayout; %layout is set in inputGrid variable. 
inputgrid = 1:128;

%set dimensions of grid according to side of leads for the rest of the script.
%This way it accounts for rectangular grid  
auxDims = dims;
if dims(1)~=dims(2)
    if (strcmp(subj_info.hemiVect.side,'l')+strcmp(subj_info.hemiVect.side,'r'))>0
        auxDims(1)=dims(2);
        auxDims(2)=dims(1);
    end
end
clear('ROI_no','contOnCommand','gridSpecs','ans') %housekeeping

save('GridInfo_SetSubjectName.mat','inputGrid','auxDims', 'dims', 'coords_ROI', 'numElec');



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
elapsedTime3=toc;

% Save the result:
save('SubjectName_ROI_11x11.mat','ROI','-v7.3') %save ROI struct after running 



%% 5 - Run model correlation to Gamma and plot results:

%housekeeping
clear('elapsedTime1','elapsedTime2','elapsedTime3','turns','rotation','projectedROIpoints') 

%Load gamma_mean and index according to gridLayout. 
load (subj_info.gamma_mean)
gamma_mean(subj_info.bad_channels) = NaN; %set bad channels to NaN
gamma_mean = gamma_mean(subj_info.chansel);
ind = indexFuncLegacy(subj_info);
gamma_mean=gamma_mean(ind); %index it according to gridLayout, values are now in order of 1:64

%Normalize gamma values (set between 0 and 1)
normGamma = rescale(gamma_mean);
normGamma = normGamma';

%Calculate correlation between normalized gamma values and the predicted
%values for all Grids on all ROI points. 
for i = 1:length(ROI)
    for ij = 1:length(ROI(1).coords)
        corrMatrix = corrcoef(ROI(i).coords(ij).weights, normGamma,'rows','complete');
        ROI(i).coords(ij).corr = corrMatrix(1,2);
    end
end

maxima = [];
%find indices of maximum correlation in matrix and set coordinates
for i = 1:length(ROI)
    maxima(i,:)= [ROI(i).coords.corr];
end

[a,b]=find(maxima(:,:) == max(maxima(:))); %find row and collumn maximum correlation
corrPred = ROI(a).coords(b).corr;

%display correlation of predicted
disp(['correlation = ' num2str(corrPred)]);



%% 6 - Display the the electrode prediction:

% Display the cortex:
figure;
colormap jet;
subaxis(4,4,[1:4, 5:8, 9:12],'SV',0,'SH',0.5);
viewBrain(cortex,subj_info.hemiVect.hemi);
axis off;
hold on;

%display the predicted locations:
coordsPred = ROI(a).coords(b).trielectrodes;
plotSpheres(coordsPred,[0.99,0.07,0.17]); 



%% 7 - Display the attenuation patterns:

%plot imagesc to compare normGamma values to predicted values from
%model
if useAngio == 1
    MvM     = ROI(a).coords(b).MvM;
end
McM     = ROI(a).coords(b).McM;
weights = ROI(a).coords(b).weights;

%Set imagescale plots 
GammaSq = reshape(normGamma,auxDims(2),auxDims(1)); %Gamma values per Elec
DepthSq = reshape(McM,auxDims(2),auxDims(1));       %Depth Values per Elec
if exist('MvM','var')==1
    AngioSq = reshape(MvM,auxDims(2),auxDims(1));   %Angio values per elec
end
CombinedSq = reshape(weights,auxDims(2),auxDims(1));%Depth+Angio per elec (if angio is not used, depth == combined).
roundCorr = round(corrPred,2); %rounded correlation for output figure

% correct imagescale plots to match view (L or R hemi)
switch subj_info.hemiVect.hemi
    case 'l'
        GammaSq = rot90((GammaSq),2);
        CombinedSq = rot90((CombinedSq),2);
        DepthSq = rot90((DepthSq),2);
        if exist('MvM','var')==1
            AngioSq=rot90((AngioSq),2);
        end
    case 'r'
        GammaSq = flipud(GammaSq);
        CombinedSq = flipud(CombinedSq);
        DepthSq = flipud(DepthSq);
        if exist('MvM','var')==1
            AngioSq=flipud(AngioSq);
        end
end

%Plot imagesc matrices
if useAngio==1 
    ax1 = subaxis(4,4,13,'Spacing',0);
    imagesc(GammaSq), colormap(ax1,parula),
    title('Gamma Mean'), axis square, axis off
    
    ax2=subaxis(4,4,14,'Spacing',0); 
    imagesc(CombinedSq), colormap(ax2,parula),
    title(['Combined Model (r = ',num2str(roundCorr),')']), axis square, axis off
    
    ax3=subaxis(4,4,16,'Spacing',0); 
    imagesc(AngioSq), colormap(ax3,parula),
    title('Angio Model'), axis square; axis off
    
    ax4=subaxis(4,4,15,'Spacing',0); 
    imagesc(DepthSq), colormap(ax4,parula), 
    title('Depth Model'), axis square, axis off
else
    ax1 = subaxis(4,4,14,'Spacing',0); 
    imagesc(GammaSq), colormap(ax1,parula), 
    title('Gamma Mean'), axis square, axis off
    
    ax2=subaxis(4,4,15,'Spacing',0); 
    imagesc(DepthSq), colormap(ax2,parula),
    title(['Depth Model (r = ',num2str(roundCorr),')']),axis square, axis off
end


%% 8 - Save electrode coordinates in correct order:
coordsPred = coordsPred(ind,:);
save('intraopXXX_projected_electrodes_gridloc.mat', 'cortex', 'coordsPred');