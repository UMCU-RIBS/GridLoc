function [coords] = getROI (cortex, hullcortex, contOnCommand, subj_info)
%Specify cortex file (struct with .vert and .tri) on input. This function
%allows a user to select a point around which a 3D mesh (points with X,Y,Z) will be
%constructed. The size of this mesh is set in the subj_info structure
%(ROIsize). 
%%
global contOnCommand; %#ok<REDEF>
global ROI_no;
global dims;

side = subj_info.hemiVect.side;
hemi = subj_info.hemiVect.hemi;
gridSize = (dims-1).*subj_info.intElec;
intElec  = subj_info.intElec;
ROIsize  = subj_info.ROIsize;
stepSize = subj_info.stepSize;

%Set variables for the while loop
contOnCommand = 1;
ROI_no = 1;
run = 1;
coords = [];

%create figure, set GUI buttons
g = figure('Name', 'Get ROI', 'CloseRequestFcn', @figCloseRequest);
viewBrain(cortex,hemi);
h = uicontrol('Position',[20 60 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
h2 = uicontrol('Position',[20 20 200 40],'String','Close',...
              'Callback',{@callback,g,contOnCommand});

%Set the datacursormode for acquiring the selected point
datacursormode on;
dcm_obj = datacursormode(g);
fprintf(' Select the middlepoint of the ROI. When done click on continue. \n To close this figure, press close. \n');
uiwait(g);     %wait for input    
          
while contOnCommand ~= 0
try
        if run >= 2
            datacursormode on;
            dcm_obj = datacursormode(g);
            uiwait(g);
        end
        if contOnCommand == 0 %break loop if variable changes
            break 
        end
        
        %gettin the info from the cursor
        f = getCursorInfo(dcm_obj);
        delete(findall(g,'Type','hggroup','HandleVisibility','off'))
        datacursormode off;
        disp('Plotting the Region of Interest and Covered Area.');

        c1 = f.Position;     %Setting c1 as the middlepoint of the grid
        
        %project the middlepoint onto the hullcortex, prevents ROI from
        %being plotted inside the cortex model
        s = [];
        s.electrodes = c1;
        [s]=projectElectrodes(hullcortex,s,40);
        c=s.trielectrodes;
        N=s.normal;
        
        %calculate coordinates ROI and coverage of used ROI
        coords_ROI = [];
        coords_ROI = calcCoords(c,stepSize,ROIsize,1);
       
        %setting side of leads
        if hemi == 'l'
            switch side
                case 'r'
                    leads = [1:floor(ROIsize(2)/stepSize)]; %right row of ROI points
                case 'l'
                    leads = [(floor(ROIsize(1)/stepSize)*floor(ROIsize(2)/stepSize))-floor(ROIsize(2)/stepSize)+1:floor(ROIsize(1)/stepSize)*floor(ROIsize(2)/stepSize)]; %left row of ROI points
                case 'u'
                    leads = [floor(ROIsize(2)/stepSize):floor(ROIsize(2)/stepSize):floor(ROIsize(1)/stepSize)*floor(ROIsize(2)/stepSize)];
                case 'd'
                    leads = [1:floor(ROIsize(2)/stepSize):floor(ROIsize(1)/stepSize)*floor(ROIsize(2)/stepSize)];
            end
        else
            switch side
                case 'r'
                    leads = [(floor(ROIsize(1)/stepSize)*floor(ROIsize(2)/stepSize))-floor(ROIsize(2)/stepSize)+1:floor(ROIsize(1)/stepSize)*floor(ROIsize(2)/stepSize)]; %left row of ROI points
                case 'l'
                    leads = [1:floor(ROIsize(2)/stepSize)]; %right row of ROI points
                case 'u'
                    leads = [floor(ROIsize(2)/stepSize):floor(ROIsize(2)/stepSize):floor(ROIsize(1)/stepSize)*floor(ROIsize(2)/stepSize)];
                case 'd'
                    leads = [1:floor(ROIsize(2)/stepSize):floor(ROIsize(1)/stepSize)*floor(ROIsize(2)/stepSize)];
            end
        end

        %Set dimensions of grids (grid starts at 0 instead of one, so gridsize +1 =
        %dims)
        altDims  = floor(ROIsize/stepSize);
        lngth = length(coords_ROI);
        
        %Calculate tangent planes
        [coords(ROI_no).ROI_Tangent] = calcTangent(hullcortex, c,coords_ROI,altDims,lngth,hemi);
        
        %plotting 
        hold on;
        plot_ROI = plot3(coords(ROI_no).ROI_Tangent(:,1),coords(ROI_no).ROI_Tangent(:,2),coords(ROI_no).ROI_Tangent(:,3),'b*');
        %plot one side of ROI yellow to indicate leads
        plot_Leads=plot3(coords(ROI_no).ROI_Tangent(leads,1),coords(ROI_no).ROI_Tangent(leads,2),coords(ROI_no).ROI_Tangent(leads,3),'y*');
        covRadius=(0.5*sqrt(gridSize(1)^2 + gridSize(2)^2))+(0.5*sqrt(ROIsize(1)^2+ROIsize(2)^2));
        plotCircle3D(c,N,covRadius);
        
        
        %Creating different UI buttons with different functions
        h1 = uicontrol('Position',[20 20 200 40],'String','Close',...
                      'Callback',{@callback,g,contOnCommand});
        j = uicontrol('Position',[20 60 200 40],'String','Save',...
                      'Callback',{@callback2,g,ROI_no}');
        k = uicontrol('Position',[20 100 200 40],'String','Reset',...
                      'Callback','uiresume(gcbf); contOnCommand = 1;');

        fprintf(' The Region Of Interest is shown in BLUE. \n The RED circle displays the Covered area. \n The yellow colored side displays the side of the leads (cables). \n');
        uiwait(g); 

        %Deleting the previous ROI and coverage
        u =findobj('type','line');
        delete(u(1));delete(u(2));; delete(h1);delete(j);delete(k);delete(u(3))

        run = run+1; %this updates the run so the ROI's are saved in different substructs
   catch
       warning('Please select a valid point on the brain!')
       warning('Please select a valid point on the brain!')
       warning('Please select a valid point on the brain!')
       run = run+1;
end
end
close(gcf)
end

function figCloseRequest(srv,event,~)
global contOnCommand;
    delete(gcf)
    clc;
    contOnCommand = 0;    
end

function callback(srv,event,g,contOnCommand)
global contOnCommand;
    contOnCommand = 0;
    uiresume(g);
end

function callback2(srv,event,g,ROI_no)
global ROI_no;
    ROI_no = ROI_no + 1;
    contOnCommand = 1;
    uiresume(g);    
end