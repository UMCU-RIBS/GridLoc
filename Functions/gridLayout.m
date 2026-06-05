classdef gridLayout < handle
    properties
        mainDialog
        controls
    end
    
    
    methods
        function obj = gridLayout
            global dims
            global numElec
            %create figure
            obj.mainDialog = figure('MenuBar','none','NumberTitle','off','Name','Set Grid Layout');
            %set gridpositions
            gridPos = [];
            i=1;
            for x = 1:dims(1)
                for y=dims(2):-1:1
                    gridPos(i,:)= [(x*50) (y*40)];
                    i = i+1;
                end
            end
            %create objects
            for k = 1:numElec
                obj.controls.box(k) = uicontrol( 'Parent', obj.mainDialog, 'Position',...
                    [gridPos(k,:) 30 30], 'Style', 'edit', 'Callback', @obj.boxCallback);
            end
            obj.controls.btn = uicontrol( 'Parent', obj.mainDialog, 'Position',...
                [400 5 150 30], 'String', 'Set Layout', 'Callback', @obj.btnCallback );
            
            obj.controls.txt = uicontrol('Style','text', 'FontSize', 12,...
                'Position',[160 5 150 20], 'String','LEADS ARE HERE','FontWeight','bold');
        end
    end
    
    methods % for callbacks
        function boxCallback( obj, ~, ~ )
            global inputGrid;
            global numElec;
            for k = 1:numElec
                inputGrid(k) = str2double(get(obj.controls.box(k), 'String'));
            end
        end
        function btnCallback(obj,~,~)
            global inputGrid;
            global numElec
            for k = 1:numElec
                if isnan(inputGrid(k))
                    errordlg('You must enter numeric values!','Invalid Input','modal')
                    uicontrol(obj.controls.box(k))
                    return
                end
            end
            errorInput = checkInput(inputGrid);
            
            if errorInput == 0
                close(gcf)
            end
        end
    end
end
