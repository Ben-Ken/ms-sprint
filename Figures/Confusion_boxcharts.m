function [] = Confusion_boxcharts(data)
%CONFUSION_BOXCHARTS Creates 2 boxcharts for every model
%   Data is either 'PP' or 'MS' and just specifies which to data to import
    
    % Creating an axis so ticks will appear between plots
    t = tiledlayout(1,1);
    x2 = linspace(0,18,19);
    y2 = x2/9;
    ax2 = axes(t);
    plot(ax2,x2,y2, 'Color', '#FFFFFF')
    xticks([3 7 11 15])

    % Creating boxcharts
    switch data
        case 'PP'
            %Import data
            Conservative = open("Results/conservative10000.mat").stats_mat;
            Default = open("Results/default10000.mat").stats_mat;
            Liberal = open("Results/liberal10000.mat").stats_mat;
            Natural = open("Results/natural10000.mat").stats_mat;
            
            %Hide second axis and set ticks
            xticklabels(["Natural" "Default" "Liberal" "Conservative"])
            ax2.Box = 'off'
            ax2.XAxisLocation = ['bottom'];
            ax2.YColor = 'none';
            ax1 = axes(t);
            ax1.XColor = 'none';
            hold on

            %Plot the boxchart for each model
            for i = 1:4
                Y = nan(length(Conservative),8);
                if i == 1
                    Y(:,1) = Natural(:,4);
                    Y(:,2) = Natural(:,5);
                    col = [0.8,0.475,0.655]; %col 4 is sens, col 5 is PPV
                elseif i == 2
                    Y(:,3) = Default(:,4);
                    Y(:,4) = Default(:,5);
                    col = [0.835,0.369,0];
                elseif i == 3
                    Y(:,5) = Liberal(:,4);
                    Y(:,6) = Liberal(:,5);
                    col = [0.941,0.894,0.259];
                else
                    Y(:,7) = Conservative(:,4);
                    Y(:,8) = Conservative(:,5);
                    col = [0.337,0.706,0.914];
                end
                %Make the outliers invisible
                boxchart(Y, 'BoxFaceColor',col, 'MarkerColor','#FFFFFF', 'MarkerSize', 0.0001);
            end

            %Clean up the figure
            ylim([0 1])
            ax1.Color = 'none';
            ax2.XTickLabelRotation = 0
      case 'MS'
        %Import data
        Conservative = open("Results/conservative10000.mat").stats_mat;
        MSdef = open("Results/MSdefault10000.mat").stats_mat;
        MScons = open("Results/MSconservative10000.mat").stats_mat;
        MS = open("Results/MSnatural10000.mat").stats_mat;
    
        %Hide second axis and set ticks
        xticklabels(["Conservative" "MS" "MS-Default" "MS-Conservative"])
        ax2.Box = 'off'
        ax2.XAxisLocation = ['bottom'];
        ax2.YColor = 'none';
        ax1 = axes(t);
        ax1.XColor = 'none';
        hold on

        %Plot the boxchart for each model
        for i = 1:4
            Y = nan(length(Conservative),8);
            if i == 1
                Y(:,1) = Conservative(:,4);
                Y(:,2) = Conservative(:,5);
                col = [0.337,0.706,0.914]; %col 4 is sens, col 5 is PPV
            elseif i == 2
                Y(:,3) = MS(:,4);
                Y(:,4) = MS(:,5);
                col = [0,0.620,0.451];
            elseif i == 3
                Y(:,5) = MSdef(:,4);
                Y(:,6) = MSdef(:,5);
                col = [0.902,0.624,0];
            else
                Y(:,7) = MScons(:,4);
                Y(:,8) = MScons(:,5);
                col = [0,0.447,0.698];
            end
            % Make the outliers invisible
            boxchart(Y, 'BoxFaceColor',col, 'MarkerColor','#FFFFFF', 'MarkerSize', 0.0001);
        end

        %Clean up the figure
        ylim([0 1])
        ax1.Color = 'none';
        ax2.XTickLabelRotation = 0
    end
end


