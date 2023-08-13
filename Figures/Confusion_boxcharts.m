
function [] = Confusion_boxcharts(data)
%CONFUSION_BOXCHARTS Creates 2 boxcharts for every model
%   Data is either 'PP' or 'MS' and just specifies which to data to import
    
    % Creating an axis so ticks will appear between plots
%     t = nexttile(2);
%     x2 = linspace(0,18,19);
%     y2 = x2/9;
%     ax2 = axes(t);
%     plot(ax2,x2,y2, 'Color', '#FFFFFF')
%     xticks([3 7 11 15])

    % Creating boxcharts
    switch data
        case 'PP'
            %Import data
            Conservative = open("Results/conservative10000.mat").stats_mat;
            Default = open("Results/default10000.mat").stats_mat;
            %Liberal = open("Results/liberal10000.mat").stats_mat;
            Natural = open("Results/natural10000.mat").stats_mat;
            
            hold on

            %Plot the boxchart for each model
            for i = 1:3
                Y = nan(length(Conservative),6);
                P = nan(length(Conservative),6);
                if i == 1
                    Y(:,1) = Natural(:,4);
                    P(:,4) = Natural(:,5);
                    col = [0.8,0.475,0.655]; %col 4 is sens, col 5 is PPV
                elseif i == 2
                    Y(:,2) = Default(:,4);
                    P(:,5) = Default(:,5);
                    col = [0.835,0.369,0];
                else
                    Y(:,3) = Conservative(:,4);
                    P(:,6) = Conservative(:,5);
                    col = [0.337,0.706,0.914];
                end
                %Make the outliers invisible
                boxchart(Y, 'BoxFaceColor',col, 'MarkerColor','#FFFFFF', 'MarkerSize', 0.0001, 'WhiskerLineColor','b');
                boxchart(P, 'BoxFaceColor',col, 'MarkerColor','#FFFFFF', 'MarkerSize', 0.0001, 'WhiskerLineColor','r');
            end

            %Clean up the figure
            yyaxis left
            ylim([-0.025 1])
            ylabel('Sensitivity')
            ax = gca;
            ax.YColor = 'b';
            
            yyaxis right
            ax = gca;
            ax.YColor = 'r';
            ylabel('Positive Predictive Value')
            ylim([-0.025 1])
            xticklabels(["Natural" "Light" "Heavy" "Natural" "Light" "Heavy"])

      case 'MS'
            %Import data
            MSConservative = open("Results/MSconservative10000.mat").stats_mat;
            Conservative = open("Results/conservative10000.mat").stats_mat;
            Default = open("Results/MSdefault10000.mat").stats_mat;
            %Liberal = open("Results/liberal10000.mat").stats_mat;
            Natural = open("Results/MSnatural10000.mat").stats_mat;
            
            hold on

            %Plot the boxchart for each model
            for i = 1:4
                Y = nan(length(Conservative),8);
                P = nan(length(Conservative),8);
                if i == 1
                    Y(:,1) = Conservative(:,4);
                    P(:,5) = Conservative(:,5);
                    col = [0.337,0.706,0.914];
                elseif i == 2
                    Y(:,2) = Natural(:,4);
                    P(:,6) = Natural(:,5);
                    col = [0,0.620,0.451]; %col 4 is sens, col 5 is PPV
                elseif i == 3
                    Y(:,3) = Default(:,4);
                    P(:,7) = Default(:,5);
                    col = [0.902,0.624,0];
                else
                    Y(:,4) = MSConservative(:,4);
                    P(:,8) = MSConservative(:,5);
                    col = [0,0.447,0.698];
                end
                %Make the outliers invisible
                boxchart(Y, 'BoxFaceColor',col, 'MarkerColor','#FFFFFF', 'MarkerSize', 0.0001, 'WhiskerLineColor','b');
                boxchart(P, 'BoxFaceColor',col, 'MarkerColor','#FFFFFF', 'MarkerSize', 0.0001, 'WhiskerLineColor','r');
            end

            %Clean up the figure
            yyaxis left
            ylim([-0.025 1])
            ylabel('Sensitivity')
            ax = gca;
            ax.YColor = 'b';
            
            yyaxis right
            ax = gca;
            ax.YColor = 'r';
            ylabel('Positive Predictive Value')
            ylim([-0.025 1])
            xticklabels(["Heavy" "MS" "MS-light" "MS-heavy" "Heavy" "MS" "MS-light" "MS-heavy"])
            xtickangle(0)
    end
end


