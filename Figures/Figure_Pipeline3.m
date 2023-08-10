%% Figure Pipeline
%% ===================================

%% Stats
figure('Position', [250 250 600 600])
color_error('PP')
title(['Stats of PP models'])
saveas(1,'Figures/PPstats', 'png')

figure('Position', [250 250 600 600])
error_plot('PP', 1)
title('Exponent error across PP models')
saveas(2,'Figures/PPerror', 'png')

figure('Position', [250 250 600 600])
color_error('MS')
title(['Stats of MS models'])
saveas(3,'Figures/MSstats', 'png')

figure('Position', [250 250 600 600])
error_plot('MS', 1)
title('Exponent error across MS models')
saveas(4,'Figures/MSerror', 'png')

close all


%% Sensitivities and error values
figure('Position', [250 250 600 600])
Sensitivity_Matrix('conservative', 1)
title('Peak probability by # peaks for PP')
saveas(1,'Figures/PPmatrix', 'png')

figure('Position', [250 250 600 600])
Sensitivity_Matrix('MSnatural', 2)
title('Peak probability by # peaks for MS')
saveas(2,'Figures/MSmatrix', 'png')

figure('Position', [250 250 600 600])
Periodic_Parameters('conservative', 'rgb', 'cf', 'PP')
title('Peak probability by CF for PP')
xlim([3,35]);
saveas(1,'Figures/PPcf_rgb2', 'png')

figure('Position', [250 250 600 600])
Periodic_Parameters('MSnatural', 'rgb', 'cf', 'MS')
title('Peak probability by CF for MS')
ylabel('')
xlim([3,35]);
saveas(2,'Figures/MScf_rgb2', 'png')

figure('Position', [250 250 600 600])
Periodic_Parameters('conservative', 'rgb', 'amp', 'PP')
title('Peak probability by AMP for PP')
ylabel('p(peak)')
saveas(5,'Figures/PPamp', 'png')

figure('Position', [250 250 600 600])
Periodic_Parameters('MSnatural', 'rgb', 'amp', 'MS')
title('Peak probability by AMP for MS')
saveas(6,'Figures/MSamp', 'png')

figure('Position', [250 250 600 600])
Error_Exponent_bypk()
title('Aperiodic exponent error by number of peaks')
ylabel('|error|')
xlabel('Aperiodic exponent (1/Hz)')
saveas(1,'Figures/Err_exp2', 'png')

figure('Position', [250 250 600 600])
Error_Frequency()
title('Absolute error by frequency')
ylabel('|error|')
xlabel('Frequency (Hz)')
saveas(1,'Figures/Err_freq2', 'png')

close all

% Supplemental
compare all the different model selection modes
sensitivity as a function of standard deviation

figure('Position', [250 250 600 600])
Stat_Frequency('conservative',1, 'PP')
title('Sens and PPV by frequency for PP')
xlabel('Frequency (Hz)')
saveas(1,'Figures/Supplemental/PPcave', 'png')

figure('Position', [250 250 600 600])
Stat_Frequency('MSnatural',2, 'MS')
title('Sens and PPV by frequency for MS')
xlabel('Frequency (Hz)')
saveas(2,'Figures/Supplemental/MScave', 'png')

close all

%% Support functions
function [] = color_error(data)
    t = tiledlayout(1,1);
    x2 = linspace(0,18,19);
    y2 = x2/9;
    ax2 = axes(t);
    plot(ax2,x2,y2, 'Color', '#FFFFFF')
    xticks([3 7 11 15])
    switch data
        case 'PP'
            Conservative = open("Results/conservative10000.mat").stats_mat;
            Default = open("Results/default10000.mat").stats_mat;
            Liberal = open("Results/liberal10000.mat").stats_mat;
            Natural = open("Results/natural10000.mat").stats_mat;
            
            xticklabels(["Natural" "Default" "Liberal" "Conservative"])
            ax2.Box = 'off'
            ax2.XAxisLocation = ['bottom'];
            ax2.YColor = 'none';
            ax1 = axes(t);
            ax1.XColor = 'none';
            hold on
            for i = 1:4
                Y = nan(length(Conservative),8);
                if i == 1
                    Y(:,1) = Natural(:,4);
                    Y(:,2) = Natural(:,5);
                    col = [0.8,0.475,0.655];
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
                % Make a bar plot
                boxchart(Y, 'BoxFaceColor',col, 'MarkerColor','#FFFFFF', 'MarkerSize', 0.0001);
            end
            ylim([0 1])
            ax1.Color = 'none';
            ax2.XTickLabelRotation = 0
      case 'MS'
        Conservative = open("Results/conservative10000.mat").stats_mat;
        MSdef = open("Results/MSdefault10000.mat").stats_mat;
        MScons = open("Results/MSconservative10000.mat").stats_mat;
        MS = open("Results/MSnatural10000.mat").stats_mat;
    
        xticklabels(["Conservative" "MS" "MS-Default" "MS-Conservative"])
        ax2.Box = 'off'
        ax2.XAxisLocation = ['bottom'];
        ax2.YColor = 'none';
        ax1 = axes(t);
        ax1.XColor = 'none';
        hold on
        for i = 1:4
            Y = nan(length(Conservative),8);
            if i == 1
                Y(:,1) = Conservative(:,4);
                Y(:,2) = Conservative(:,5);
                col = [0.337,0.706,0.914];
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
            % Make a bar plot
            boxchart(Y, 'BoxFaceColor',col, 'MarkerColor','#FFFFFF', 'MarkerSize', 0.0001);
        end
        ylim([0 1])
        ax1.Color = 'none';
        ax2.XTickLabelRotation = 0
    end
end

function [] = error_plot(data, p)
    % import the mAE of the models we want to plot
    cons_mAE = open('Results/conservative10000.mat').mAE_cell;
    if data == 'PP'
        lib_mAE = open('Results/liberal10000.mat').mAE_cell;
        nat_mAE = open('Results/natural10000.mat').mAE_cell;
        def_mAE = open('Results/default10000.mat').mAE_cell;
    else
        MS_mAE = open('Results/MSnatural10000.mat').mAE_cell;
        MSdef_mAE = open('Results/MSdefault10000.mat').mAE_cell;
        MScons_mAE = open('Results/MSconservative10000.mat').mAE_cell;
    end

    T = {'Exponent', 'Offset', 'Center frequency', 'Amplitude','Standard Deviation'};
    
    
    %figure('Position', [100 100 750 300])
    switch data
        case 'PP'
            mAE_cell = [nat_mAE(:,p) def_mAE(:,p) lib_mAE(:,p) cons_mAE(:,p)];
        case 'MS'
            mAE_cell = [cons_mAE(:,p) MS_mAE(:,p) MSdef_mAE(:,p) MScons_mAE(:,p)];
    end
    % % Panel B: error by spectral parameter
    swatch = [52 72 148]./255;
    violin(mAE_cell,data,'facecolor',swatch,'facealpha',1);
    yticks([-4 -3 -2 -1 0 1])
    yticklabels({'10^{-4}','10^{-3}','10^{-2}','0.1','1','10'})
    ylim([-3.2 1])
    ylabel('|error|')
    xticks(1:6)
    ax = gca(); 
    switch data
        case 'PP'
            row1 = {'Natural','Default','Liberal','Conservative'};
        case 'MS'
            row1 = {'Conservative','MS','MS-def','MS-cons'};  
    end
    %row2 = {'units','more units'};
    labelArray = [row1]; 
    %labelArray = [row1; row2]; 
    %tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
    xticklabels(labelArray);
end

function[h,L,MX,MED,bw]=violin(Y,col,varargin)
% This function is modified from the following source:
% Hoffmann H, 2015: violin.m - Simple violin plot using matlab default kernel
% density estimation. INRES (University of Bonn), Katzenburgweg 5, 53115 Germany.
% hhoffmann@uni-bonn.de
% INPUT
%
% Y:     Data to be plotted, being either
%        a) n x m matrix. A 'violin' is plotted for each column m, OR
%        b) 1 x m Cellarry with elements being numerical colums of nx1 length.
%
% varargin:
% xlabel:    xlabel. Set either [] or in the form {'txt1','txt2','txt3',...}
% facecolor: FaceColor. (default [1 0.5 0]); Specify abbrev. or m x 3 matrix (e.g. [1 0 0])
% edgecolor: LineColor. (default 'k'); Specify abbrev. (e.g. 'k' for black); set either [],'' or 'none' if the mean should not be plotted
% facealpha: Alpha value (transparency). default: 0.5
% mc:        Color of the bars indicating the mean. (default 'k'); set either [],'' or 'none' if the mean should not be plotted
% medc:      Color of the bars indicating the median. (default 'r'); set either [],'' or 'none' if the mean should not be plotted
% bw:        Kernel bandwidth. (default []); prescribe if wanted as follows:
%            a) if bw is a single number, bw will be applied to all
%            columns or cells
%            b) if bw is an array of 1xm or mx1, bw(i) will be applied to cell or column (i).
%            c) if bw is empty (default []), the optimal bandwidth for
%            gaussian kernel is used (see Matlab documentation for
%            ksdensity()
%
% OUTPUT
%
% h:     figure handle
% L:     Legend handle
% MX:    Means of groups
% MED:   Medians of groups
% bw:    bandwidth of kernel

%defaults:
%_____________________
xL=[];
fc=[205 205 205]./255;
lc='k';
mc=[];%'k';
medc='r';
b=[]; %bandwidth
plotlegend=0;
plotmean=0;
plotmedian=1;
x = [];
%_____________________
%convert single columns to cells:
if iscell(Y)==0
    Y = num2cell(Y,1);
end
%get additional input parameters (varargin)
if isempty(find(strcmp(varargin,'xlabel')))==0
    xL = varargin{find(strcmp(varargin,'xlabel'))+1};
end
if isempty(find(strcmp(varargin,'facecolor')))==0
    fc = varargin{find(strcmp(varargin,'facecolor'))+1};
end
if isempty(find(strcmp(varargin,'edgecolor')))==0
    lc = varargin{find(strcmp(varargin,'edgecolor'))+1};
end
if isempty(find(strcmp(varargin,'facealpha')))==0
    alp = varargin{find(strcmp(varargin,'facealpha'))+1};
end
if isempty(find(strcmp(varargin,'mc')))==0
    if isempty(varargin{find(strcmp(varargin,'mc'))+1})==0
        mc = varargin{find(strcmp(varargin,'mc'))+1};
        plotmean = 1;
    else
        plotmean = 0;
    end
end
if isempty(find(strcmp(varargin,'medc')))==0
    if isempty(varargin{find(strcmp(varargin,'medc'))+1})==0
        medc = varargin{find(strcmp(varargin,'medc'))+1};
        plotmedian = 1;
    else
        plotmedian = 0;
    end
end
if isempty(find(strcmp(varargin,'bw')))==0
    b = varargin{find(strcmp(varargin,'bw'))+1};
    if length(b)==1
        disp(['same bandwidth bw = ',num2str(b),' used for all cols'])
        b=repmat(b,size(Y,2),1);
    elseif length(b)~=size(Y,2)
        warning('length(b)~=size(Y,2)')
        error('please provide only one bandwidth or an array of b with same length as columns in the data set')
    end
end
if isempty(find(strcmp(varargin,'plotlegend')))==0
    plotlegend = varargin{find(strcmp(varargin,'plotlegend'))+1};
end
if isempty(find(strcmp(varargin,'x')))==0
    x = varargin{find(strcmp(varargin,'x'))+1};
end
%%
if size(fc,1)==1
    %fc=repmat(fc,size(Y,2),1);
    switch col
        case 'PP'
            fc = [0.8,0.475,0.655;0.835,0.369,0;0.941,0.894,0.259;0.337,0.706,0.914];
        case 'MS'
            fc = [0.337,0.706,0.914;0,0.620,0.451;0.902,0.624,0;0,0.447,0.698];
    end
end
%% Calculate the kernel density
i=1;
%fa = linspace(1,1./size(Y,2),size(Y,2));
fa = repmat(0.40,size(Y,2),1)'
for i=1:size(Y,2)
    
    if isempty(b)==0
        [f, u, bb]=ksdensity(Y{i},linspace(min(Y{i}),max(Y{i}),800),'bandwidth',b(i));
    elseif isempty(b)
        [f, u, bb]=ksdensity(Y{i},linspace(min(Y{i}),max(Y{i}),800));
        %so for some reason, after running thru 5, on the second run,
        %number 5 has an empty vector instead of values...
    end
    
    f=f/max(f)*0.3; %normalize
    F(:,i)=f;
    U(:,i)=u;
    MED(:,i)=nanmedian(Y{i});
    MX(:,i)=nanmean(Y{i});
    bw(:,i)=bb;
    prc25(:,i) = prctile(Y{i},25);
    prc75(:,i) = prctile(Y{i},75);
    IQR(:,i) = diff([prc25(:,i) prc75(:,i)]);
    
end
%%
%-------------------------------------------------------------------------
% Put the figure automatically on a second monitor
% mp = get(0, 'MonitorPositions');
% set(gcf,'Color','w','Position',[mp(end,1)+50 mp(end,2)+50 800 600])
%-------------------------------------------------------------------------
%Check x-value options
if isempty(x)
    x = zeros(size(Y,2));
    setX = 0;
else
    setX = 1;
    if isempty(xL)==0
        disp('_________________________________________________________________')
        warning('Function is not designed for x-axis specification with string label')
        warning('when providing x, xlabel can be set later anyway')
        error('please provide either x or xlabel. not both.')
    end
end
%% Plot the violins
i=1;
for i=i:size(Y,2)
    if isempty(lc) == 1
        if setX == 0
            h(i)=fill([F(:,i)+i;flipud(i-F(:,i))],[U(:,i);flipud(U(:,i))],fc(i,:),'FaceAlpha',fa(i),'EdgeColor','none');
        else
            h(i)=fill([F(:,i)+x(i);flipud(x(i)-F(:,i))],[U(:,i);flipud(U(:,i))],fc(i,:),'FaceAlpha',fa(i),'EdgeColor','none');
        end
    else
        if setX == 0
            h(i)=fill([F(:,i)+i;flipud(i-F(:,i))],[U(:,i);flipud(U(:,i))],fc(i,:),'FaceAlpha',fa(i),'EdgeColor',lc);
        else
            h(i)=fill([F(:,i)+x(i);flipud(x(i)-F(:,i))],[U(:,i);flipud(U(:,i))],fc(i,:),'FaceAlpha',fa(i),'EdgeColor',lc);
        end
    end
    hold on
    if setX == 0
        if plotmean == 1
            p(1)=plot([interp1(U(:,i),F(:,i)+i,MX(:,i)), interp1(flipud(U(:,i)),flipud(i-F(:,i)),MX(:,i)) ],[MX(:,i) MX(:,i)],mc,'LineWidth',2);
        end
        if plotmedian == 1
            plot(ones(2,1).*mean([interp1(U(:,i),F(:,i)+i,MED(:,i)), interp1(flipud(U(:,i)),flipud(i-F(:,i)),MED(:,i))]),...
                [max([min(U(:,i)) prc25(:,i)-1.5.*IQR(:,i)]) min([max(U(:,i)) prc75(:,i)+1.5.*IQR(:,i)])],'Color',[0.2 0.2 0.2])
            plot(ones(2,1).*mean([interp1(U(:,i),F(:,i)+i,MED(:,i)), interp1(flipud(U(:,i)),flipud(i-F(:,i)),MED(:,i))]),...
                [prc25(:,i) prc75(:,i)],'Color',[0.2 0.2 0.2],'LineWidth',7)
            scatter(mean([interp1(U(:,i),F(:,i)+i,MED(:,i)), interp1(flipud(U(:,i)),flipud(i-F(:,i)),MED(:,i))]),MED(:,i),...
                abs(diff([interp1(U(:,i),F(:,i)+i,MED(:,i)), interp1(flipud(U(:,i)),flipud(i-F(:,i)),MED(:,i))])).*40,[1 1 1],'Filled');
            plot([F(:,i)+i;flipud(i-F(:,i))],[U(:,i);flipud(U(:,i))],'k')            
        end
    elseif setX == 1
        if plotmean == 1
            p(1)=plot([interp1(U(:,i),F(:,i)+i,MX(:,i))+x(i)-i, interp1(flipud(U(:,i)),flipud(i-F(:,i)),MX(:,i))+x(i)-i],[MX(:,i) MX(:,i)],mc,'LineWidth',2);
        end
        if plotmedian == 1
            p(2)=plot([interp1(U(:,i),F(:,i)+i,MED(:,i))+x(i)-i, interp1(flipud(U(:,i)),flipud(i-F(:,i)),MED(:,i))+x(i)-i],[MED(:,i) MED(:,i)],medc,'LineWidth',2);
        end
    end
end
%% Add legend if requested
if plotlegend==1 & plotmean==1 | plotlegend==1 & plotmedian==1
    
    if plotmean==1 & plotmedian==1
        L=legend([p(1) p(2)],'Mean','Median');
    elseif plotmean==0 & plotmedian==1
        L=legend([p(2)],'Median');
    elseif plotmean==1 & plotmedian==0
        L=legend([p(1)],'Mean');
    end
    
    set(L,'box','off','FontSize',14)
else
    L=[];
end
%% Set axis
if setX == 0
    axis([0.5 size(Y,2)+0.5, min(U(:)) max(U(:))]);
elseif setX == 1
    axis([min(x)-0.05*range(x) max(x)+0.05*range(x), min(U(:)) max(U(:))]);
end
%% Set x-labels
xL2={''};
i=1;
for i=1:size(xL,2)
    xL2=[xL2,xL{i},{''}];
end
% set(gca,'TickLength',[0 0],'FontSize',12)
box on
if isempty(xL)==0
    set(gca,'XtickLabel',xL2)
end
%-------------------------------------------------------------------------
end %of function

function [] = Sensitivity_Matrix(filename, col)
    load(['Results/' filename '10000.mat'])
    ts(:) = linspace(1.5, 58.5, 115);
    npeak_bic_mat = zeros(5,7);
    det_v_sim_mat = zeros(length(rstmp), length(ts), 2);
    
    for sim = 1:length(rstmp)
        for time = 1:length(ts)
            if isempty(rstmp(sim).pk_det)
                det_v_sim_mat(sim,time,1) = sum([rstmp(sim).resid_peaks(:).time] == ts(time)); %the detected peaks
            else
                det_v_sim_mat(sim,time,1) = sum([rstmp(sim).resid_peaks(:).time] == ts(time)) + sum(rstmp(sim).pk_det(:,time), 'omitnan'); %the detected peaks
            end
            if ~isempty(vstmp(sim).pk_exist)
                det_v_sim_mat(sim,time,2) = sum(vstmp(sim).pk_exist(:,time)); %the simulated peaks
            else
                det_v_sim_mat(sim,time,2) = 0;
            end
        end
    end
    
    for simulated = 1:size(npeak_bic_mat,1)
        for detected = 1:size(npeak_bic_mat,2)
            npeak_bic_mat(simulated,detected) = sum(logical(logical(det_v_sim_mat(:,:,1) == (detected - 1) & logical(det_v_sim_mat(:,:,2) == (simulated - 1)))),'all');
        end
    end
    
    data = nan(5,5);
    for row = 1:size(npeak_bic_mat,1);
        data(row,:) = npeak_bic_mat(row,1:5)./sum(npeak_bic_mat(row,:))*1000;
    end
    % c3 = [21 102 5]./255;
    c3 = {[0,0.447,0.698], [0.835,0.369,0]};
    c3 = c3{col};
    hold on
    xs = [0 1 2 3 4];
    ys = [0 1 2 3 4];
    imagesc(xs,ys, data')
    for x = 1:5
      for y = 1:5
        text(xs(x),ys(y),num2str(round(data(x,y)/1000,3)*100),'HorizontalAlignment','center','VerticalAlignment','middle') 
      end
    end
    yticks(0:4)
    xticks(0:4)
    xlim([-0.5 4.5])
    ylim([-0.5 4.5])
    ylabel('Number of detected peaks')
    xlabel('Number of simulated peaks')
    colormap([linspace(1,c3(1),1000)' linspace(1,c3(2),1000)' linspace(1,c3(3),1000)']);
    caxis([0 1000])
    c = colorbar('Ticks',[0,250,500,750,1000],...
         'TickLabels',{'0','25','50','75','100'});
    c.Label.String = 'Density (%)';
    %set(gca,'FontSize',14)
end

function [] = Periodic_Parameters(filename, cf_mode, parameter, type)
    load(['Results/' filename '10000.mat'])
    ts(:) = linspace(1.5, 58.5, 115);
    
    c = 1;
    for s = 1:length(vstmp)
        for i = 1:size(vstmp(s).pk_cf,1)
        c = c+1;
        err_timefreq(c,1) = ostmp(s).pk_times(i,2); % duration, in s
        err_timefreq(c,2) = ostmp(s).pk_cf(i); % centre frequency
        err_timefreq(c,3) = ostmp(s).pk_sd(i); % sd
        err_timefreq(c,4) = ostmp(s).pk_amp(i); % amplitude
        %err_timefreq(c,5) = sum(~isnan(rstmp(s).pk_cf(i,:)))./sum(~isnan(vstmp(s).pk_cf(i,:))); % peak detection rate
        err_timefreq(c,5) = sum(~isnan(rstmp(s).pk_cf(i,~isnan(vstmp(s).pk_cf(i,:)))))./sum(~isnan(vstmp(s).pk_cf(i,:))); % p.d.r fig. 3
        end
    end
    
    %by center frequency
    switch parameter
        case 'cf'
            switch cf_mode  %alternatives are 'bw' (black and white) or 'rgb' (color)
                case 'bw'
                    switch type
                        case 'PP'
                            swatch = [0,0.447,0.698];
                        case'MS'
                            swatch = [0.835,0.369,0];
                    end
                    scatter(err_timefreq(:,2),err_timefreq(:,5),6,swatch,'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
                case 'rgb'
                    pars = [0.85 4];
                    q = 0;
                    swatch = [108 156 167;236 189 78;198 97 43;162 55 27;91 61 28;86 34 98]./255;
                    hold on
%                     scatter((err_timefreq((err_timefreq(:,2)<8),2)),err_timefreq((err_timefreq(:,2)<8),5),6,swatch(1,:),'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
%                     scatter((err_timefreq((8<=err_timefreq(:,2) & err_timefreq(:,2)<13),2)),err_timefreq((8<=err_timefreq(:,2) & err_timefreq(:,2)<13),5),6,swatch(2,:),'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
%                     scatter((err_timefreq((13<=err_timefreq(:,2) & err_timefreq(:,2)<18),2)),err_timefreq((13<=err_timefreq(:,2) & err_timefreq(:,2)<18),5),6,swatch(3,:),'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
%                     scatter((err_timefreq((18<=err_timefreq(:,2)),2)),err_timefreq(18<=(err_timefreq(:,2)),5),6,swatch(4,:),'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
                    scatter(err_timefreq(err_timefreq(:,2)<8 & err_timefreq(:,1)>q,2),err_timefreq(err_timefreq(:,2)<8 & err_timefreq(:,1)>q,5),6,swatch(1,:),'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3)
                    scatter(err_timefreq(abs(err_timefreq(:,2)-10.5)<2.5 & err_timefreq(:,1)>q,2),err_timefreq(abs(err_timefreq(:,2)-10.5)<2.5 & err_timefreq(:,1)>q,5),6,swatch(2,:),'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3)
                    scatter(err_timefreq(abs(err_timefreq(:,2)-15.5)<2.5 & err_timefreq(:,1)>q,2),err_timefreq(abs(err_timefreq(:,2)-15.5)<2.5 & err_timefreq(:,1)>q,5),6,swatch(3,:),'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3)
                    scatter(err_timefreq(err_timefreq(:,2)>18 & err_timefreq(:,1)>q,2),err_timefreq(err_timefreq(:,2)>18 & err_timefreq(:,1)>q,5),6,[0.5 0.5 0.5],'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3)
                    hold off
            end
            xlabel('Center Frequency (Hz)')
            ylabel('p(peak)')
        case 'amp'
            % Peak detection rate by peak amplitude
            switch type
                case 'PP'
                    swatch = [0,0.447,0.698];
                case'MS'
                    swatch = [0.835,0.369,0];
            end
            scatter(err_timefreq(:,4),err_timefreq(:,5),6,swatch,'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
            xlabel('Amplitude (a.u.)')
            %ylabel('p(peak)')
            xlim([0.6 1.6])
    
        case 'sd'
            switch type
                case 'PP'
                    swatch = [0,0.447,0.698];
                case'MS'
                    swatch = [0.835,0.369,0];
            end
            scatter(err_timefreq(:,3),err_timefreq(:,5),6,swatch,'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
            xlabel('Standard Deviation (Hz)')
            %ylabel('p(peak)')
            xlim([1 2])
    end
end

function [] = Stat_Frequency(filename, i, type)
    load(['Results/' filename '10000.mat'])
    freqs_mat = nan(40, 3);
    buckets = linspace(1, 40, 40);
    
    for chan = 1:length(rstmp)
        for freq = 1:length(buckets)
            %true positives
            if isnan(freqs_mat(freq, 1))
                freqs_mat(freq,1) = sum(buckets(freq) == round(vstmp(chan).pk_cf,0).*rstmp(chan).pk_det, 'all');
            else
                freqs_mat(freq, 1) = freqs_mat(freq, 1) + sum(buckets(freq) == round(vstmp(chan).pk_cf,0).*rstmp(chan).pk_det, 'all');
            end
            %false positives
            if isnan(freqs_mat(freq,2))
                freqs_mat(freq,2) = sum(round([rstmp(chan).resid_peaks(:).center_frequency],0) == buckets(freq), 'all');
            else
                freqs_mat(freq,2) = freqs_mat(freq,2) + sum(round([rstmp(chan).resid_peaks(:).center_frequency],0) == buckets(freq), 'all');
            end
            %false negatives
            if isnan(freqs_mat(freq, 3))
                freqs_mat(freq,3) = sum(buckets(freq) == round(vstmp(chan).pk_cf,0).*(1-rstmp(chan).pk_det), 'all');
            else
                freqs_mat(freq,3) = freqs_mat(freq, 3) + sum(buckets(freq) == round(vstmp(chan).pk_cf,0).*(1-rstmp(chan).pk_det), 'all');
            end        
        end
    end
    
    sens(:) = [freqs_mat(:,1)./(freqs_mat(:,1) + freqs_mat(:,3))];
    PPV(:) = [freqs_mat(:,1)./(freqs_mat(:,1) + freqs_mat(:,2))];
    
    hold on
    
    bottom = sens;    
    top = PPV;
    
    yyaxis left
    a = bar(linspace(1,40,40), bottom);
    switch type
        case 'PP'
           a.FaceColor = "#0072B2";
        case 'MS'
           a.FaceColor = "#D55E00";
    end
    a.EdgeAlpha = 0.5;
    set( get(subplot(5,2,i),'YLabel'), 'String', 'Sensitivity' );
    ylim([0 2])
    switch type
        case 'PP'
            set(gca,'YColor','#0072B2');        
        case 'MS'
            set(gca,'YColor','#D55E00');
    end
    yyaxis right
    b = bar(linspace(1,40,40),top);
    switch type
        case 'PP'
           b.FaceColor = "#56B4E9";
        case 'MS'
           b.FaceColor = "#E69F00";
    end   
    b.EdgeAlpha = 0.5;
    set( get(subplot(5,2,i),'YLabel'), 'String', 'PPV' );
    switch type
        case 'PP'
            set(gca,'YColor','#56B4E9');        
        case 'MS'
            set(gca,'YColor','#E69F00');
    end
    set(gca, 'YDir','reverse')
    ylim([0 2])
    
    hold off
    
    
    % %adapted from this: https://www.mathworks.com/matlabcentral/answers/316005-plot-and-bar-
    % % graphs-with-independent-axis
end

function [] = Error_Exponent_bypk()
    files = {'conservative' 'MSnatural2' 'MSconservative'};
    name = {'Conservative' 'MS' 'MS-conservative'};
    nchan = 10000;
    col = [0.337,0.706,0.914,0,0.620,0.451,0,0.447,0.698];
    hold on
    for t = 1:3
        clear per_peak_mat
        load(['Results/' files{t} '10000.mat'])
        expE = abs([rstmp.ap_exp] - [vstmp.ap_exp]);
        expE = reshape(expE, [115,nchan]).';
        expE(:,116) = [rstmp(:).n_peaks];
        per_peak_mat = nan(5,1150000);
        for n = 0:max(expE(:,116))
            clear per_peak_tmp
            rows = size([expE((expE(:,116) == n), :)],1);
            per_peak_tmp(1:rows,1:115) = [expE((expE(:,116) == n), 1:115)];
            per_peak_mat(n+1,1:length([per_peak_tmp(:)])) = per_peak_tmp(:);
        end
        per_peak_mat(per_peak_mat == 0) = NaN;
        per_peak_mat = per_peak_mat';
        column_mat = nan(1150000,15);
        column_mat(:,t) = per_peak_mat(:,1);
        column_mat(:,3+t) = per_peak_mat(:,2);
        column_mat(:,6+t) = per_peak_mat(:,3);
        column_mat(:,9+t) = per_peak_mat(:,4);
        column_mat(:,12+t) = per_peak_mat(:,5);
        boxchart(column_mat, 'BoxFaceColor', col(1,(3*t-2:3*t)), 'MarkerColor','#FFFFFF', 'MarkerSize',0.01);
        ylim([0 0.5])
    end
    xticklabels(["" "0 peaks" "" "" "1 peak" "" "" "2 peaks" "" "" "3 peaks" "" "" "4 peaks" ""])
    xtickangle(0)
end

function [] = Error_Frequency()
    hold on
    files = {'conservative' 'MSnatural2' 'MSconservative'};
    name = {'Conservative' 'MS' 'MS-conservative'};
    col = [0.337,0.706,0.914,0,0.620,0.451,0,0.447,0.698];
    legend('FontSize', 6)
    for i = 1:3
        load(['Results/SPECTRA' files{i} '10000.mat'])
        freq_error = nan(length(spectra),40);
        for sim = 1:length(spectra)
            freq_error(sim, :) = spectra(sim).frequency_wise_error;
        end
        freq_error = freq_error';
        
        %taking the average
        freq_error = mean(freq_error, 2);
        plot(1:40,freq_error, 'DisplayName', name{i},'Color',[col(((i-1)*3+1):i*3)], 'LineWidth', 1.5)
        ylim([0.1 0.2])
    end
end
