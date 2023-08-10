function [] = Error_violins(data,p)
%ERROR_VIOLINS Creates a violin plot of a parameter 'p's error for each
%model
%   Data is either 'PP' or 'MS', p is 1: exp, 2: off, 3: cf, 4: amp, 5: sd.

    T = {'Exponent', 'Offset', 'Center frequency', 'Amplitude','Standard Deviation'};

    % Import the mAE of the models we want to plot, define color scheme
    switch data
        case 'PP'
            cons_mAE = open('Results/conservative10000.mat').mAE_cell;
            lib_mAE = open('Results/liberal10000.mat').mAE_cell;
            nat_mAE = open('Results/natural10000.mat').mAE_cell;
            def_mAE = open('Results/default10000.mat').mAE_cell;
            mAE_cell = [nat_mAE(:,p) def_mAE(:,p) lib_mAE(:,p) cons_mAE(:,p)];
            swatch = [0.8,0.475,0.655;0.835,0.369,0;0.941,0.894,0.259;0.337,0.706,0.914];
        case 'MS'
            cons_mAE = open('Results/conservative10000.mat').mAE_cell;
            MS_mAE = open('Results/MSnatural10000.mat').mAE_cell;
            MSdef_mAE = open('Results/MSdefault10000.mat').mAE_cell;
            MScons_mAE = open('Results/MSconservative10000.mat').mAE_cell;
            mAE_cell = [cons_mAE(:,p) MS_mAE(:,p) MSdef_mAE(:,p) MScons_mAE(:,p)];
            swatch = [0.337,0.706,0.914;0,0.620,0.451;0.902,0.624,0;0,0.447,0.698];
    end

    %Plot violins
    violin(mAE_cell,swatch,0.4);

    %Clean up figure
    yticks([-4 -3 -2 -1 0 1])
    yticklabels({'10^{-4}','10^{-3}','10^{-2}','0.1','1','10'})
    ylim([-3.2 1])
    ylabel('|error|')
    xticks(1:6)
    ax = gca(); 
    switch data
        case 'PP'
            xticklabels([{'Natural','Default','Liberal','Conservative'}]);
        case 'MS'
            xticklabels([{'Conservative','MS','MS-def','MS-cons'}]);  
    end
end


function[h,L,MX,MED,bw]=violin(Y,swatch,fa,varargin)
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
fc = swatch;

%% Calculate the kernel density
i=1;
%fa = linspace(1,1./size(Y,2),size(Y,2));
fa = repmat(fa,size(Y,2),1)'
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
