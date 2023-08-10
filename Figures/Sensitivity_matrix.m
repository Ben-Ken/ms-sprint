function [] = Sensitivity_matrix(filename,col)
%SENSITIVITY_MATRIX Creates a matrix with percentage chances of detecting y
%peaks when x peaks are simulated

    %Load relevant file
    load(['Results/' filename '10000.mat'])

    %Create a matrix with every timepoint
    ts(:) = linspace(1.5, 58.5, 115);

    %Initializing matrices
    npeak_bic_mat = zeros(5,7);
    det_v_sim_mat = zeros(length(rstmp), length(ts), 2);
    
    %Populate det_v_sim_mat with the detected v. simulated peaks for every
    %channel-time point
    for sim = 1:length(rstmp)
        for time = 1:length(ts)
            if isempty(rstmp(sim).pk_det) %if no peaks were simulated
                det_v_sim_mat(sim,time,1) = sum([rstmp(sim).resid_peaks(:).time] == ts(time)); %if no peaks were simulated
                det_v_sim_mat(sim,time,2) = 0; 
            else %if peaks were simulated
                det_v_sim_mat(sim,time,1) = sum([rstmp(sim).resid_peaks(:).time] == ts(time)) + sum(vstmp(sim).pk_exist(:,time)); %if peaks were simulated
                det_v_sim_mat(sim,time,2) = sum(vstmp(sim).pk_exist(:,time)); 
            end
        end
    end
    
    %Populate npeak_bic_mat with the sum of channel-time points that
    %sort into the cell with the given amount of simulated-detected peaks
    for simulated = 1:size(npeak_bic_mat,1)
        for detected = 1:size(npeak_bic_mat,2)
            npeak_bic_mat(simulated,detected) = sum(logical(logical(det_v_sim_mat(:,:,1) == (detected - 1) & logical(det_v_sim_mat(:,:,2) == (simulated - 1)))),'all');
        end
    end
    
    %turn into ratio row-wise (multiplied by 1000 just for proper
    %color-mapping)
    data = nan(5,5);
    for row = 1:size(npeak_bic_mat,1);
        data(row,:) = npeak_bic_mat(row,1:5)./sum(npeak_bic_mat(row,:))*1000;
    end

    %colors
    c3 = {[0,0.447,0.698], [0.835,0.369,0]};
    c3 = c3{col};
    hold on

    %define dimensions
    xs = [0 1 2 3 4];
    ys = [0 1 2 3 4];

    %fill in cell colors and text
    imagesc(xs,ys, data')
    for x = 1:5
      for y = 1:5
        text(xs(x),ys(y),num2str(round(data(x,y)/1000,3)*100),'HorizontalAlignment','center','VerticalAlignment','middle') 
      end
    end

    %clean up figure
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
end
