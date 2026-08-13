function [] = plot_profile_likelihood_nL_nV()

    clc;
    close all;

    n_L_span = 12:1:24;
    MLvalues = load('MLvalues.mat'); 
    MLvalues = MLvalues.MLvalues;

    output = [];
    for ii = 1:13
        surfvals = load(strcat('./PL_nL_nV/nL_',num2str(n_L_span(ii)),'_nV.mat'));
        surfvals = surfvals.funval;
        output = [output; surfvals];
    end

    value = 9999.9;
    for ii=1:size(output,1)
        if output(ii,3) < value
            value = output(ii,3);
            pars = output(ii,1:2);
        end
    end
    disp(strcat('value = ', num2str(value),...
                ' n_L = ', num2str(pars(1)),...
                ' n_V = ', num2str(pars(2))))


    figure();

    xlin = linspace(min(output(:,1)), max(output(:,1)), 13);
    ylin = linspace(min(output(:,2)), max(output(:,2)), 13);

    [X,Y] = meshgrid(xlin, ylin);
    Z = griddata(output(:,1),output(:,2),output(:,3),X,Y,'v4');
    S = mesh(X,Y,Z);
    axis tight; 
    box on;
    hold on;
    plot3(MLvalues(1,end-1),MLvalues(1,end),-MLvalues(1,1),'Marker','o',...
                                              'MarkerSize',8,...
                                              'LineStyle','none',...
                                              'MarkerFaceColor','red',...
                                              'MarkerEdgeColor','red');
    plot3([MLvalues(1,end-1) MLvalues(1,end-1)],[MLvalues(1,end) MLvalues(1,end)],[min(output(:,3)) max(output(:,3))],'LineStyle','-','Color','red','LineWidth',2);
    hold on;
    xlabel('n_L'); 
    ylabel('n_V');
    zlabel('Negative profile likelihood');
    S.FaceColor = 'flat';
    S(1).EdgeColor = 'black';   
    % Rotate the camera: 45 degrees azimuth, 30 degrees elevation
    view(-45, 30);

    exportgraphics(gcf,strcat('../figures/profile_likelihood_nL_nV.png'),'Resolution',600);
    exportgraphics(gcf,strcat('../../../LaTeX/figures/profile_likelihood_nL_nV.eps'),'ContentType','vector');

end