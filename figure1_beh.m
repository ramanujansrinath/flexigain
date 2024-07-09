%% load parameters for behavior in figure 1
clc; clear; close all;
load('data/figure1_beh_data.mat','figure1_beh_data')
sessIds = unique([figure1_beh_data.session]);

%% plot choice behavior for selected sessions 
figure('color','w','position',[88,506,1311,315]);
for ss=1:length(sessIds)
    sData = figure1_beh_data([figure1_beh_data.session] == sessIds(ss));
    sets = [sData.set];
    cols = lines(length(sets));
    legendStr_shape = cell(1,length(sets)*2);
    subplot(1,4,ss); hold on;
    for ii=1:length(sets)
        curv = sData([sData.set]==sets(ii)).curv;
        choice = sData([sData.set]==sets(ii)).choice;

        % choices
        plotTrendPatchLine(gca,curv,choice,cols(ii,:),false);
        
        % shape ids
        legendStr_shape{2*ii-1} = '';
        legendStr_shape{ii*2} = ['shape ' num2str(sets(ii))];
    end

    fixPlot(gca,[-0.1 1.1],[-0.1 1.1],'curvature','choice',0:0.25:1,0:0.25:1,ss,legendStr_shape)
    legend('location','southeast')
end

%%
function [xx_u,yy_m] = plotTrendPatchLine(h,xVar,yVar,col,doSmooth)
    [xx_u,~,grp] = unique(xVar);
    yy_m = groupsummary(yVar',grp,'mean'); if doSmooth; yy_m = smooth(yy_m); end
    yy_s = groupsummary(yVar',grp,'std')./sqrt(groupsummary(yVar',grp,'nnz'));

    if ~isempty(h)
        hold(h,'on')
        patch([xx_u fliplr(xx_u)]',[yy_m-yy_s/2; flipud(yy_m+yy_s/2)],col,'edgecolor','none','facealpha',0.5,'parent',h)
        plot(h,xx_u,yy_m,'linewidth',2,'color',col);
    end
end