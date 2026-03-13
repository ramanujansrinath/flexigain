%%
clc; clear; close all;
load('data/figure2_error_eg.mat','rec_v4','rec_v1')

idx = [rec_v4.exptType]==1;
rec_v1 = rec_v1(idx);
rec_v4 = rec_v4(idx);

%% example session plots
cols = lines(10);
sessIds = unique([rec_v4.sessId]);
fig1 = figure('color','w','position',[-1074,129,665,1089]);
% fig2 = figure('color','w','position',[-395,952,322,267]);
for ss=1 % 1:length(sessIds)
    figure(fig1); clf(fig1)
    exptIds = find([rec_v4.sessId]==sessIds(ss));
    setIds = [rec_v4(exptIds).setId];
    for recNum=exptIds
        col = cols(recNum-exptIds(1)+1,:);

        % for v4
        [~,~,grp] = unique([rec_v4(recNum).curv rec_v4(recNum).sel],'rows'); % mean across folds
        a4 = rec_v4(recNum).sel-rec_v4(recNum).curv; 
        a4 = groupsummary(a4,grp,'mean'); a4 = a4-mean(a4);
        b4 = rec_v4(recNum).predcc_gen; b4 = groupsummary(b4,grp,'mean');
        c4 = rec_v4(recNum).predcc_gen-rec_v4(recNum).curv;  c4 = groupsummary(c4,grp,'mean');
        d4 = rec_v4(recNum).curv; d4 = groupsummary(d4,grp,'mean');
        e4 = rec_v4(recNum).sel; e4 = groupsummary(e4,grp,'mean');
        f4 = rec_v4(recNum).predcc_spe-rec_v4(recNum).curv;  f4 = groupsummary(f4,grp,'mean');
        g4 = rec_v4(recNum).predcc_spe; g4 = groupsummary(g4,grp,'mean');

        % for v1
        [~,~,grp] = unique([rec_v1(recNum).curv rec_v1(recNum).sel],'rows'); % mean across folds
        a1 = rec_v1(recNum).sel-rec_v1(recNum).curv; 
        a1 = groupsummary(a1,grp,'mean'); a1 = a1-mean(a1);
        b1 = rec_v1(recNum).predcc_gen; b1 = groupsummary(b1,grp,'mean');
        c1 = rec_v1(recNum).predcc_gen-rec_v1(recNum).curv;  c1 = groupsummary(c1,grp,'mean');
        d1 = rec_v1(recNum).curv; d1 = groupsummary(d1,grp,'mean');
        e1 = rec_v1(recNum).sel; e1 = groupsummary(e1,grp,'mean');
        f1 = rec_v1(recNum).predcc_spe-rec_v1(recNum).curv;  f1 = groupsummary(f1,grp,'mean');
        g1 = rec_v1(recNum).predcc_spe; g1 = groupsummary(g1,grp,'mean');

        % for v4
        subplot(434); hold on;
        plotTrendPatchLine(gca,d4',b4',col,[0 1],9); co4(recNum-exptIds(1)+1) = corr(d4,b4);
        fixPlot(gca,[-0.1 1.1],[-0.1 1.1],'curvature','dec',-1:0.25:1,-1:0.25:1,'v4 gen dec')

        subplot(435); hold on;
        plotTrendPatchLine(gca,e4',b4',col,[0 1],9); co5(recNum-exptIds(1)+1) = corr(b4,e4);
        fixPlot(gca,[-0.1 1.1],[-0.1 1.1],'choices','dec',-1:0.25:1,-1:0.25:1,'v4 gen dec v choice')

        subplot(436); hold on;
        plotTrendPatchLine(gca,e4',g4',col,[0 1],9); co6(recNum-exptIds(1)+1) = corr(g4,e4);
        fixPlot(gca,[-0.1 1.1],[-0.1 1.1],'choices','dec',-1:0.25:1,-1:0.25:1,'v4 spe dec v choice')

        subplot(4,3,10); hold on;
        plotTrendPatchLine(gca,a4',c4',col,[-0.5 0.5],9); co10(recNum-exptIds(1)+1) = corr(a4,c4);
        fixPlot(gca,[-0.6 0.6],[-0.6 0.6],'beh error','dec error',-1:0.25:1,-1:0.25:1,'v4 gen dec error')

        subplot(4,3,11); hold on;
        plotTrendPatchLine(gca,a4',f4',col,[-0.5 0.5],9); co11(recNum-exptIds(1)+1) = corr(a4,f4);
        fixPlot(gca,[-0.6 0.6],[-0.6 0.6],'beh error','dec error',-1:0.25:1,-1:0.25:1,'v4 spe dec error')
    
        
        % for v1
        subplot(431); hold on;
        plotTrendPatchLine(gca,d1',b1',col,[0 1],9); co1(recNum-exptIds(1)+1) = corr(d1,b1);
        fixPlot(gca,[-0.1 1.1],[-0.1 1.1],'curvature','dec',-1:0.25:1,-1:0.25:1,'v1 gen dec')

        subplot(432); hold on;
        plotTrendPatchLine(gca,e1',b1',col,[0 1],9); co2(recNum-exptIds(1)+1) = corr(b1,e1);
        fixPlot(gca,[-0.1 1.1],[-0.1 1.1],'choices','dec',-1:0.25:1,-1:0.25:1,'v1 gen dec v choice')

        subplot(433); hold on;
        plotTrendPatchLine(gca,e1',g1',col,[0 1],9); co3(recNum-exptIds(1)+1) = corr(e1,g1);
        fixPlot(gca,[-0.1 1.1],[-0.1 1.1],'choices','dec',-1:0.25:1,-1:0.25:1,'v1 spe dec v choice')

        subplot(437); hold on;
        plotTrendPatchLine(gca,a1',c1',col,[-0.5 0.5],9); co7(recNum-exptIds(1)+1) = corr(a1,c1);
        fixPlot(gca,[-0.6 0.6],[-0.6 0.6],'beh error','dec error',-1:0.25:1,-1:0.25:1,'v1 gen dec error')

        subplot(438); hold on;
        plotTrendPatchLine(gca,a1',f1',col,[-0.5 0.5],9); co8(recNum-exptIds(1)+1) = corr(a1,f1);
        fixPlot(gca,[-0.6 0.6],[-0.6 0.6],'beh error','dec error',-1:0.25:1,-1:0.25:1,'v1 spe dec error')

    end
    sgtitle([num2str(sessIds(ss)) ': ' num2str(setIds)])

    for ii=[1 2 3 4 5 6 7 8 10 11]
        subplot(4,3,ii);
        legend(num2str(round(eval(['co' num2str(ii)]),2)))
    end

    % figure(fig2); clf(fig2)
    % for recNum=exptIds
    %     col = cols(recNum-exptIds(1)+1,:);
    % 
    %     [~,~,grp] = unique([rec_v4(recNum).curv rec_v4(recNum).sel],'rows'); % mean across folds
    %     a4 = rec_v4(recNum).sel-rec_v4(recNum).curv; 
    %     a4 = groupsummary(a4,grp,'mean'); a4 = a4-mean(a4);
    %     d4 = rec_v4(recNum).curv; d4 = groupsummary(d4,grp,'mean');
    % 
    %     % behavior
    %     hold on;
    %     plotTrendPatchLine(gca,d4',a4',col,[0 1],9);
    %     fixPlot(gca,[-0.1 1.1],[-0.6 0.6],'curvature','beh error',0:0.1:1,-1:0.25:1,'beh error')
    % end
end



%%
function [xx_u,yy_m,grp] = plotTrendPatchLine(h,xVar,yVar,col,xLims,nBin)
    % [xx_u,~,grp] = unique(xVar);
    [~,binE,grp] = histcounts(xVar,linspace(xLims(1),xLims(2),nBin+1));
    grp(grp<1) = 1; grp(grp>nBin) = nBin;
    xx_u = (binE+circshift(binE,-1))/2; xx_u = xx_u(1:end-1)';
    yy_m = groupsummary(yVar',grp','mean'); % yy_m = smooth(yy_m);
    yy_s = groupsummary(yVar',grp','std')./sqrt(groupsummary(yVar',grp','nnz'));

    % fix for when there are empty bins
    if length(yy_m)~=length(xx_u)
        % yy_m(unique(grp)) = yy_m; yy_m(~ismember(1:nBin,unique(grp))) = nan;
        % yy_s(unique(grp)) = yy_s; yy_s(~ismember(1:nBin,unique(grp))) = nan;
        emptyBins = 1:nBin;
        emptyBins(ismember(emptyBins,unique(grp))) = [];
        xx_u(emptyBins) = [];
    end

    if ~isempty(h)
        hold(h,'on')
        patch([xx_u;flipud(xx_u)],[yy_m-yy_s/2; flipud(yy_m+yy_s/2)],col,'edgecolor','none','facealpha',0.5,'parent',h)
    
        beta = regress(yy_m(:),[ones(length(xx_u),1) xx_u(:)]);
        yFit = beta(1) + beta(2)*xLims;
        line(xLims,yFit,'linestyle','--','color',col,'linewidth',2)
    
        plot(h,xx_u,yy_m,'-','linewidth',2,'color',col); % 'MarkerFaceColor','w','markersize',10
    end
end