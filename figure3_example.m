%% load a random session
clc; close all; clear
load('data/figure3_example.mat');

%% V1
figure('color','w','position',[100,272,879,528],'name','V1');

resp_stim_mid = figure3_example_v1_mid.resp_stim;
resp_arc_mid = figure3_example_v1_mid.resp_arc;
p_set_mid = figure3_example_v1_mid.params;;
plot_curvDec([subplot(231) subplot(232)],p_set_mid,resp_stim_mid,resp_arc_mid)
plot_saccMidDec(subplot(233),p_set_mid,resp_arc_mid)

resp_stim_len = figure3_example_v1_len.resp_stim;
resp_arc_len = figure3_example_v1_len.resp_arc;
p_set_len = figure3_example_v1_len.params;
plot_curvDec([subplot(234) subplot(235)],p_set_len,resp_stim_len,resp_arc_len)
plot_saccLenDec(subplot(236),p_set_len,resp_arc_len)

%%
figure('color','w','position',[997,272,879,528],'name','V4');

resp_stim_mid = figure3_example_v4_mid.resp_stim;
resp_arc_mid = figure3_example_v4_mid.resp_arc;
p_set_mid = figure3_example_v4_mid.params;;
plot_curvDec([subplot(231) subplot(232)],p_set_mid,resp_stim_mid,resp_arc_mid)
plot_saccMidDec(subplot(233),p_set_mid,resp_arc_mid)

resp_stim_len = figure3_example_v4_len.resp_stim;
resp_arc_len = figure3_example_v4_len.resp_arc;
p_set_len = figure3_example_v4_len.params;
plot_curvDec([subplot(234) subplot(235)],p_set_len,resp_stim_len,resp_arc_len)
plot_saccLenDec(subplot(236),p_set_len,resp_arc_len)

%%
function plot_curvDec(ha,p,resp_stim,resp_arc)
    % stim/arc time curvature decoding
    stim_col = [103 221 103]/255;
    pred_curv = decodeFn([p.curv]',resp_stim);
    axes(ha(1));
    jitter_c = [p.curv] + 0.1*(-0.5+rand(size([p.curv])));
    scatter(jitter_c,pred_curv,20,stim_col,'filled','MarkerFaceAlpha',0.3); hold on;
    [xx_u,yy_m] = plotTrendPatchLine(gca,[p.curv]',pred_curv,stim_col,[0 1],9);
    plot(gca,xx_u,yy_m,'o','linewidth',2,'color',stim_col,'MarkerFaceColor','w','markersize',10)
    
    stim_col = [221 103 103]/255;
    pred_curv = decodeFn([p.curv]',resp_arc);
    scatter(jitter_c,pred_curv,20,stim_col,'filled','MarkerFaceAlpha',0.3); hold on;
    [xx_u,yy_m] = plotTrendPatchLine(gca,[p.curv]',pred_curv,stim_col,[0 1],9);
    plot(gca,xx_u,yy_m,'o','linewidth',2,'color',stim_col,'MarkerFaceColor','w','markersize',10)
    
    fixPlot(gca,[-0.1 1.1],[-0.1 1.1],'curvature',{'pred curvature' '(stim time)'},0:0.25:1,0:0.25:1,'curvature decoding',{'' '' '' 'stim time' '' '' '' '' 'arc time' ''})
    legend('location','northwest')
    
    % arc time saccade decoding
    axes(ha(2));
    stim_col = [221 103 103]/255;
    pred_sacc = decodeFn([p.sacc_th]',resp_arc);
    scatter([p.sacc_th],pred_sacc,20,stim_col,'filled','MarkerFaceAlpha',0.3); hold on;
    [xx_u,yy_m] = plotTrendPatchLine(gca,[p.sacc_th]',pred_sacc,stim_col,[-65 103],10);
    plot(gca,xx_u,yy_m,'o','linewidth',2,'color',stim_col,'MarkerFaceColor','w','markersize',10)
    fixPlot(gca,[-70 105],[-70 105],'saccade th',{'pred saccade th' '(arc time)'},-180:45:180,-180:45:180,'sacc dec')
end

function plot_saccMidDec(h,p_set_mid,resp_arc_mid)
    % arc time saccade decoding (per arc mid)
    pred_sacc = decodeFn([p_set_mid.sacc_th]',resp_arc_mid);
    sacc_lim_x = [-70 90];
    sacc_lim_y = [-50 50];
    
    % arc mids
    uu = unique([p_set_mid.arc_st]);
    cols = winter(length(uu));
    axes(h); hold on;
    for ii=1:length(uu)
        idx_uu = [p_set_mid.arc_st]==uu(ii) & [p_set_mid.arc_len]== 140;
        sacc_uu = [p_set_mid(idx_uu).sacc_th];
        pred_sacc_uu = pred_sacc(idx_uu);
    
        localLims = [min([p_set_mid(idx_uu).sacc_th_correct]) max([p_set_mid(idx_uu).sacc_th_correct])];
        plotTrendPatchLine(gca,sacc_uu',pred_sacc_uu,cols(ii,:),localLims,10);
    end
    fixPlot(gca,sacc_lim_x,sacc_lim_y,'saccade th',{'pred saccade th' '(arc time)'},-180:45:180,-180:45:180,'sacc dec')
end

function plot_saccLenDec(h,p_set_len,resp_arc_len)
    % arc time saccade decoding (per arc len)
    pred_sacc = decodeFn([p_set_len.sacc_th]',resp_arc_len);
    sacc_lim_x = [-70 90];
    sacc_lim_y = [-50 50];
    
    % arc lens
    uu = unique([p_set_len.arc_len]);
    cols = [copper(4);winter(3)]; cols = cols([3 6],:);
    axes(h); hold on;
    for ii=1:length(uu)
        idx_uu = [p_set_len.arc_len]==uu(ii) & [p_set_len.arc_mid]== 20;
        sacc_uu = [p_set_len(idx_uu).sacc_th];
        pred_sacc_uu = pred_sacc(idx_uu);
    
        localLims = [min([p_set_len(idx_uu).sacc_th_correct]) max([p_set_len(idx_uu).sacc_th_correct])];
        plotTrendPatchLine(gca,sacc_uu',pred_sacc_uu,cols(ii,:),localLims,10);
    end
    fixPlot(gca,sacc_lim_x,sacc_lim_y,'saccade th',{'pred saccade th' '(arc time)'},-180:45:180,-180:45:180,'sacc dec')
end

function [xx_u,yy_m] = plotTrendPatchLine(h,xVar,yVar,col,xLims,nBin)
    % [xx_u,~,grp] = unique(xVar);
    [~,binE,grp] = histcounts(xVar,linspace(xLims(1),xLims(2),nBin+1));
    grp(grp<1) = 1; grp(grp>nBin) = nBin;
    xx_u = (binE+circshift(binE,-1))/2; xx_u = xx_u(1:end-1)';
    yy_m = groupsummary(yVar,grp,'mean'); yy_m = smooth(yy_m);
    yy_s = groupsummary(yVar,grp,'std')./sqrt(groupsummary(yVar,grp,'nnz'));

    % fix for when there are empty bins
    if length(yy_m)~=length(xx_u)
        yy_m(unique(grp)) = yy_m; yy_m(~ismember(1:nBin,unique(grp))) = nan;
        yy_s(unique(grp)) = yy_s; yy_s(~ismember(1:nBin,unique(grp))) = nan;
    end

    hold(h,'on')
    patch([xx_u;flipud(xx_u)],[yy_m-yy_s/2; flipud(yy_m+yy_s/2)],col,'edgecolor','none','facealpha',0.5,'parent',h)

    beta = regress(yy_m(:),[ones(length(xx_u),1) xx_u(:)]);
    yFit = beta(1) + beta(2)*xLims;
    line(xLims,yFit,'linestyle','--','color',col,'linewidth',2)

    plot(h,xx_u,yy_m,'-','linewidth',2,'color',col); % 'MarkerFaceColor','w','markersize',10
end