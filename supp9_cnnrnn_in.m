clc; close all; clear;
load('data/supp9_cnnrnn_in_model1.mat')
load('data/supp9_cnnrnn_in_model2.mat')

%%
figure('color','w','position',[476,211,560,569]); 

% image set 1 training data curvature decoding
plotCurvatureDec(subplot(221),in_model1.train.cc,in_model1.train.act_red,'set 1 train')

% image set 1 test data curvature decoding
plotCurvatureDec(subplot(222),in_model1.test.cc,in_model1.test.act_red,'set 1 test')

% image set 2 training data curvature decoding
plotCurvatureDec(subplot(223),in_model2.train.cc,in_model2.train.act_red,'set 2 train')

% image set 2 test data curvature decoding
plotCurvatureDec(subplot(224),in_model2.test.cc,in_model2.test.act_red,'set 2 test')

%%
function plotCurvatureDec(h,cc,act_red,titleStr)
    hold(h,'on'); 
    dec = decodeFn(cc',act_red); dec(dec>2) = nan; dec(dec<-2) = nan;
    plot(h,cc,dec,'.','markersize',10,'color',[0.5 0.5 0.5]);
    plotTrendPatchLine(h,cc',dec,'k',[0 1],21);
    fixPlot(h,[-0.1 1.1],[-0.1 1.1],'curvature','decoded curvature',0:0.25:1,0:0.25:1,titleStr)
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