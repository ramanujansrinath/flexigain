%%
clc; clear; close all;
load('data/figure2_cross.mat')

% this code also generates supp5

%% one session
setList = 1:120;
setGroups_ori = arrayfun(@(ii) ii:20:80,1:20,'UniformOutput',false);
setGroups_col = arrayfun(@(ii) [ii ii+80 ii+100],1:20,'UniformOutput',false);

red_v4 = pca(resp_v4','NumComponents',10);
red_v1 = pca(resp_v1','NumComponents',10);

sets = [17 16];
% sets = unique([params.set]); sets = randsample(sets,2,false);
cols = lines(length(sets));
legendStr_shape = cell(1,length(sets)*2);
figure('color','w','position',[1,345,1512,521]);
for ii=1:length(sets)
    % get this set
    setP = params([params.set] == sets(ii));
    cc = [setP.curv];
    red_v1_set = red_v1([params.set] == sets(ii),:);
    red_v4_set = red_v4([params.set] == sets(ii),:);
    
    % get opposite set
    xsetP = params([params.set] == sets(3-ii));
    xcc = [xsetP.curv];
    red_v1_xset = red_v1([params.set] == sets(3-ii),:); 
    red_v4_xset = red_v4([params.set] == sets(3-ii),:); 
    
    % choices across sets
    if ii==1
        ss = [setP.sel];
        [yy1_m,yy1_s] = plotTrendPatchLine([],cc,ss,[],[0 1],19);
        xss = [xsetP.sel];
        [yy2_m,yy2_s] = plotTrendPatchLine([],xcc,xss,[],[0 1],19);

        arrayfun(@(jj) cla(subplot(2,6,jj)),[1 2 3 7 8]);
        beta = regress(yy2_m,[ones(length(yy1_m),1) yy1_m]);
        line(subplot(263),[0 1],[beta(1) sum(beta)],'linestyle','--','linewidth',2,'color','k'); hold on;
        line(subplot(263),[0 1],[0 1],'linestyle','--','linewidth',2,'color','k'); hold on;

        arrayfun(@(kk) line([yy1_m(kk)-yy1_s(kk)/2 yy1_m(kk)+yy1_s(kk)/2],[yy2_m(kk) yy2_m(kk)],'linewidth',2,'color','k'),1:length(yy1_m));
        arrayfun(@(kk) line([yy1_m(kk) yy1_m(kk)],[yy2_m(kk)-yy2_s(kk)/2 yy2_m(kk)+yy2_s(kk)/2],'linewidth',2,'color','k'),1:length(yy1_m));
        plot(subplot(263),yy1_m,yy2_m,'ko','markerfacecolor','w','linewidth',2);
    end

    % mean subtract?
    % red_v1_set = red_v1_set-repmat(mean(red_v1_set),size(red_v1_set,1),1);
    % red_v4_set = red_v4_set-repmat(mean(red_v4_set),size(red_v4_set,1),1);
    % red_v1_xset = red_v1_xset-repmat(mean(red_v1_xset),size(red_v1_xset,1),1);
    % red_v4_xset = red_v4_xset-repmat(mean(red_v4_xset),size(red_v4_xset,1),1);

    % dec v4
    [dec_cc_v4,dec_xcc_v4] = cv_xset_decode(cc,red_v4_set,xcc,red_v4_xset);
    plotTrendPatchLine(subplot(261),cc,dec_cc_v4,cols(ii,:),[0 1],19); mean((cc-dec_cc_v4).^2)
    plotTrendPatchLine(subplot(262),xcc,dec_xcc_v4,cols(ii,:)/2,[0 1],19); mean((xcc-dec_xcc_v4).^2)
    
    % dec v1
    [dec_cc_v1,dec_xcc_v1] = cv_xset_decode(cc,red_v1_set,xcc,red_v1_xset);
    plotTrendPatchLine(subplot(267),cc,dec_cc_v1,cols(ii,:),[0 1],19); mean((cc-dec_cc_v1).^2)
    plotTrendPatchLine(subplot(268),xcc,dec_xcc_v1,cols(ii,:)/2,[0 1],19); mean((xcc-dec_xcc_v1).^2)

    % shape ids
    legendStr_shape{2*ii-1} = '';
    legendStr_shape{ii*2} = ['shape ' num2str(sets(ii))];
end

fixPlot(subplot(263),[-0.1 1.1],[-0.1 1.1],'choice set 1','choice set 2',0:0.25:1,0:0.25:1,'choice correlation')
fixPlot(subplot(261),[-0.1 1.1],[-0.1 1.1],'curvature','decoded curvature',0:0.25:1,0:0.25:1,'v4 decoding self',legendStr_shape)
legend('location','southeast')
fixPlot(subplot(262),[-0.1 1.1],[-0.1 1.1],'curvature','decoded curvature',0:0.25:1,0:0.25:1,'v4 decoding cross')
fixPlot(subplot(267),[-0.1 1.1],[-0.1 1.1],'curvature','decoded curvature',0:0.25:1,0:0.25:1,'v1 decoding self')
fixPlot(subplot(268),[-0.1 1.1],[-0.1 1.1],'curvature','decoded curvature',0:0.25:1,0:0.25:1,'v1 decoding cross')


%% across sessions
% choice similarity
subplot(269); cla;
line([mean(choice_corr) mean(choice_corr)],[0 400],'linestyle','--','color','k','linewidth',2);  hold on;
histogram(choice_corr(exptType==1),linspace(0,1,25),'displayStyle','stairs','LineWidth',2,'edgecolor','k');
histogram(choice_corr(exptType==2),linspace(0,1,25),'displayStyle','stairs','LineWidth',2,'edgecolor','b');
histogram(choice_corr(exptType==3),linspace(0,1,25),'displayStyle','stairs','LineWidth',2,'edgecolor','r');
fixPlot(gca,[-0.1 1],[0 200],'corr between choices','shape pairs',0:0.25:1,0:50:400,{'beh similarity' 'across shapes'})

% v4
subplot(264); cla
line([0 1],[0 1],'linestyle','--','color','k','linewidth',2); hold on;
% plot(mse_v4_same,mse_v4_opp,'ko','MarkerSize',6,'MarkerFaceColor','w','LineWidth',2); 
scatter(mse_v4_same(exptType==1),mse_v4_opp(exptType==1),20,'k','filled');
scatter(mse_v4_same(exptType==2),mse_v4_opp(exptType==2),20,'b','filled');
scatter(mse_v4_same(exptType==3),mse_v4_opp(exptType==3),20,'r','filled');
fixPlot(gca,[-0.1 1],[-0.1 1],'same set','other set',0:0.25:1,0:0.25:1,'v4')
% fixPlot(gca,[0 0.1],[0 0.1],'same set','other set',0:0.25:1,0:0.25:1,'v4')

subplot(265); cla; hold on;
% line([mean(mse_v4_same) mean(mse_v4_same)],[0 2000],'linestyle','--','color','k','linewidth',2);  hold on;
histogram(mse_v4_same(exptType==1),linspace(0,1.01,25),'displayStyle','stairs','LineWidth',2,'edgecolor','k');
histogram(mse_v4_same(exptType==2),linspace(0,1.01,25),'displayStyle','stairs','LineWidth',2,'edgecolor','b');
histogram(mse_v4_same(exptType==3),linspace(0,1.01,25),'displayStyle','stairs','LineWidth',2,'edgecolor','r');
fixPlot(gca,[-0.1 1],[0 420],'same set','count',0:0.25:1,0:400:1200)

subplot(266); cla; hold on;
% line([mean(mse_v4_opp) mean(mse_v4_opp)],[0 400],'linestyle','--','color','k','linewidth',2); hold on;
histogram(mse_v4_opp(exptType==1),linspace(0,1.01,25),'displayStyle','stairs','LineWidth',2,'edgecolor','k');
histogram(mse_v4_opp(exptType==2),linspace(0,1.01,25),'displayStyle','stairs','LineWidth',2,'edgecolor','b');
histogram(mse_v4_opp(exptType==3),linspace(0,1.01,25),'displayStyle','stairs','LineWidth',2,'edgecolor','r');
fixPlot(gca,[-0.1 1],[0 210],'other set','count',0:0.25:1,0:50:300)

% v1
subplot(2,6,10); cla 
line([0 1],[0 1],'linestyle','--','color','k','linewidth',2); hold on;
% plot(mse_v1_same,mse_v1_opp,'ko','MarkerSize',6,'MarkerFaceColor','w','LineWidth',2); 
scatter(mse_v1_same(exptType==1),mse_v1_opp(exptType==1),20,'k','filled');
scatter(mse_v1_same(exptType==2),mse_v1_opp(exptType==2),20,'b','filled');
scatter(mse_v1_same(exptType==3),mse_v1_opp(exptType==3),20,'r','filled');
fixPlot(gca,[-0.1 1],[-0.1 1],'same set','other set',0:0.25:1,0:0.25:1,'v1')
% fixPlot(gca,[0 0.1],[0 0.1],'same set','other set',0:0.25:1,0:0.25:1,'v1')

subplot(2,6,11); cla; hold on;
% line([mean(mse_v1_same) mean(mse_v1_same)],[0 2000],'linestyle','--','color','k','linewidth',2);  hold on;
histogram(mse_v1_same(exptType==1),linspace(0,1.01,25),'displayStyle','stairs','LineWidth',2,'edgecolor','k');
histogram(mse_v1_same(exptType==2),linspace(0,1.01,25),'displayStyle','stairs','LineWidth',2,'edgecolor','b');
histogram(mse_v1_same(exptType==3),linspace(0,1.01,25),'displayStyle','stairs','LineWidth',2,'edgecolor','r');
fixPlot(gca,[-0.1 1],[0 420],'same set','count',0:0.25:1,0:400:1200)

subplot(2,6,12); cla; hold on;
% line([mean(mse_v1_opp) mean(mse_v1_opp)],[0 400],'linestyle','--','color','k','linewidth',2); hold on;
histogram(mse_v1_opp(exptType==1),linspace(0,1.01,25),'displayStyle','stairs','LineWidth',2,'edgecolor','k');
histogram(mse_v1_opp(exptType==2),linspace(0,1.01,25),'displayStyle','stairs','LineWidth',2,'edgecolor','b');
histogram(mse_v1_opp(exptType==3),linspace(0,1.01,25),'displayStyle','stairs','LineWidth',2,'edgecolor','r');
fixPlot(gca,[-0.1 1],[0 210],'other set','count',0:0.25:1,0:50:300,'',{'shape' 'ori' 'color'})
legend('Location','northeast')

%%
function [yy_m,yy_s] = plotTrendPatchLine(h,xVar,yVar,col,xLims,nBin)
    % [xx_u,~,grp] = unique(xVar);
    [~,binE,grp] = histcounts(xVar,linspace(xLims(1),xLims(2),nBin+1));
    grp(grp<1) = 1; grp(grp>nBin) = nBin;
    xx_u = (binE+circshift(binE,-1))/2; xx_u = xx_u(1:end-1)';
    yy_m = groupsummary(yVar',grp','mean');
    yy_s = groupsummary(yVar',grp','std')./sqrt(groupsummary(yVar',grp','nnz'));

    % fix for when there are empty bins
    if length(yy_m)~=length(xx_u)
        yy_m(unique(grp)) = yy_m; yy_m(~ismember(1:nBin,unique(grp))) = nan;
        yy_s(unique(grp)) = yy_s; yy_s(~ismember(1:nBin,unique(grp))) = nan;

        yy_m = fillmissing(yy_m,'movmean',3);  
    end

    if ~isempty(h)
        hold(h,'on')
        patch([xx_u;flipud(xx_u)],[yy_m-yy_s/2; flipud(yy_m+yy_s/2)],col,'edgecolor','none','facealpha',0.5,'parent',h)
        plot(h,xx_u,yy_m,'linewidth',2,'color',col);
    
        beta = regress(yy_m(:),[ones(length(xx_u),1) xx_u(:)]);
        yFit = beta(1) + beta(2)*xLims;
        line(xLims,yFit,'linestyle','--','color',col,'linewidth',2)
    end
end

function [cc_pred,xcc_pred] = cv_xset_decode(cc_set1,rr_set1,cc_set2,rr_set2)
    nFold = 1000;
    cc_pred = nan(length(cc_set1),nFold);
    xcc_pred = nan(length(cc_set2),nFold);
    for ii=1:nFold
        tt_train = true(1,length(cc_set1));
        tt_train(randperm(length(tt_train),round(length(tt_train)/2))) = false;
        tt_test = ~tt_train;
        
        xtt_test = false(1,length(cc_set2));
        xtt_test(randperm(length(xtt_test),min(length(xtt_test),sum(tt_test)))) = true;
        
        cc_pred(tt_test,ii) = rr_set1(tt_test,:) * regress(cc_set1(tt_train)',rr_set1(tt_train,:));
        
        xcc_pred(xtt_test,ii) = rr_set2(xtt_test,:) * regress(cc_set1(tt_train)',rr_set1(tt_train,:));
    end
    cc_pred = nanmean(cc_pred,2)'; %#ok<NANMEAN> 
    xcc_pred = nanmean(xcc_pred,2)'; %#ok<NANMEAN> 
end