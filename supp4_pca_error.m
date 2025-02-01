%%
clc; clear; close all;
plotSets = {4; [12 13 14]; [5 25]; [13 93 113]};
load('data/supp4_pca_data.mat')

%% shape-specific decoders
sessId = 4;
trialMarkerAl = 0.2; % draw trial dots in low alpha
sets = plotSets{sessId};
resp_v4 = supp4_pca_data(sessId).resp_v4;
resp_v1 = supp4_pca_data(sessId).resp_v1;

red_v4 = pca(resp_v4','NumComponents',10);
red_v1 = pca(resp_v1','NumComponents',10);

cols = lines(length(sets));
legendStr_shape = cell(1,length(sets)*2);
legendStr_v4_dec = cell(1,length(sets)*2);
legendStr_v1_dec = cell(1,length(sets)*2);

figure('color','w','position',[-1306,146,508,844]);
for ii=1:length(sets)
    idx = [supp4_pca_data(sessId).set] == sets(ii);
    cc = supp4_pca_data(sessId).curv(idx);
    ss = supp4_pca_data(sessId).sel(idx);

    red_v1_set = red_v1(idx,:);
    red_v4_set = red_v4(idx,:);

    % pca or tuning v4
    plotPCfn(subplot(321),red_v4_set,cc,cols(ii,:),'V4',trialMarkerAl);
    
    % pca v1
    plotPCfn(subplot(322),red_v1_set(:,1:3),cc,cols(ii,:),'V1',trialMarkerAl);

    % dec v4
    dec_cc_v4 = decodeFn(cc',red_v4_set);
    plotTrendPatchLine(subplot(323),cc,dec_cc_v4',cols(ii,:),[0 1],10);

    % dec v1
    dec_cc_v1 = decodeFn(cc',red_v1_set);
    plotTrendPatchLine(subplot(324),cc,dec_cc_v1',cols(ii,:),[0 1],10);
    
    % err v4
    dec_cc_v4 = decodeFn(cc',red_v4_set);
    cc_err = (dec_cc_v4' - cc); cc_err = cc_err-mean(cc_err);
    ss_err = (dec_cc_v4' - ss); ss_err = ss_err-mean(ss_err);
    plotTrendPatchLine(subplot(325),cc_err,ss-cc,cols(ii,:),[-0.5 0.5],10);

    % err v1
    dec_cc_v1 = decodeFn(cc',red_v1_set);
    cc_err = (dec_cc_v1' - cc); cc_err = cc_err-mean(cc_err);
    ss_err = (dec_cc_v1' - ss); ss_err = ss_err-mean(ss_err);
    plotTrendPatchLine(subplot(326),cc_err,ss-cc,cols(ii,:),[-0.5 0.5],10);

    % shape ids
    legendStr_shape{3*ii-2} = ''; legendStr_shape{3*ii-1} = '';
    legendStr_shape{ii*3} = ['shape ' num2str(sets(ii))];

    % dec acc v4
    legendStr_v4_dec{3*ii-2} = ''; legendStr_v4_dec{3*ii-1} = '';
    legendStr_v4_dec{ii*3} = ['acc = ' num2str(round(corr(dec_cc_v4,cc'),2))];

    % dec acc v1
    legendStr_v1_dec{3*ii-2} = ''; legendStr_v1_dec{3*ii-1} = '';
    legendStr_v1_dec{ii*3} = ['acc = ' num2str(round(corr(dec_cc_v1,cc'),2))];
end

fixPlot(subplot(323),[-0.1 1.1],[-0.1 1.1],'curvature','decoded curvature',0:0.25:1,0:0.25:1,'v4 decoding',legendStr_v4_dec)
legend('location','southeast')
fixPlot(subplot(324),[-0.1 1.1],[-0.1 1.1],'curvature','decoded curvature',0:0.25:1,0:0.25:1,'v1 decoding',legendStr_v1_dec)
legend('location','southeast')

fixPlot(subplot(325),[-0.6 0.6],[-0.6 0.6],'curvature dec error','beh error',-1:0.25:1,-1:0.25:1,'v4 decoding')
fixPlot(subplot(326),[-0.6 0.6],[-0.6 0.6],'curvature dec error','beh error',-1:0.25:1,-1:0.25:1,'v1 decoding')

ht = sgtitle(supp4_pca_data(sessId).name);
ht.Interpreter = 'none';
ht.Color = 'k'; ht.FontSize = 26; ht.FontName = 'Lato';

%% shape-general decoders
trialMarkerAl = 0.2; % draw trial dots in low alpha
sets = plotSets{sessId};
resp_v4 = supp4_pca_data(sessId).resp_v4;
resp_v1 = supp4_pca_data(sessId).resp_v1;

nPCs = 10;
red_v4 = pca(resp_v4','NumComponents',nPCs);
red_v1 = pca(resp_v1','NumComponents',nPCs);

% get nFold general decoders
nFold = 5;
nCurv = length(supp4_pca_data(sessId).curv);
beta_gen_v4 = nan(nFold,nPCs);
beta_gen_v1 = nan(nFold,nPCs);
for ff=1:nFold
    trainIdx = 1:nCurv;
    testIdx = randperm(nCurv,floor(nCurv/2));
    trainIdx(ismember(trainIdx,testIdx)) = [];
    beta_gen_v4(ff,:) = regress(supp4_pca_data(sessId).curv(trainIdx)',red_v4(trainIdx,:));
    beta_gen_v1(ff,:) = regress(supp4_pca_data(sessId).curv(trainIdx)',red_v1(trainIdx,:));
    idx_gen(ff,:) = testIdx;
end

cols = lines(length(sets));
legendStr_shape = cell(1,length(sets)*2);
legendStr_v4_dec = cell(1,length(sets)*2);
legendStr_v1_dec = cell(1,length(sets)*2);

figure('color','w','position',[-751,146,508,844]);
for ii=1:length(sets)
    idx_set = [supp4_pca_data(sessId).set] == sets(ii);

    cc = supp4_pca_data(sessId).curv(idx_set);
    ss = supp4_pca_data(sessId).sel(idx_set);

    red_v1_set = red_v1(idx_set,:);
    red_v4_set = red_v4(idx_set,:);

    % pca or tuning v4
    plotPCfn(subplot(321),red_v4_set,cc,cols(ii,:),'V4',trialMarkerAl);
    
    % pca v1
    plotPCfn(subplot(322),red_v1_set(:,1:3),cc,cols(ii,:),'V1',trialMarkerAl);

    % dec v4
    % testcc = nan(nFold,)
    dec_cc_v1 = []; dec_cc_v4 = [];
    testcc = []; testss = [];
    for ff=1:nFold
        testIdx = false(1,nCurv);
        testIdx(idx_gen(ff,:)) = true;
        testIdx = testIdx & idx_set;

        testcc = [testcc; supp4_pca_data(sessId).curv(testIdx)'];
        testss = [testss; supp4_pca_data(sessId).sel(testIdx)'];

        red_v1_set = red_v1(testIdx,:);
        red_v4_set = red_v4(testIdx,:);

        dec_cc_v1 = [dec_cc_v1; red_v1_set*beta_gen_v1(ff,:)'];
        dec_cc_v4 = [dec_cc_v4; red_v4_set*beta_gen_v4(ff,:)'];
    end
    plotTrendPatchLine(subplot(323),testcc',dec_cc_v4',cols(ii,:),[0 1],10);
    plotTrendPatchLine(subplot(324),testcc',dec_cc_v1',cols(ii,:),[0 1],10);

    % err v4
    cc_err = (dec_cc_v4 - testcc); cc_err = cc_err-mean(cc_err);
    ss_err = (dec_cc_v4 - testss); ss_err = ss_err-mean(ss_err);
    plotTrendPatchLine(subplot(325),cc_err',testss'-testcc',cols(ii,:),[-0.5 0.5],10);

    % err v1
    cc_err = (dec_cc_v1 - testcc); cc_err = cc_err-mean(cc_err);
    ss_err = (dec_cc_v1 - testss); ss_err = ss_err-mean(ss_err);
    plotTrendPatchLine(subplot(326),cc_err',testss'-testcc',cols(ii,:),[-0.5 0.5],10);

    % shape ids
    legendStr_shape{3*ii-2} = ''; legendStr_shape{3*ii-1} = '';
    legendStr_shape{ii*3} = ['shape ' num2str(sets(ii))];

    % dec acc v4
    legendStr_v4_dec{3*ii-2} = ''; legendStr_v4_dec{3*ii-1} = '';
    legendStr_v4_dec{ii*3} = ['acc = ' num2str(round(corr(dec_cc_v4,testcc),2))];

    % dec acc v1
    legendStr_v1_dec{3*ii-2} = ''; legendStr_v1_dec{3*ii-1} = '';
    legendStr_v1_dec{ii*3} = ['acc = ' num2str(round(corr(dec_cc_v1,testcc),2))];
end

fixPlot(subplot(323),[-0.1 1.1],[-0.1 1.1],'curvature','decoded curvature',0:0.25:1,0:0.25:1,'v4 decoding',legendStr_v4_dec)
legend('location','southeast')
fixPlot(subplot(324),[-0.1 1.1],[-0.1 1.1],'curvature','decoded curvature',0:0.25:1,0:0.25:1,'v1 decoding',legendStr_v1_dec)
legend('location','southeast')

fixPlot(subplot(325),[-0.6 0.6],[-0.6 0.6],'curvature dec error','beh error',-1:0.25:1,-1:0.25:1,'v4 decoding')
fixPlot(subplot(326),[-0.6 0.6],[-0.6 0.6],'curvature dec error','beh error',-1:0.25:1,-1:0.25:1,'v1 decoding')

ht = sgtitle(supp4_pca_data(sessId).name);
ht.Interpreter = 'none';
ht.Color = 'k'; ht.FontSize = 26; ht.FontName = 'Lato';


%% functions
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

function plotPCfn(h,redResp,param,colmap,titleStr,trialMarkerAl)
    [uParam,~,paramG] = unique(param);
    cols = [linspace(0.3,colmap(1),max(paramG))' linspace(0.3,colmap(2),max(paramG))' linspace(0.3,colmap(3),max(paramG))'];
    colsF = cols(paramG,:);

    mResp = groupsummary(redResp,paramG,'mean');
    respForFit = mResp;
    paramForFit = uParam;
    if trialMarkerAl > 0 % if you don't do this, then the axis gets scaled/zoomed out for no reason
        scatter3(h,redResp(:,1),redResp(:,2),redResp(:,3),50,colsF,'filled','MarkerFaceAlpha',trialMarkerAl); 
        respForFit = redResp;
        paramForFit = param;
    end
    hold(h,'on'); axis(h,'equal','square'); 
    
    scatter3(h,mResp(:,1),mResp(:,2),mResp(:,3),80,cols,'filled'); 

    % fit data parametrically, 3rd order polynomial
    t = paramForFit;
    polyOrder = [2 2 2];
    xcoefs=polyfit(t,respForFit(:,1),polyOrder(1)); x=polyval(xcoefs,t);
    ycoefs=polyfit(t,respForFit(:,2),polyOrder(2)); y=polyval(ycoefs,t);
    zcoefs=polyfit(t,respForFit(:,3),polyOrder(3)); z=polyval(zcoefs,t);
    plot3(h,x,y,z,'color',colmap,'linestyle','-','linewidth',2); 

    fix3dPlot(h,[],[],[],'PC1','PC2','PC3',[],[],[],titleStr)
    set(h,'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{})
    grid(h,'off')
end

