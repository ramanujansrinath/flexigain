clear; clf; close all
addpath('dep/')

nTrials = 2000;
uParams = initUnits_21dg();
sParams = rand(nTrials,2);
resp = getRespMat_21dg(uParams,sParams);
resp_red = pca(resp','NumComponents',3);
% figure('color','w','Position',[368,227,582,997])

cols = [linspace(253,31,nTrials)' linspace(244,181,nTrials)' linspace(104,236,nTrials)']/255;
[~,idx] = sort(sParams(:,1));
cols = cols.*repmat(0.3+sParams(idx,2)*0.7,1,3);

%% plot tuning functions
figure('color','w','Position',[368,705,985,519])
ha = arrayfun(@(ii) subplot(2,4,ii),1:4);
[ss,idx] = sort(sParams(:,1));
[~,uidx] = sort(uParams(:,2));
plotTuningFns(ha,uParams(uidx,1:3),ss,resp(idx,uidx))

ha = arrayfun(@(ii) subplot(2,4,ii),5:8);
[ss,idx] = sort(sParams(:,2));
[~,uidx] = sort(uParams(:,5));
plotTuningFns(ha,uParams(uidx,4:6),ss,resp(idx,uidx))

%% first 3 and next 3 PCs
figure('color','w','Position',[368,286,653,938])
rr = resp(idx,uidx);
redResp = pca(rr','NumComponents',5);
scatter3(subplot(321),redResp(:,1),redResp(:,2),redResp(:,3),20,'k','filled','MarkerFaceAlpha',1); 
set(gca,'xticklabel',{},'yticklabel',{},'zticklabel',{})
fix3dPlot(gca,[],[],[],'PC1','PC2','PC3',[],[],[],'PCA 1-3')
axis(gca,'square','equal')

scatter3(subplot(322),redResp(:,3),redResp(:,4),redResp(:,5),20,'k','filled','MarkerFaceAlpha',1); 
set(gca,'xticklabel',{},'yticklabel',{},'zticklabel',{})
fix3dPlot(gca,[],[],[],'PC3','PC4','PC5',[],[],[],'PCA 3-5')
axis(gca,'square','equal')

% PCA visualization of both features for PC1-3, and QR
subplot(323); scatter3(resp_red(:,1),resp_red(:,2),resp_red(:,3),20,sParams(:,1),'filled');
set(gca,'xticklabel',{},'yticklabel',{},'zticklabel',{})
fix3dPlot(gca,[],[],[],'PC1','PC2','PC3',[],[],[],'Feature 1: PCA 1-3')
axis(gca,'square','equal')
subplot(324); scatter3(resp_red(:,1),resp_red(:,2),resp_red(:,3),20,sParams(:,2),'filled');
set(gca,'xticklabel',{},'yticklabel',{},'zticklabel',{})
fix3dPlot(gca,[],[],[],'PC1','PC2','PC3',[],[],[],'Feature 2: PCA 1-3')
axis(gca,'square','equal')

beta_1 = regress(sParams(:,1),resp);
beta_2 = regress(sParams(:,2),resp);
[beta_orth,rr] = qr([beta_1 beta_2],'econ');
resp_qr = resp*beta_orth;
subplot(325); scatter(resp_qr(:,1),resp_qr(:,2),20,sParams(:,1),'filled');
fixPlot(gca,[],[],'feature 1','feature 2')
subplot(326); scatter(resp_qr(:,1),resp_qr(:,2),20,sParams(:,2),'filled');
fixPlot(gca,[],[],'feature 1','feature 2')


%% here, i assert the readout axis and learn the gains per task per neuron based on selectivity
trainIdx = sort(randperm(nTrials,nTrials*0.8)); testIdx = 1:nTrials; testIdx(trainIdx) = [];
pSel = (uParams(:,1)-uParams(:,4))./(uParams(:,1)+uParams(:,4));

% weights from the previous piece of code
% w_fixed = datasample([w_t1;w_t2],125,'Replace',false); save('dualTuningSim_w_fixed.mat','w_fixed');
% load('tuning_to_arc_sim.mat','w_fixed');

% random output weights?
% w_fixed = rand(125,1);

% assign fixed weights aligned with one of the variables
w_fixed = regress(sParams(trainIdx,1),resp(trainIdx,:));

% assign fixed weights from the weights shuffled from the axes aligned with both variables
% y = sParams(:,1); w_fixed_temp(:,1) = regress(y(trainIdx),resp(trainIdx,:));
% y = sParams(:,2); w_fixed_temp(:,2) = regress(y(trainIdx),resp(trainIdx,:));
% w_fixed = arrayfun(@(ii) w_fixed_temp(ii,randi(2)),1:125)';

tic
% task 1
% mult = double(pSel'>0); 
sacc_1 = -70+sParams(:,1)*140;
[task_1_gains,task_1_perf_train,task_1_perf_test] = ...
        getTaskGains(resp,sacc_1,pSel,w_fixed,trainIdx,testIdx);
sacc_pred_test_1 = resp(testIdx,:) * (task_1_gains .* w_fixed);
sacc_test_1 = sacc_1(testIdx);

toc
% task 2
% mult = double(pSel'<0); 
sacc_2 = -70+sParams(:,2)*140;
[task_2_gains,task_2_perf_train,task_2_perf_test] = ...
        getTaskGains(resp,sacc_2,-pSel,w_fixed,trainIdx,testIdx);
sacc_pred_test_2 = resp(testIdx,:) * (task_2_gains .* w_fixed);
sacc_test_2 = sacc_2(testIdx);
toc

%% shuffle gains for each task and check performance on test set and also plot gains
task_1_gains_shuff = shuffPerSel(task_1_gains,pSel);
task_2_gains_shuff = shuffPerSel(task_2_gains,pSel);

sacc_pred_test_1_shuff = resp(testIdx,:) * (task_1_gains_shuff .* w_fixed);
sacc_pred_test_2_shuff = resp(testIdx,:) * (task_2_gains_shuff .* w_fixed);

%% plot
gainRatio = (task_1_gains-task_2_gains)./(task_1_gains+task_2_gains);
gainRatio_shuf = (task_1_gains_shuff-task_2_gains_shuff)./(task_1_gains_shuff+task_2_gains_shuff);

figure('pos',[368,866,608,358],'color','w')
subplot(121);
h1 = plot(fitlm(pSel,gainRatio));
set(h1(1),'Marker','o','LineWidth',1.5,'MarkerFaceColor','w','color',[0 0 0]); hold on
set(h1(2),'LineWidth',2,'color',[0 0 0]);
set(h1(3),'LineWidth',1,'color',[0 0 0],'LineStyle','-');
set(h1(4),'LineWidth',1,'color',[0 0 0],'LineStyle','-');
h1 = plot(fitlm(pSel,gainRatio_shuf));
set(h1(1),'Marker','o','LineWidth',1.5,'MarkerFaceColor','w','color',[0.5 0.5 0.5]); hold on
set(h1(2),'LineWidth',2,'color',[0.5 0.5 0.5]);
set(h1(3),'LineWidth',1,'color',[0.5 0.5 0.5],'LineStyle','-');
set(h1(4),'LineWidth',1,'color',[0.5 0.5 0.5],'LineStyle','-');
legend off
fixPlot(gca,[-0.33 0.33],[-0.2 0.2],'parameter selectivity',{'task-dependent' 'response modulation'},-0.4:0.2:0.4,-2:0.1:2,'')

subplot(122);
yy = sacc_pred_test_1-mean(sacc_pred_test_1); yy = 70*yy./max(abs(yy));
plotTrendPatchLine(gca,sacc_test_1,yy,'k',[-70 70],21);
yy = sacc_pred_test_1_shuff-mean(sacc_pred_test_1_shuff); yy = 70*yy./max(abs(yy));
plotTrendPatchLine(gca,sacc_test_1,yy,[0.5 0.5 0.5],[-70 70],21);
fixPlot(gca,[-90 90],[-90 90],'correct saccade','learned saccade',-70:35:70,-70:35:70,'task 1 performance',{'' '' 'learned gains' '' '' 'shuffled gains'})

yy = sacc_pred_test_2-mean(sacc_pred_test_2); yy = 70*yy./max(abs(yy));
plotTrendPatchLine(gca,sacc_test_2,yy,'k',[-70 70],21);
yy = sacc_pred_test_2_shuff-mean(sacc_pred_test_2_shuff); yy = 70*yy./max(abs(yy));
plotTrendPatchLine(gca,sacc_test_2,yy,[0.5 0.5 0.5],[-70 70],21);
fixPlot(gca,[-90 90],[-90 90],'correct saccade','learned saccade',-70:35:70,-70:35:70,'task 2 performance',{'' '' 'learned gains' '' '' 'shuffled gains'})
legend('location','southeast')

%% decoding predictions
resp_gain_1 = resp.*repmat(task_1_gains',nTrials,1);
resp_gain_2 = resp.*repmat(task_2_gains',nTrials,1);
sacc_1 = -70+sParams(:,1)*140;
sacc_2 = -70+sParams(:,2)*140;

axis_1 = regress(sParams(trainIdx,1),resp_gain_1(trainIdx,:));
gain1_axis1 = [corr(sParams(testIdx,1),resp_gain_1(testIdx,:)*axis_1) ...
    corr(sParams(testIdx,2),resp_gain_1(testIdx,:)*axis_1) ...
    corr(sacc_1(testIdx),resp_gain_1(testIdx,:)*axis_1)];

axis_1 = regress(sParams(trainIdx,1),resp_gain_2(trainIdx,:));
gain2_axis1 = [corr(sParams(testIdx,1),resp_gain_2(testIdx,:)*axis_1) ...
    corr(sParams(testIdx,2),resp_gain_2(testIdx,:)*axis_1) ...
    corr(sacc_2(testIdx),resp_gain_2(testIdx,:)*axis_1)];

axis_2 = regress(sParams(trainIdx,2),resp_gain_1(trainIdx,:));
gain1_axis2 = [corr(sParams(testIdx,1),resp_gain_1(testIdx,:)*axis_2) ...
    corr(sParams(testIdx,2),resp_gain_1(testIdx,:)*axis_2) ...
    corr(sacc_1(testIdx),resp_gain_1(testIdx,:)*axis_2)];

axis_2 = regress(sParams(trainIdx,2),resp_gain_2(trainIdx,:));
gain2_axis2 = [corr(sParams(testIdx,1),resp_gain_2(testIdx,:)*axis_2) ...
    corr(sParams(testIdx,2),resp_gain_2(testIdx,:)*axis_2) ...
    corr(sacc_2(testIdx),resp_gain_2(testIdx,:)*axis_2)];

%%
function uParams = initUnits_21dg()
    % 5 levels of means of gaussian tuning for both parameter
    % [m1,m2] = meshgrid(linspace(0,1,5),linspace(0,1,5));
    % m1 = repmat(m1(:),5,1);
    % m2 = repmat(m2(:),5,1);

    % random means
    m1 = rand(125,1);
    m2 = rand(125,1);

    s1 = 0.6+0.3*rand(125,1);
    s2 = 0.6+0.3*rand(125,1);
    a1 = 0.5+0.5*rand(125,1);
    a2 = 0.5+0.5*rand(125,1);

    uParams = [a1 m1 s1 a2 m2 s2];
end

function resp = getRespMat_21dg(uParams,sParams)
    resp = nan(size(sParams,1),size(uParams,1));
    for uu=1:size(uParams,1)
        r1 = getGaussian(uParams(uu,1:3),sParams(:,1));
        r2 = getGaussian(uParams(uu,4:6),sParams(:,2));
        resp(:,uu) = mean([r1 r2],2);
        % resp(:,uu) = max([r1 r2],[],2);
    end
end

function num = getBoundedRandn_n(i_mean,i_std,lowLim,highLim)
    num = i_mean+i_std.*randn(length(i_mean),1);
    while sum(num > highLim | num < lowLim)
        idx = num > highLim | num < lowLim;
        num_i = i_mean(idx)+i_std(idx).*randn(sum(idx),1);
        num(idx) = num_i;
    end
end

function [task_gains,task_perf_train,task_perf_test] = ...
        getTaskGains(resp,sacc,pSel,w_fixed,trainIdx,testIdx)

    nFolds = 10000;
    cc = nan(1,nFolds);
    gains = nan(length(pSel),nFolds);
    X = resp(trainIdx,:);
    y = sacc(trainIdx);
    gainDist_mean = 1+0.05*(pSel./max(abs(pSel)));
    gainDist_std = 0.05*ones(length(pSel),1);
    for ff=1:nFolds
        % gg = pSel';
        % gg(pSel==1) = getBoundedRandn(1.05,0.5,sum(pSel),0.95,1.2);
        % gg(pSel==0) = getBoundedRandn(0.95,0.5,sum(~pSel),0.8,1.05);
        gg = getBoundedRandn_n(gainDist_mean,gainDist_std,0,2);
        % y_pred = X * (gg .* w_fixed);
        y_pred = (X .* repmat(gg',size(X,1),1)) * w_fixed;
        cc(ff) = corr(y,y_pred);
        gains(:,ff) = gg;
    end
    [task_perf_train,idx] = max(cc);
    task_gains = gains(:,idx);
    y_pred_test = resp(testIdx,:) * (task_gains .* w_fixed);
    task_perf_test = corr(sacc(testIdx),y_pred_test);
end

function task_gains_shuff = shuffPerSel(task_gains,pSel)
    [~,~,idx] = histcounts(pSel,linspace(-0.3,0.3,11));
    task_gains_shuff = task_gains;
    for ii=1:10
        gg = task_gains(idx==ii);
        task_gains_shuff(idx==ii) = gg(randperm(length(gg)));
    end
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

function plotTuningFns(ha,uParams,sParams,resp)
    nUnit = size(uParams,1);
    unitA = uParams(:,1);
    unitM = uParams(:,2);
    unitS = uParams(:,3);
    
    ss = linspace(0,1,100);
    % plot(ha(1),1:length(unitA),unitA,1:length(unitA),unitM,1:length(unitA),unitS)
    % axis(ha(1),[0 length(unitA)+1 -0.1 1.1])
    % set(ha(1),'xtick',[]); axis(ha(1),'square');
    histogram(unitA,linspace(0,1,21),'DisplayStyle','stairs','Parent',ha(1)); hold(ha(1),'on')
    histogram(unitM,linspace(0,1,21),'DisplayStyle','stairs','Parent',ha(1));
    histogram(unitS,linspace(0,1,21),'DisplayStyle','stairs','Parent',ha(1));
    fixPlot(ha(1),[0 1],[0 35],'paramter value','cell count',0:0.25:1,0:10:50,{'distribution of' 'simulated neuron params'},{'amp' 'pref' 'width'})
    legend('Location','northwest')
    for ii=1:nUnit
        % if mod(ii,nUnit/10)==0
        %     plot(ha(2),cc,getGaussian([unitA(ii) unitM(ii) unitS(ii)],cc),'k','linewidth',1);
        % else
            plot(ha(2),ss,getGaussian([unitA(ii) unitM(ii) unitS(ii)],ss),'color',[0.7 0.7 0.7],'linewidth',0.5);
        % end
        hold(ha(2),'on');
    end
    % axis(ha(2),'square','off')
    fixPlot(ha(2),[0 1],[0 1],'paramter value','firing rate',0:0.25:1,0:25:1,'tuning functions')
    
    imagesc(resp','parent',ha(3))
    % axis(ha(3),'square','off')
    fixPlot(ha(3),[0.5 2000.5],[0.5 125.5],'stimuli','neurons',0:500:2000,0:25:150,'responses')
    
    colspace = [zeros(256,1) linspace(0,1,256)' zeros(256,1)];
    redResp = pca(resp','NumComponents',3);
    scatter3(ha(4),redResp(:,1),redResp(:,2),redResp(:,3),50,sParams,'filled','MarkerFaceAlpha',1); 
    % scatter(redResp(:,1),redResp(:,2),50,stimC,'filled','MarkerFaceAlpha',0.3); 
    axis(ha(4),'square','equal')
    colormap(ha(4),colspace)
    set(ha(4),'xticklabel',{},'yticklabel',{},'zticklabel',{})
end