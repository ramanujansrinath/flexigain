clear; clc;
addpath('dep/')

nTrials = 60;
uParams = initUnits();
sParams(:,1) = repmat(linspace(0,1,nTrials/2)',2,1);
sParams(:,2) = [0.2*ones(nTrials/2,1);0.8*ones(nTrials/2,1)];

resp = getRespMat(uParams,sParams);
resp_red = pca(resp','NumComponents',3);
% resp_red = resp_red(:,[1 3]); % for 2d pca plots

%% plot the tuning functions
figure('color','w','Position',[368,778,985,446])
ha = arrayfun(@(ii) subplot(2,4,ii),1:4);
plotTuningFns(ha,uParams(:,1:3),sParams(:,1),resp)

ha = arrayfun(@(ii) subplot(2,4,ii),5:8);
plotTuningFns(ha,uParams(:,[1 4 5]),sParams(:,2),resp)

%% first, plot the space
cols1 = [linspace(0.3,0.9,nTrials/2)' linspace(0.3,0.5,nTrials/2)' linspace(0.3,0.3,nTrials/2)'];
cols2 = [linspace(0.3,0.3,nTrials/2)' linspace(0.3,0.5,nTrials/2)' linspace(0.3,0.9,nTrials/2)'];
cols = [cols1;cols2];

clf; set(gcf,'color','w','Position',[368,840,951,384])
subplot(131);
% scatter3(resp_red(:,1),resp_red(:,2),resp_red(:,3),100,cols,'filled'); axis equal
scatter(resp_red(:,1),resp_red(:,2),35,cols,'filled'); axis equal
axis square; hold on;

%% now, general and specific decoders
ll_g = plotPCfn_2d(gca,resp_red,sParams(:,1),[0.3 0.3 0.3],1);
sets = unique(sParams(:,2));

cols = [0.9 0.5 0.3; 0.3 0.5 0.9];
for ii=1:length(sets)
    cc = sParams(sParams(:,2)==sets(ii),1);
    resp_red_set = resp_red(sParams(:,2)==sets(ii),:);

    ll = plotPCfn_2d(gca,resp_red_set,cc,cols(ii,:),1);
end
% if 2d plot only
fixPlot(gca,[-0.2 0.2],[-0.3 0.3],'pc1','pc2',[],[])

%% assert the readout axis and learn the gains per task per neuron
clearvars sParams
nTrials = 500;
sParams(:,1) = repmat(linspace(0,1,nTrials/2)',2,1);
sParams(:,2) = [0.2*ones(nTrials/2,1);0.8*ones(nTrials/2,1)];
resp = getRespMat(uParams,sParams);
% [~,~,~,~,ex] = pca(resp,'NumComponents',20); ex = ex(1:20); plot(ex); hold on; find(cumsum(ex)>95,1,'first')-1
% [~,~,~,~,ex] = pca(resp.*repmat(task_1_gains',500,1),'NumComponents',20); ex = ex(1:20); plot(ex); find(cumsum(ex)>95,1,'first')-1
% [~,~,~,~,ex] = pca(resp.*repmat(task_2_gains',500,1),'NumComponents',20); ex = ex(1:20); plot(ex); find(cumsum(ex)>95,1,'first')-1

trainIdx = sort(randperm(nTrials,nTrials*0.8)); testIdx = 1:nTrials; testIdx(trainIdx) = [];

% you can pick random weights, weights from a different simulation, or, as
% in here, you can make it easier to find weights from arbitrary
% distributions by asserting that the fixed readout is aligned to one of
% the arc representations
w_fixed = regress(sParams(trainIdx,1),resp(trainIdx,:));

% arc 1
sacc = -90+sParams(:,1)*140;
[task_1_gains,task_1_perf_train,task_1_perf_test] = ...
        getTaskGains(resp,sacc,w_fixed,trainIdx,testIdx);
sacc_pred_test_1 = resp(testIdx,:) * (task_1_gains .* w_fixed);
sacc_test_1 = sacc(testIdx);

% arc 2
sacc = -50+sParams(:,1)*140;
[task_2_gains,task_2_perf_train,task_2_perf_test] = ...
        getTaskGains(resp,sacc,w_fixed,trainIdx,testIdx);
sacc_pred_test_2 = resp(testIdx,:) * (task_2_gains .* w_fixed);
sacc_test_2 = sacc(testIdx);

% shuffle gains for each task and check performance on test set and also plot gains
task_1_gains_shuff = task_1_gains(randperm(length(task_1_gains)));
task_2_gains_shuff = task_2_gains(randperm(length(task_2_gains)));

sacc_pred_test_1_shuff = resp(testIdx,:) * (task_1_gains_shuff .* w_fixed);
sacc_pred_test_2_shuff = resp(testIdx,:) * (task_2_gains_shuff .* w_fixed);

% for rescaling
minSacc = resp(1,:) * (task_1_gains .* w_fixed);
maxSacc = resp(end,:) * (task_2_gains .* w_fixed);
meanSacc = (minSacc+maxSacc)/2;

% now plot ratio of gains and task performance
gainRatio = (task_1_gains-task_2_gains)./(task_1_gains+task_2_gains);
gainRatio_shuf = (task_1_gains_shuff-task_2_gains_shuff)./(task_1_gains_shuff+task_2_gains_shuff);
subplot(132); cla; hold on
histogram(gainRatio,linspace(-0.2,0.2,21),'DisplayStyle','stairs','edgecolor','k','linewidth',2)
histogram(gainRatio_shuf,linspace(-0.2,0.2,20),'DisplayStyle','stairs','edgecolor',[0.5 0.5 0.5],'linewidth',2)
fixPlot(gca,[-0.25 0.25],[0 50],{'task-dependent' 'response modulation'},'# neurons',-0.2:0.1:0.2,0:20:100,'',{'learned','shuffled'})
legend('Location','northwest')

subplot(133); cla; hold on
yy = sacc_pred_test_1-meanSacc; yy = 140*yy./(maxSacc-minSacc);
plotTrendPatchLine(gca,sacc_test_1,yy,'b',[-90 50],21);

yy = sacc_pred_test_2-meanSacc; yy = 140*yy./(maxSacc-minSacc);
plotTrendPatchLine(gca,sacc_test_2,yy,'g',[-50 90],21);

yy = sacc_pred_test_1_shuff-meanSacc; yy = 140*yy./(maxSacc-minSacc);
plotTrendPatchLine(gca,sacc_test_1,yy,[0.5 0.5 0.5],[-90 50],21);
% yy = sacc_pred_test_2_shuff-meanSacc; yy = 70*yy./(maxSacc-minSacc);
% plotTrendPatchLine(gca,sacc_test_2,yy,[0.5 0.5 0.5],[-50 90],21);
fixPlot(gca,[-100 100],[-100 100],'correct saccade','learned saccade',-70:35:70,-70:35:70,'',{'' '' 'arc 1' '' '' 'arc 2' '' '' 'shuffle'})


%% helpers
function uParams = initUnits()
    % 5 levels of means of gaussian tuning for both parameter
    [m1,m2] = meshgrid(linspace(0,1,5),linspace(0,1,5));
    m1 = repmat(m1(:),5,1);
    m2 = repmat(m2(:),5,1);

    s1 = 0.3+0.4*rand(125,1);
    s2 = 0.3+0.4*rand(125,1);
    a = 0.5+0.5*rand(125,1);

    uParams = [a m1 s1 m2 s2];
end

function resp = getRespMat(uParams,sParams)
    resp = nan(size(sParams,1),size(uParams,1));
    for uu=1:size(uParams,1)
        resp(:,uu) = getGaussian2d(uParams(uu,:),sParams(:,1),sParams(:,2));
    end
end


function num = getBoundedRandn(i_mean,i_std,lowLim,highLim)
    num = i_mean+i_std.*randn(length(i_mean),1);
    while sum(num > highLim | num < lowLim)
        idx = num > highLim | num < lowLim;
        num_i = i_mean(idx)+i_std(idx).*randn(sum(idx),1);
        num(idx) = num_i;
    end
end

function [task_gains,task_perf_train,task_perf_test] = ...
        getTaskGains(resp,sacc,w_fixed,trainIdx,testIdx)

    nFold = 100000;
    nUnit = size(resp,2);
    cc = nan(1,nFold);
    gains = nan(nUnit,nFold);
    X = resp(trainIdx,:);
    y = sacc(trainIdx);
    gainDist_mean = ones(nUnit,1);
    gainDist_std = 0.05*ones(nUnit,1);
    parfor ff=1:nFold
        gg = getBoundedRandn(gainDist_mean,gainDist_std,0,2);
        y_pred = (X .* repmat(gg',size(X,1),1)) * w_fixed;
        cc(ff) = corr(y,y_pred);
        gains(:,ff) = gg;
    end
    [task_perf_train,idx] = max(cc);
    task_gains = gains(:,idx);
    y_pred_test = resp(testIdx,:) * (task_gains .* w_fixed);
    task_perf_test = corr(sacc(testIdx),y_pred_test);
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

function fit = plotPCfn_2d(h,redResp,param,colmap,polyOrder)
    xcoefs=polyfit(param,redResp(:,1),polyOrder); x=polyval(xcoefs,[0 1]);
    ycoefs=polyfit(param,redResp(:,2),polyOrder); y=polyval(ycoefs,[0 1]);
    plot(h,x,y,'color',colmap,'linestyle','-','linewidth',2); 
    fit = [x;y]';
end

function plotTuningFns(ha,uParams,sParams,resp)
    nUnit = size(uParams,1);
    unitA = uParams(:,1);
    unitM = uParams(:,2);
    unitS = uParams(:,3);
    
    ss = linspace(0,1,100);
    plot(ha(1),1:length(unitA),unitA,1:length(unitA),unitM,1:length(unitA),unitS)
    axis(ha(1),[0 length(unitA)+1 -0.1 1.1])
    set(ha(1),'xtick',[]); axis(ha(1),'square');
    for ii=1:nUnit
        % if mod(ii,nUnit/10)==0
        %     plot(ha(2),cc,getGaussian([unitA(ii) unitM(ii) unitS(ii)],cc),'k','linewidth',1);
        % else
            plot(ha(2),ss,getGaussian([unitA(ii) unitM(ii) unitS(ii)],ss),'color',[0.7 0.7 0.7],'linewidth',0.5);
        % end
        hold(ha(2),'on');
    end
    axis(ha(2),'square','off')
    
    imagesc(resp','parent',ha(3))
    axis(ha(3),'square','off')
    
    colspace = [zeros(256,1) linspace(0,1,256)' zeros(256,1)];
    redResp = pca(resp','NumComponents',3);
    scatter3(ha(4),redResp(:,1),redResp(:,2),redResp(:,3),50,sParams,'filled','MarkerFaceAlpha',1); 
    % scatter(redResp(:,1),redResp(:,2),50,stimC,'filled','MarkerFaceAlpha',0.3); 
    axis(ha(4),'square','equal')
    colormap(ha(4),colspace)
    set(ha(4),'xticklabel',{},'yticklabel',{},'zticklabel',{})
end