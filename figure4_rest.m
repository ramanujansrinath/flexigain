clear; close all; clc
figure('pos',[285,722,1158,522],'color','w');

%% plot behavior and representation for an example session
load('data/figure4_example','params','resp','resp_base')
plot_psychometric(subplot(141),params)
plot_PCA_QR([subplot(243) subplot(244) subplot(247) subplot(248)],params,resp,resp_base)

%% plot gains vs selectivity across sessions
load('data/figure4_gains.mat','pSel','task_gains')
% here, i'm loading the gains and selectivity for each session because i
% would otherwise have to upload all the responses. that
% said, the code to calculate those two values for each unit for each
% session is below ("function getGains")
plot_gainsSelectivity(subplot(142),pSel,task_gains)

%%
function plot_psychometric(h,params)
    colorDiff = (mod([params.stimRF_num],5)-mod([params.stimOpp_num],5)) == 0;
    stimVals = [mat2vec(repmat(1:5,5,1)) mat2vec(repmat(5:-1:1,1,5))];
    colIds_1 = stimVals([params.stimRF_num],1);
    colIds_2 = stimVals([params.stimOpp_num],1);
    
    shpIds_1 = stimVals([params.stimRF_num],2);
    shpIds_2 = stimVals([params.stimOpp_num],2);
    choices = [params.selected]';
    
    hold(h,'on');
    diff_col = colIds_1(colorDiff)-colIds_2(colorDiff);
    xx = unique(diff_col);
    mm = groupsummary(choices(colorDiff),diff_col,'mean')-1;
    ss = groupsummary(choices(colorDiff),diff_col,'std')./sqrt(groupsummary(choices(colorDiff),diff_col,'nnz'));
    patch(h,[xx;flipud(xx)],[mm-ss/2;flipud(mm+ss/2)],[0.9 0.5 0.2],'edgecolor','none','facealpha',0.5)
    plot(h,xx,mm,'-o','color',[0.9 0.5 0.2],'LineWidth',2,'MarkerFaceColor','w');
    
    diff_shp = shpIds_1(~colorDiff)-shpIds_2(~colorDiff);
    xx = unique(diff_shp);
    mm = groupsummary(choices(~colorDiff),diff_shp,'mean')-1;
    ss = groupsummary(choices(~colorDiff),diff_shp,'std')./sqrt(groupsummary(choices(~colorDiff),diff_shp,'nnz'));
    patch(h,[xx;flipud(xx)],[mm-ss/2;flipud(mm+ss/2)],[0.2 0.5 0.9],'edgecolor','none','facealpha',0.5)
    plot(h,xx,mm,'-o','color',[0.2 0.5 0.9],'LineWidth',2,'MarkerFaceColor','w');
    
    fixPlot(h,[-4.8 4.8],[-0.1 1.1],'stim difference','probability of saccde to RF',-4:4,0:0.25:1)
end


function plot_PCA_QR(ha,params,resp,resp_base)
    eid = [3*ones(1,32) ones(1,64) 3*ones(1,32)];
    resp = resp(:,eid==1);
    resp_base = resp_base(:,eid==1);
    good = mean(resp)>1.1*mean(resp_base); % true(1,size(resp,2));
    resp = resp(:,good);
    
    stimVals = [mat2vec(repmat(1:5,5,1)) mat2vec(repmat(1:5,1,5))];
    colIds = stimVals([params.stimRF_num],1);
    shpIds = stimVals([params.stimRF_num],2);
    
    % PCA
    resp_pca = pca(resp','NumComponents',2);
    
    scatter(ha(1),resp_pca(:,1),resp_pca(:,2),20,colIds,'filled'); title('pca color')
    fixPlot(ha(1),[],[],'pc1','pc2')
    colormap(ha(1),round([linspace(247,87,5)' linspace(240,181,5)' linspace(107,229,5)'])/255)
    
    scatter(ha(2),resp_pca(:,1),resp_pca(:,2),20,shpIds,'filled'); title('pca shape')
    fixPlot(ha(2),[],[],'pc1','pc2')
    colormap(ha(2),round([linspace(62,188,5)' linspace(62,188,5)' linspace(62,188,5)'])/255)
    
    % QR
    beta_col = regress(colIds,resp);
    beta_shp = regress(shpIds,resp);
    [beta_orth,rr] = qr([beta_col beta_shp],'econ'); %#ok<ASGLU>
    resp_qr = resp*beta_orth;
    resp_qr = resp_qr - repmat(mean(resp_qr),size(resp_qr,1),1);
    resp_qr = resp_qr./max(resp_qr(:));
    
    scatter(ha(3),resp_qr(:,1),resp_qr(:,2),20,colIds,'filled'); title('qr color')
    fixPlot(ha(3),[-1 1],[-1 1],'color axis','shape axis',-1:0.5:1,-1:0.5:1)
    colormap(ha(3), round([linspace(247,87,5)' linspace(240,181,5)' linspace(107,229,5)'])/255)
    
    scatter(ha(4),resp_qr(:,1),resp_qr(:,2),20,shpIds,'filled'); title('qr shape')
    fixPlot(ha(4),[-1 1],[-1 1],'color axis','shape axis',-1:0.5:1,-1:0.5:1)
    colormap(ha(4),round([linspace(62,188,5)' linspace(62,188,5)' linspace(62,188,5)'])/255)
end

% this is the function to calculate 
% - response modulation ("gain"): the modulation ratio of max vs min
% response between two tasks for both features
% - feature selectivity: the modulation ratio of the max vs min response
% between two features for both tasks
function [pSel,task_gains] = getGains(params,resp,resp_base) %#ok<DEFNU>
    colorDiff = (mod([params.stimRF_num],5)-mod([params.stimOpp_num],5)) == 0;
    eid = [3*ones(1,32) ones(1,64) 3*ones(1,32)];
    resp = resp(:,eid==1);
    resp_base = resp_base(:,eid==1);
    good = mean(resp)>1.1*mean(resp_base); % true(1,size(resp,2));
    resp = resp(:,good);
    
    stimVals = [mat2vec(repmat(1:5,5,1)) mat2vec(repmat(1:5,1,5))];
    colIds = stimVals([params.stimRF_num],1);
    shpIds = stimVals([params.stimRF_num],2);
    
    pSel = nan(1,size(resp,2));
    task_gains = nan(1,size(resp,2));
    for uu=1:size(resp,2)
        a = max(groupsummary(resp(:,uu),colIds,'mean'));
        b = max(groupsummary(resp(:,uu),shpIds,'mean'));
        pSel(uu) = (a-b)./(a+b);
        
        a = max(groupsummary(resp(colorDiff,uu),colIds(colorDiff),'mean'));
        b = max(groupsummary(resp(~colorDiff,uu),colIds(~colorDiff),'mean'));
        % b = max(groupsummary(resp(~colorDiff,uu),colIds(~colorDiff),'mean'));
        task_gains(uu) = (a-b)./(a+b);
    end
end

function plot_gainsSelectivity(h,pSel,task_gains)
    pSel = cell2mat(pSel);
    task_gains = cell2mat(task_gains);

    h1 = plot(h,fitlm(pSel,task_gains));
    set(h1(1),'Marker','o','LineWidth',1.5,'MarkerFaceColor','w','color',[0 0 0]); hold on
    set(h1(2),'LineWidth',2,'color',[0 0 0]);
    set(h1(3),'LineWidth',1,'color',[0 0 0],'LineStyle','-');
    set(h1(4),'LineWidth',1,'color',[0 0 0],'LineStyle','-');
    fixPlot(h,[-0.2 0.2],[-0.2 0.2],'parameter selectivity',{'task-dependent' 'response modulation'},-0.4:0.05:0.4,-2:0.05:2,'',{})
end