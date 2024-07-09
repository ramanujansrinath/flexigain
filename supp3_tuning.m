clc; close all; clear;
load('data/supp3_tuning.mat')

%%
figure('color','w','Position',[94,316,895,508])
% ====== random trial raster ================
randTrial = 291; % randi(length(params));
% set 14, curv 0.9074, sel 0.1853
hRaster = subplot(4,4,[1 2 5 6 9 10]);
hPsth = subplot(4,4,[13 14]);
plotRaster(hRaster,hPsth,params(randTrial),eid)

% ====== random unit tuning example v4 ======
sets = unique([params.set]); % sets = sets(2:3);
cols = [187 135 134; 142 230 232; 158 208 150; ]/255; % autumn(length(sets));

subplot(244); hold on; cla;
randUnitId = 46; % randi(size(v4sites_stim,2));
for ii=1:length(sets)
    cc = [params([params.set] == sets(ii)).curv];
    
    nnresp_stim = v4sites_stim([params.set] == sets(ii),randUnitId);
    % nnresp_stim = nnresp_stim-mean(respBase.v4sites(:,randUnitId));
    plotTrendPatchLine(gca,cc,nnresp_stim',cols(ii,:),true);

    nnresp_arc = v4sites_arc([params.set] == sets(ii),randUnitId);
    % nnresp_arc = nnresp_arc-mean(respBase.v4sites(:,randUnitId));
    plotTrendPatchLine(gca,cc,nnresp_arc',0.8*cols(ii,:),true);
end
fixPlot(gca,[-0.1 1.1],[30 270],'curvature','response',0:0.25:1,0:50:250,'example v4 unit'); 
% {'' 'shape 13' '' 'shape 13+arc' '' 'shape 14' '' 'shape 14+arc'})
% legend('box','on','location','northwest')

% ====== random unit tuning example v1 ======
subplot(243); hold on; cla;
randUnitId = 11; % randi(size(v1sites_stim,2));
for ii=1:length(sets)
    cc = [params([params.set] == sets(ii)).curv];
    
    nnresp_stim = v1sites_stim([params.set] == sets(ii),randUnitId);
    % nnresp_stim = nnresp_stim-mean(respBase.v1sites(:,randUnitId));
    plotTrendPatchLine(gca,cc,nnresp_stim',cols(ii,:),true);

    nnresp_arc = v1sites_arc([params.set] == sets(ii),randUnitId);
    % nnresp_arc = nnresp_arc-mean(respBase.v1sites(:,randUnitId));
    plotTrendPatchLine(gca,cc,nnresp_arc',0.8*cols(ii,:),true);
end
fixPlot(gca,[-0.1 1.1],[35 215],'curvature','response',0:0.25:1,0:50:250,'example v1 unit')
% legend('box','on','location','northwest')

% ====== selectivity for this session ======
sets = unique([params.set]); % sets = sets(2:3);
selId_v4_stim = nan(size(v4sites_stim,2),length(sets));
selId_v4_arc = nan(size(v4sites_stim,2),length(sets));
selId_v1_stim = nan(size(v1sites_stim,2),length(sets));
selId_v1_arc = nan(size(v1sites_stim,2),length(sets));
for ii=1:length(sets)
    cc = [params([params.set] == sets(ii)).curv];
    for nn=1:size(v4sites_stim,2)
        nnresp_stim = v4sites_stim([params.set] == sets(ii),nn);
        % nnresp_stim = nnresp_stim-mean(respBase.v4sites(:,nn));
        [~,nnresp_stim] = plotTrendPatchLine([],cc,nnresp_stim',[],true);
        selId_v4_stim(nn,ii) = (max(nnresp_stim)-min(nnresp_stim))./(max(nnresp_stim)+min(nnresp_stim));
        [nnresp_stim_max,idx] = max(nnresp_stim);
        
        nnresp_arc = v4sites_arc([params.set] == sets(ii),nn);
        % nnresp_arc = nnresp_arc-mean(respBase.v4sites(:,nn));
        [~,nnresp_arc] = plotTrendPatchLine([],cc,nnresp_arc',[],true);
        % selId_v4_arc(nn,ii) = (max(nnresp_arc)-min(nnresp_arc))./(max(nnresp_arc)+min(nnresp_arc));
        % selId_v4_arc(nn,ii) = (max(nnresp_arc)-max(nnresp_stim))./(max(nnresp_arc)+max(nnresp_stim));
        selId_v4_arc(nn,ii) = (nnresp_stim_max-nnresp_arc(idx))./(nnresp_stim_max+nnresp_arc(idx));
    end

    for nn=1:size(v1sites_stim,2)
        nnresp_stim = v1sites_stim([params.set] == sets(ii),nn);
        % nnresp_stim = nnresp_stim-mean(respBase.v1sites(:,nn));
        [~,nnresp_stim] = plotTrendPatchLine([],cc,nnresp_stim',[],true);
        selId_v1_stim(nn,ii) = (max(nnresp_stim)-min(nnresp_stim))./(max(nnresp_stim)+min(nnresp_stim));
        [nnresp_stim_max,idx] = max(nnresp_stim);
    
        nnresp_arc = v1sites_arc([params.set] == sets(ii),nn);
        % nnresp_arc = nnresp_arc-mean(respBase.v1sites(:,nn));
        [~,nnresp_arc] = plotTrendPatchLine([],cc,nnresp_arc',[],true);
        % selId_v1_arc(nn,ii) = (max(nnresp_arc)-min(nnresp_arc))./(max(nnresp_arc)+min(nnresp_arc));
        % selId_v1_arc(nn,ii) = (max(nnresp_arc)-max(nnresp_stim))./(max(nnresp_arc)+max(nnresp_stim));
        selId_v1_arc(nn,ii) = (nnresp_stim_max-nnresp_arc(idx))./(nnresp_stim_max+nnresp_arc(idx));
    end
end
subplot(247); hold on;
% plot(selId_v1_stim(:),abs(selId_v1_arc(:)),'.','color',[0.2 0.5 0.9],'markersize',12)
% plot(selId_v4_stim(:),abs(selId_v4_arc(:)),'.','color',[0.9 0.5 0.2],'markersize',12)
scatter(abs(selId_v1_stim(:)),abs(selId_v1_arc(:)),20,[0.2 0.5 0.9],'filled','MarkerFaceAlpha',1)
fixPlot(gca,[-0.04 0.44],[-0.04 0.44],'curvature selectivity','curvature+arc selectivity',0:0.2:1,0:0.2:1,'v1')
subplot(248); hold on;
scatter(abs(selId_v4_stim(:)),abs(selId_v4_arc(:)),20,[0.9 0.5 0.2],'filled','MarkerFaceAlpha',1)
fixPlot(gca,[-0.04 0.44],[-0.04 0.44],'curvature selectivity','curvature+arc selectivity',0:0.2:1,0:0.2:1,'v4')

%% selectivity across sessions
load('data/supp3_selectivity.mat','selId_v4_arc_all','selId_v1_arc_all','selId_v4_stim_all','selId_v1_stim_all');

idx = (selId_v1_stim_all==1 | selId_v1_arc_all==1);
selId_v1_stim_all(idx) = [];
selId_v1_arc_all(idx) = [];
idx = (selId_v4_stim_all==1 | selId_v4_arc_all==1);
selId_v4_stim_all(idx) = [];
selId_v4_arc_all(idx) = [];

figure('color','w','pos',[390,410,794,329])
subplot(245); hold on;
plot(selId_v1_stim_all,selId_v1_arc_all,'.','color',[0.2 0.5 0.9],'markersize',10)
plot(nanmean(selId_v1_stim_all),nanmean(selId_v1_arc_all),'+','markersize',15,'linewidth',1,'color',0.3*[0.2 0.5 0.9]);
line([0 1],[0 1],'linestyle','--','color','k','linewidth',2);

beta = regress(selId_v1_arc_all,[ones(size(selId_v1_stim_all,1),1) selId_v1_stim_all]);
line([0 1],[beta(1) sum(beta)],'linestyle','--','color','r','linewidth',2);

fixPlot(gca,[0 1],[-0.5 1],'curvature selectivity','curvature+arc selectivity',0:0.5:1,-0.5:0.5:1,'v1')

% binE = -0.0125:0.025:1.0125;
% binC = (binE+circshift(binE,-1))/2; binC = binC(1:end-1);
% [cx,cy] = meshgrid(binC,binC);
% clevels = 50:100:1000;
% a = histcounts2(selId_v1_stim_all,selId_v1_arc_all,binE,binE);
% [~,ha] = contour(cx,cy,a,clevels);

subplot(241);
histogram(selId_v1_stim_all,linspace(0,1,25),'DisplayStyle','stairs','EdgeColor',[0.2 0.5 0.9],'linewidth',2)
fixPlot(gca,[0 1],[0 12000],'','',0:0.5:1,[])

subplot(246);
histogram(selId_v1_arc_all,linspace(-0.5,1,25),'DisplayStyle','stairs','EdgeColor',[0.2 0.5 0.9],'linewidth',2)
fixPlot(gca,[-0.5 1],[0 12000],'','',-0.5:0.5:1,[])

subplot(247); hold on;
plot(selId_v4_stim_all,selId_v4_arc_all,'.','color',[0.9 0.5 0.2],'markersize',10)
plot(nanmean(selId_v4_stim_all),nanmean(selId_v4_arc_all),'+','markersize',15,'linewidth',1,'color',0.3*[0.9 0.5 0.2]);
line([0 1],[0 1],'linestyle','--','color','k','linewidth',2);
beta = regress(selId_v4_arc_all,[ones(size(selId_v4_stim_all,1),1) selId_v4_stim_all]);
line([0 1],[beta(1) sum(beta)],'linestyle','--','color','r','linewidth',2);
fixPlot(gca,[0 1],[-0.5 1],'curvature selectivity','curvature+arc selectivity',0:0.5:1,-0.5:0.5:1,'v4')

% a = histcounts2(abs(selId_v4_stim_all),abs(selId_v4_arc_all),binE,binE);
% contour(cx,cy,a,clevels)

subplot(243);
histogram(selId_v4_stim_all,linspace(0,1,25),'DisplayStyle','stairs','EdgeColor',[0.9 0.5 0.2],'linewidth',2)
fixPlot(gca,[0 1],[0 12000],'','',0:0.5:1,[])

subplot(248);
histogram(selId_v4_arc_all,linspace(-0.5,1,25),'DisplayStyle','stairs','EdgeColor',[0.9 0.5 0.2],'linewidth',2)
fixPlot(gca,[-0.5 1],[0 12000],'','',-0.5:0.5:1,[])

%% functions
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

function plotRaster(hRaster,hPsth,param_trial,eid)
    tMarks = [0 param_trial.time.stimOn param_trial.time.fixStayAfterTargOn param_trial.time.targStay param_trial.time.corrStay];
    tMarks = [-param_trial.time.fixStay  cumsum(tMarks)];
    
    cols = [0.9 0.5 0.2; 0.2 0.5 0.9];
    hold(hRaster,'on');
    
    gSites = find(eid(:,1) == 1);
    pl = ismember(param_trial.spikes(:,2),gSites);
    v4s = param_trial.spikes(pl,:); 
    [~,~,a] = unique(v4s(:,2)); v4s(:,2) = a;
    % plot(v4s(:,1),v4s(:,2),'b.','color',cols(1,:))
    scatter(hRaster,v4s(:,1),v4s(:,2),7,cols(1,:),'filled','markerfacealpha',1);

    gSites = find(eid(:,1) == 2);
    pl = ismember(param_trial.spikes(:,2),gSites);
    v1s = param_trial.spikes(pl,:);
    [~,~,a] = unique(v1s(:,2)); v1s(:,2) = a+max(v4s(:,2));
    scatter(hRaster,v1s(:,1),v1s(:,2),7,cols(2,:),'filled','markerfacealpha',1);
    
    arrayfun(@(jj) line(hRaster,[tMarks(jj) tMarks(jj)],[1 256],'linewidth',2,'color','k'),1:length(tMarks))
    tMax = round(max(tMarks),-2)+200;
    tMax = 950;
    
    % t = [floor(min([v4s(:,1);v1s(:,1);sas(:,1)])) ceil(max([v4s(:,1);v1s(:,1);sas(:,1)]))];
    t = (param_trial.lims_trial - param_trial.lims_stim(1))/30;
    tLims = [t(1) tMax];
    fixPlot(hRaster,tLims,[0 193],'','channel',0:250:t(2),0:96:256) % 257 for 7a
    set(hRaster,'xtick',[]);
    axis(hRaster,'normal'); grid(hRaster,'off')

    % PSTH
    hold(hPsth,'on')
    allS = [v1s;v4s];
    % t = floor(min(allS(:,1))) :1: ceil(max(allS(:,1)));
    t = linspace(tLims(1),tLims(2),200);
    psths = zeros(192,length(t));
    for ii=1:192
        sp = allS(allS(:,2)==ii,1);
        
        psth = zeros(1,length(t));
        for s=1:length(sp)
            temp1 = getGaussian([1,sp(s),20],t);
            temp2 = getGaussian([1,sp(s),90],t);
            [~,idx] = max(temp1);
            temp = [temp1(1:idx) temp2(idx+1:end)];
            psth = psth + temp;
        end
        psth = psth-mean(psth(t<0));
        psth = psth * 1;
        % if ii<=96
        %     col = cols(2,:);
        % else
        %     col = cols(1,:);
        % end
        % plot(hPsth,t,psth,'color',col); hold on;
        psths(ii,:) = psth;
    end
    psth_v1 = mean(psths(1:96,:));
    psth_v4 = mean(psths(97:end,:));
    plot(hPsth,t,psth_v1,'color',cols(2,:),'linewidth',2);
    plot(hPsth,t,psth_v4,'color',cols(1,:),'linewidth',2);
    line(hPsth,[50 550],[0 0],'linewidth',5,'color','k')
    arrayfun(@(jj) line(hPsth,[tMarks(jj) tMarks(jj)],[-1 10],'linewidth',2,'color','k'),1:length(tMarks))
    fixPlot(hPsth,tLims,[-1 10],'time (ms)','response (a.u.)',0:250:max(t),0:96:256)
    axis(hPsth,'normal'); grid(hRaster,'off')
    % psths = flipud(psths);
end
