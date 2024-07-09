clc; clear
load('data/supp8_rnn.mat')

%%
figure('color','w','pos',[476,521,891,345])
params = cell2mat(cellfun(@(x) [x.s double(x.b) double(x.a) double(x.o)]',trial_params,'UniformOutput',false))';
arc_conds = unique(params(:,2:3),'rows');
cols = [copper(5);winter(3)]; cols(1:2,:) = []; % cols = reshape(cols',3,6)';
clf;
for ii=1:size(arc_conds,1)
    curv = params(params(:,2)==arc_conds(ii,1) & params(:,3)==arc_conds(ii,2),1);
    correct_sacc = params(params(:,2)==arc_conds(ii,1) & params(:,3)==arc_conds(ii,2),4);
    actual_sacc = outputs(:,end);
    actual_sacc = actual_sacc(params(:,2)==arc_conds(ii,1) & params(:,3)==arc_conds(ii,2));

    subplot(132); hold on;
    scatter(curv,correct_sacc,12,repmat(cols(ii,:),length(curv),1).*(0.5+0.5*repmat(curv,1,3)),'filled');

    subplot(133); hold on;
    plotTrendPatchLine(gca,curv,actual_sacc,cols(ii,:),[0 1],11);
    scatter(curv,actual_sacc,12,repmat(cols(ii,:),length(curv),1).*(0.5+0.5*repmat(curv,1,3)),'filled');
end
fixPlot(subplot(132),[-0.1 1.1],[-90 90],'curvature','saccade target',0:0.25:1,-90:45:90,'inference to expected behaviour')
fixPlot(subplot(133),[-0.1 1.1],[-90 90],'curvature','actual saccade',0:0.25:1,-90:45:90,'inference to actual behaviour')

%%
% relevant time points
son = round(mean(cellfun(@(x) x.s_onset,trial_params))/10);
abon = round(mean(cellfun(@(x) x.ab_onset,trial_params))/10);
smid = (son+abon)/2;
fixoff = round(mean(cellfun(@(x) x.fix_offset,trial_params))/10);
abmid = (abon+fixoff)/2;

nCurvBin = 8;
curv = cellfun(@(x) x.s,trial_params);
targ = cellfun(@(x) x.o,trial_params);
binE = linspace(0,1,nCurvBin+1);
binC = (binE+circshift(binE,-1))/2; binC = binC(1:end-1);

% qr
dimIDs = [1 2 3]; nComp = 3;
[~,~,aConds] = unique(params(:,2:3),'rows');
beta_curv = regress(curv',squeeze(state_vars(:,fixoff,:)));
beta_acond = regress(aConds,squeeze(state_vars(:,fixoff,:)));
beta_sacc = regress(outputs(:,fixoff),squeeze(state_vars(:,fixoff,:)));
[beta_orth,rr] = qr([beta_curv beta_acond beta_sacc],'econ');
traj = cell2mat(arrayfun(@(ii) reshape(squeeze(state_vars(:,ii,:))*beta_orth,120,1,nComp),1:180,'UniformOutput',false));

% pca
% nComp = 10; dimIDs = [2 3 4];
% [~,beta] = pca(squeeze(state_vars(:,fixoff,:))','NumComponents',nComp);
% traj = cell2mat(arrayfun(@(ii) reshape(squeeze(state_vars(:,ii,:))*beta,1000,1,nComp),1:180,'UniformOutput',false));
% traj = reshape(pca(reshape(state_vars,[1000*180 50])','NumComponents',nComp),[1000 180 nComp]);

traj_ii_jj_cc = nan(size(traj,2),nComp,6,nCurvBin);
subplot(131); hold on;
for ii=1:size(arc_conds,1)
    curv_ii = params(params(:,2)==arc_conds(ii,1) & params(:,3)==arc_conds(ii,2),1);
    traj_ii = traj(params(:,2)==arc_conds(ii,1) & params(:,3)==arc_conds(ii,2),:,:);
    traj_ii = traj_ii(:,:,dimIDs);

    % plot after binning by curvatures
    [~,~,idx] = histcounts(curv_ii,binE);
    for jj=1:length(binC)
        traj_jj = squeeze(mean(traj_ii(idx==jj,:,:)));
        % plot3(traj_jj(:,1),traj_jj(:,2),traj_jj(:,3),...
        %     'color',binC(jj)*cols(ii,:));
        % 
        % plot3(traj_jj(son,1),traj_jj(son,2),traj_jj(son,3),...
        %     'o-','color',binC(jj)*cols(ii,:));
        % plot3(traj_jj(smid,1),traj_jj(smid,2),traj_jj(smid,3),...
        %     'o-','color',binC(jj)*cols(ii,:));
        % plot3(traj_jj(abon,1),traj_jj(abon,2),traj_jj(abon,3),...
        %     'o-','color',binC(jj)*cols(ii,:));
        % plot3(traj_jj(fixoff,1),traj_jj(fixoff,2),traj_jj(fixoff,3),...
        %     'o-','color',binC(jj)*cols(ii,:));

        traj_ii_jj_cc(:,:,ii,jj) = traj_jj;
    end

    if ii==1
        plot3(squeeze(traj_ii_jj_cc(smid,1,ii,:)),...
              squeeze(traj_ii_jj_cc(smid,2,ii,:)),...
              squeeze(traj_ii_jj_cc(smid,3,ii,:)),...
                'color',[1 0 0],'linewidth',2); hold on;
        scatter3(squeeze(traj_ii_jj_cc(smid,1,ii,:)),...
              squeeze(traj_ii_jj_cc(smid,2,ii,:)),...
              squeeze(traj_ii_jj_cc(smid,3,ii,:)),...
                40,repmat([1 0 0],nCurvBin,1).*repmat(binC',1,3),'filled');

        % plot(squeeze(traj_ii_jj_cc(smid,1,ii,:)),...
        %       squeeze(traj_ii_jj_cc(smid,3,ii,:)),...
        %         'color',[1 0 0],'linewidth',2); hold on;
        % scatter(squeeze(traj_ii_jj_cc(smid,1,ii,:)),...
        %       squeeze(traj_ii_jj_cc(smid,3,ii,:)),...
        %         40,repmat([1 0 0],nCurvBin,1).*(0.3+0.7*repmat(binC',1,3)),'filled');
    end
    
    plot3(squeeze(traj_ii_jj_cc(fixoff,1,ii,:)),...
          squeeze(traj_ii_jj_cc(fixoff,2,ii,:)),...
          squeeze(traj_ii_jj_cc(fixoff,3,ii,:)),...
            'color',cols(ii,:),'linewidth',2); hold on;
    scatter3(squeeze(traj_ii_jj_cc(fixoff,1,ii,:)),...
          squeeze(traj_ii_jj_cc(fixoff,2,ii,:)),...
          squeeze(traj_ii_jj_cc(fixoff,3,ii,:)),...
            40,repmat(cols(ii,:),nCurvBin,1).*repmat(binC',1,3),'filled');

    % plot(squeeze(traj_ii_jj_cc(fixoff,1,ii,:)),...
    %       squeeze(traj_ii_jj_cc(fixoff,3,ii,:)),...
    %         'color',cols(ii,:),'linewidth',2); hold on;
    % scatter(squeeze(traj_ii_jj_cc(fixoff,1,ii,:)),...
    %       squeeze(traj_ii_jj_cc(fixoff,3,ii,:)),...
    %         40,repmat(cols(ii,:),nCurvBin,1).*(0.5+0.5*repmat(binC',1,3)),'filled');
    
end

fix3dPlot(gca,[],[],[],'curvature dim','arc condition dim','sacc dim')
set(gca,'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{}); grid(gca,'off')

%%
function [xx_u,yy_m] = plotTrendPatchLine(h,xVar,yVar,col,xLims,nBin)
    % [xx_u,~,grp] = unique(xVar);
    [~,binE,grp] = histcounts(xVar,linspace(xLims(1),xLims(2),nBin+1));
    grp(grp<1) = 1; grp(grp>nBin) = nBin;
    xx_u = (binE+circshift(binE,-1))/2; xx_u = xx_u(1:end-1)';
    yy_m = groupsummary(yVar,grp,'mean'); % yy_m = smooth(yy_m);
    yy_s = groupsummary(yVar,grp,'std')./sqrt(groupsummary(yVar,grp,'nnz'));

    % fix for when there are empty bins
    if length(yy_m)~=length(xx_u)
        yy_m(unique(grp)) = yy_m; yy_m(~ismember(1:nBin,unique(grp))) = nan; yy_m = smooth(yy_m);
        yy_s(unique(grp)) = yy_s; yy_s(~ismember(1:nBin,unique(grp))) = nan; yy_s = smooth(yy_s);
    end

    hold(h,'on')
%     patch([xx_u;flipud(xx_u)],[yy_m-yy_s/2; flipud(yy_m+yy_s/2)],col,'edgecolor','none','facealpha',0.5,'parent',h)

%     beta = regress(yy_m(:),[ones(length(xx_u),1) xx_u(:)]);
%     yFit = beta(1) + beta(2)*xLims;
%     line(xLims,yFit,'linestyle','--','color',col,'linewidth',2)

    plot(h,xx_u,yy_m,'-','linewidth',2,'color',col); % 'MarkerFaceColor','w','markersize',10
end