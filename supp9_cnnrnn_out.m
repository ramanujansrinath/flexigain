clc; clear; close all
load('data/supp9_cnnrnn_out_model1.mat')

%%
figure('color','w','pos',[476,521,891,345])
params = cell2mat(cellfun(@(x) [x.s double(x.b) double(x.a) double(x.o)]',trial_params,'UniformOutput',false))';
arc_conds = unique(params(:,2:3),'rows');
cols = [copper(5);winter(3)]; cols(1:2,:) = []; % cols = reshape(cols',3,6)';
clf;
for ii=1:size(arc_conds,1)
    idx = params(:,2)==arc_conds(ii,1) & params(:,3)==arc_conds(ii,2);
    curv = params(idx,1);
    correct_sacc = params(idx,4);

    % last time point
    actual_sacc = outputs(:,end);
    actual_sacc = actual_sacc(idx);

    % mean after fix off
    % actual_sacc = outputs(:,trial_params{1}.fix_offset/10 : end);
    % actual_sacc = actual_sacc(params(:,2)==arc_conds(ii,1) & params(:,3)==arc_conds(ii,2),:);
    % actual_sacc = mean(actual_sacc,2);

    subplot(132); hold on;
    scatter(curv,correct_sacc,12,repmat(cols(ii,:),length(curv),1).*(0.5+0.5*repmat(curv,1,3)),'filled');

    subplot(133); hold on;
    plotTrendPatchLine(gca,curv,actual_sacc,cols(ii,:),[0 1],11);
end
fixPlot(subplot(132),[-0.1 1.1],[-90 90],'curvature','saccade target',0:0.25:1,-90:45:90,'inference to expected behaviour')
fixPlot(subplot(133),[-0.1 1.1],[-90 90],'curvature','actual saccade',0:0.25:1,-90:45:90,'inference to actual behaviour')

%%
% relevant time points
son = trial_params{1}.s_onset/10;
abon = trial_params{1}.ab_onset/10;
smid = (son+abon)/2;
fixoff = trial_params{1}.fix_offset/10;

nCurvBin = 8;
curv = cellfun(@(x) x.s,trial_params);
targ = cellfun(@(x) x.o,trial_params);
binE = linspace(0,1,nCurvBin+1);
binC = (binE+circshift(binE,-1))/2; binC = binC(1:end-1);

% qr
beta_curv = regress(curv',squeeze(state_vars(:,fixoff,:)));
beta_sacc = regress(outputs(:,fixoff),squeeze(state_vars(:,fixoff,:)));
[beta_orth,rr] = qr([beta_curv beta_sacc],'econ');
stateAtGo = squeeze(state_vars(:,fixoff,:))*beta_orth;
% stateAtStim = squeeze(state_vars(:,fixoff,:))*beta_orth;

subplot(131); hold on;
for ii=1:size(arc_conds,1)
    idx = params(:,2)==arc_conds(ii,1) & params(:,3)==arc_conds(ii,2);
    curv_ii = curv(idx);
    state_go_ii = stateAtGo(idx,:);
    % state_stim_ii = stateAtStim(idx,:);

    % plot after binning by curvatures
    [~,~,idx] = histcounts(curv_ii,binE);
    state_go_jj = nan(nCurvBin,2);
    state_stim_jj = nan(nCurvBin,2);
    for jj=1:length(binC)
        state_go_jj(jj,:) = squeeze(mean(state_go_ii(idx==jj,:)));
        % state_stim_jj(jj,:) = squeeze(mean(state_stim_ii(idx==jj,:)));
    end

    % if ii==1
    %     plot(state_stim_jj(:,1),state_stim_jj(:,2),...
    %             'color',[1 0 0],'linewidth',2); hold on;
    %     scatter(state_stim_jj(:,1),state_stim_jj(:,2),...
    %             40,repmat([1 0 0],nCurvBin,1).*(0.3+0.7*repmat(binC',1,3)),'filled');
    % end

    plot(state_go_jj(:,1),state_go_jj(:,2),...
            'color',cols(ii,:),'linewidth',2); hold on;
    scatter(state_go_jj(:,1),state_go_jj(:,2),...
            40,repmat(cols(ii,:),nCurvBin,1).*(0.5+0.5*repmat(binC',1,3)),'filled');
    
end

fixPlot(gca,[],[],'curvature dim','sacc dim')
set(gca,'XTickLabel',{},'YTickLabel',{}); grid(gca,'off')

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

    beta = regress(yy_m(:),[ones(length(xx_u),1) xx_u(:)]);
    yFit = beta(1) + beta(2)*xLims;
    line(xLims,yFit,'linestyle','--','color',col,'linewidth',2)

%     plot(h,xx_u,yy_m,'.','linewidth',2,'color',col,'MarkerFaceColor','w','markersize',7); % 'MarkerFaceColor','w','markersize',10
    scatter(h,xx_u,yy_m,20,repmat(col,length(xx_u),1).*(0.5+0.5*repmat(xx_u,1,3)),'filled'); % 'MarkerFaceColor','w','markersize',10
end