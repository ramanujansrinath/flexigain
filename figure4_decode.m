%% prelim
close all; clear; clc

% i'm going to load all the decoder performances for the curvature, color,
% and choice decoding here. i can send you the raw data against which you
% can run the code below ("function decode_all_perSession") to get the
% decoder performances.
load('data/figure4_decode.mat','r_sess')

%% boxplots
figure('pos',[1000,319,601,919],'color','w'); clf
titleStr = {{'shape decoder' 'across tasks'} {'color decoder' 'across tasks'} {'shape decoder' 'shape task'} {'color decoder' 'shape task'} {'shape decoder' 'color task'} {'color decoder' 'color task'}};
for ii=1:6
    subplot(3,2,ii); 
    notBoxPlot([squeeze(r_sess(ii,1,:)),squeeze(r_sess(ii,2,:)),squeeze(r_sess(ii,3,:))]); 
    fixPlot(gca,[0 4],[-0.2 1],'','decoding accuracy',1:3,0:0.25:1,titleStr{ii})
    set(gca,'XTickLabel',{'shape' 'color' 'choice'})
end

%% scatters
col = lines(4);
figure('pos',[1000,319,601,919],'color','w'); clf
subplot(321)
line([0 1],[0 1],'linestyle','--','color','k','linewidth',2); hold on;
plot(squeeze(r_sess(3,1,:)),squeeze(r_sess(5,1,:)),'o','markerfacecolor','w','markersize',8,'linewidth',2,'color',col(2,:));
plot(squeeze(r_sess(3,2,:)),squeeze(r_sess(5,2,:)),'o','markerfacecolor','w','markersize',8,'linewidth',2,'color',col(3,:));
plot(squeeze(r_sess(3,3,:)),squeeze(r_sess(5,3,:)),'o','markerfacecolor','w','markersize',8,'linewidth',2,'color',col(4,:));
plot(mean(squeeze(r_sess(3,1,:))),mean(squeeze(r_sess(5,1,:))),'+','markersize',30,'color',col(2,:),'LineWidth',2);
plot(mean(squeeze(r_sess(3,2,:))),mean(squeeze(r_sess(5,2,:))),'+','markersize',30,'color',col(3,:),'LineWidth',2);
plot(mean(squeeze(r_sess(3,3,:))),mean(squeeze(r_sess(5,3,:))),'+','markersize',30,'color',col(4,:),'LineWidth',2);

fixPlot(gca,[-0.2 1],[-0.2 1],'shape task trials','color task trials',0:0.25:1,0:0.25:1,'shape decoder performance')
% legend('Location','northwest')

subplot(322)
line([0 1],[0 1],'linestyle','--','color','k','linewidth',2); hold on;
plot(squeeze(r_sess(4,1,:)),squeeze(r_sess(6,1,:)),'o','markerfacecolor','w','markersize',8,'linewidth',2,'color',col(2,:));
plot(squeeze(r_sess(4,2,:)),squeeze(r_sess(6,2,:)),'o','markerfacecolor','w','markersize',8,'linewidth',2,'color',col(3,:));
plot(squeeze(r_sess(4,3,:)),squeeze(r_sess(6,3,:)),'o','markerfacecolor','w','markersize',8,'linewidth',2,'color',col(4,:));
plot(mean(squeeze(r_sess(4,1,:))),mean(squeeze(r_sess(6,1,:))),'+','markersize',30,'color',col(2,:),'LineWidth',2);
plot(mean(squeeze(r_sess(4,2,:))),mean(squeeze(r_sess(6,2,:))),'+','markersize',30,'color',col(3,:),'LineWidth',2);
plot(mean(squeeze(r_sess(4,3,:))),mean(squeeze(r_sess(6,3,:))),'+','markersize',30,'color',col(4,:),'LineWidth',2);
fixPlot(gca,[-0.2 1],[-0.2 1],'shape task trials','color task trials',0:0.25:1,0:0.25:1,'color decoder performance',{'' 'shape decoding' 'color decoding' 'choice decoding'})
legend('Location','southeast')

subplot(323); hold on; % marg dist for x axis for left scatter
histogram(squeeze(r_sess(3,1,:)),linspace(-0.2,1,15),'DisplayStyle','stairs','edgecolor',col(2,:),'linewidth',2)
histogram(squeeze(r_sess(3,2,:)),linspace(-0.2,1,15),'DisplayStyle','stairs','edgecolor',col(3,:),'linewidth',2)
histogram(squeeze(r_sess(3,3,:)),linspace(-0.2,1,15),'DisplayStyle','stairs','edgecolor',col(4,:),'linewidth',2)
fixPlot(gca,[-0.2 1],[0 20],{'shape decoder perf','shape task trials'},'count',0:0.25:1,0:10:20)

subplot(324); hold on; % marg dist for x axis for right scatter
histogram(squeeze(r_sess(4,1,:)),linspace(-0.2,1,15),'DisplayStyle','stairs','edgecolor',col(2,:),'linewidth',2)
histogram(squeeze(r_sess(4,2,:)),linspace(-0.2,1,15),'DisplayStyle','stairs','edgecolor',col(3,:),'linewidth',2)
histogram(squeeze(r_sess(4,3,:)),linspace(-0.2,1,15),'DisplayStyle','stairs','edgecolor',col(4,:),'linewidth',2)
fixPlot(gca,[-0.2 1],[0 20],{'color decoder perf','shape task trials'},'count',0:0.25:1,0:10:20)

subplot(325); hold on; % marg dist for y axis for left scatter
histogram(squeeze(r_sess(5,1,:)),linspace(-0.2,1,15),'DisplayStyle','stairs','edgecolor',col(2,:),'linewidth',2)
histogram(squeeze(r_sess(5,2,:)),linspace(-0.2,1,15),'DisplayStyle','stairs','edgecolor',col(3,:),'linewidth',2)
histogram(squeeze(r_sess(5,3,:)),linspace(-0.2,1,15),'DisplayStyle','stairs','edgecolor',col(4,:),'linewidth',2)
fixPlot(gca,[-0.2 1],[0 20],{'shape decoder perf','color task trials'},'count',0:0.25:1,0:10:20)

subplot(326); hold on; % marg dist for y axis for right scatter
histogram(squeeze(r_sess(6,1,:)),linspace(-0.2,1,15),'DisplayStyle','stairs','edgecolor',col(2,:),'linewidth',2)
histogram(squeeze(r_sess(6,2,:)),linspace(-0.2,1,15),'DisplayStyle','stairs','edgecolor',col(3,:),'linewidth',2)
histogram(squeeze(r_sess(6,3,:)),linspace(-0.2,1,15),'DisplayStyle','stairs','edgecolor',col(4,:),'linewidth',2)
fixPlot(gca,[-0.2 1],[0 20],{'color decoder perf','color task trials'},'count',0:0.25:1,0:10:20)

%%
% this is the code to find the curvature and color decoders and try to
% decode the curvature, color, and choices on them
% there are six decoders:
%    1. curvature decoder trained on all trials
%    2. color decoder trained on all trials
%    3. curvature decoder trained on curvature discrimination trials
%    4. color decoder trained on curvature discrimination trials
%    5. curvature decoder trained on color discrimination trials
%    6. color decoder trained on color discrimination trials
% since i have trained the decoder in an odd order, i reorder the
% correlation values such that they match up with the order in the box
% plots. this is a sort of silly way to do it.
function r_sess = decode_all_perSession(params,resp,resp_base) %#ok<DEFNU>
    colorDiff = (mod([params.stimRF_num],5)-mod([params.stimOpp_num],5)) == 0;
    eid = [3*ones(1,32) ones(1,64) 3*ones(1,32)];
    resp = resp(:,eid==1);
    resp_base = resp_base(:,eid==1);
    good = mean(resp)>1.1*mean(resp_base); % true(1,size(resp,2));
    resp = resp(:,good);

    resp_red = pca(resp','NumComponents',10);
    
    stimVals = [mat2vec(repmat(1:5,5,1)) mat2vec(repmat(1:5,1,5))];
    colIds = stimVals([params.stimRF_num],1);
    shpIds = stimVals([params.stimRF_num],2);
    choices = [params.selected]';
    
    % the seemingly random shuffling after r gets assigned
    % is to format it in the way that i will plot it later.
    % it's silly.

    % all trials
    r = get_dec(resp_red,shpIds,resp_red,colIds,resp_red,choices); r_sess([1 2 3]) = r;
    r = get_dec(resp_red,colIds,resp_red,shpIds,resp_red,choices); r_sess([5 4 6]) = r;
    
    % shape task
    r = get_dec(resp_red(~colorDiff,:),shpIds(~colorDiff),resp_red(~colorDiff,:),colIds(~colorDiff),resp_red(~colorDiff,:),choices(~colorDiff)); r_sess([7 8 9]) = r;
    r = get_dec(resp_red(~colorDiff,:),colIds(~colorDiff),resp_red(~colorDiff,:),shpIds(~colorDiff),resp_red(~colorDiff,:),choices(~colorDiff)); r_sess([11 10 12]) = r;

    % color task
    r = get_dec(resp_red(colorDiff,:),shpIds(colorDiff),resp_red(colorDiff,:),colIds(colorDiff),resp_red(colorDiff,:),choices(colorDiff)); r_sess([13 14 15]) = r;
    r = get_dec(resp_red(colorDiff,:),colIds(colorDiff),resp_red(colorDiff,:),shpIds(colorDiff),resp_red(colorDiff,:),choices(colorDiff)); r_sess([17 16 18]) = r;

    r_sess = reshape(r_sess,[3 6])';
end

% train on x1, test on x1, x2, and x3
function r = get_dec(x1,y1,x2,y2,x3,y3)
    nFold = 100;
    r_fold = nan(nFold,3);
    nTrainTestTrials = floor(length(y1)/2);
    for ff=1:nFold
        train_1_idx = sort(randperm(length(y1),nTrainTestTrials));
        test_1_idx = 1:length(y1); 
        test_1_idx(ismember(test_1_idx,train_1_idx)) = [];

        beta = regress(y1(train_1_idx),x1(train_1_idx,:));

        r_fold(ff,1) = corr(x1(test_1_idx,:)*beta,y1(test_1_idx));
        r_fold(ff,2) = corr(x2*beta,y2);
        r_fold(ff,3) = corr(x3*beta,y3);
    end
    mult =  repmat(2*(0.5-double(sum(r_fold<0)==nFold)),nFold,1);
    r_fold = r_fold.*mult;
    r = mean(r_fold);
end