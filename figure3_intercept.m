%% get all data ids
clc; clear; close all;
load('data/figure3_intercept.mat','v1','v4')
intercept_mid = v4.intercept_mid;
intercept_len = v4.intercept_len;

%% plot scatter
intercept_mid(~ismember([intercept_mid.mid_val],-20:20:20)) = [];
intercept_len([intercept_len.sess_mark]==1) = [];

figure('color','w','pos',[476,305,971,561])
subplot(241); hold on;
a = intercept_mid;
% a = stats_mid([stats_mid.len_val]==100);
% b = stats_mid([stats_mid.len_val]==140);
plot(-2+4*rand(size(a))+[a.mid_val],[a.beta_stim_inter],'k.','markersize',6,'linewidth',1.5,'markerfacecolor','w'); hold on;
% plot(-2+4*rand(size(b))+[b.mid_val],[b.beta_stim_inter],'r.','markersize',6,'linewidth',1.5,'markerfacecolor','w');
plot(unique([a.mid_val]),groupsummary([a.beta_stim_inter]',[a.mid_val]','median'),'ko-','markersize',10,'linewidth',1.5,'markerfacecolor','w');
% plot(unique([b.mid_val]),groupsummary([b.beta_stim_inter]',[b.mid_val]','median'),'ro-','markersize',10,'linewidth',1.5,'markerfacecolor','w');
fixPlot(gca,[-25 25],[-35 35],'arc mid value','intercept of saccade decoder',-20:20:40,-20:20:40,'stim time')

% histograms of the stuff in the previous plot
subplot(242); hold on;
cla; hold on;
histogram([a([a.mid_val]==-20).beta_stim_inter],linspace(-35,35,13),'DisplayStyle','stairs','linewidth',2,'edgecolor',[62 83 159]/255)
histogram([a([a.mid_val]==0).beta_stim_inter],linspace(-35,35,13),'DisplayStyle','stairs','linewidth',2,'edgecolor',[53 126 187]/255)
histogram([a([a.mid_val]==20).beta_stim_inter],linspace(-35,35,13),'DisplayStyle','stairs','linewidth',2,'edgecolor',[109 159 96]/255)
fixPlot(gca,[-35 35],[0 35],'intercept of saccade decoder','count',-30:15:30,0:10:30,'stim time')

subplot(243); hold on;
plot(-2+4*rand(size(a))+[a.mid_val],[a.beta_arc_inter],'k.','markersize',6,'linewidth',1.5,'markerfacecolor','w'); hold on;
% plot(-2+4*rand(size(b))+[b.mid_val],[b.beta_arc_inter],'r.','markersize',6,'linewidth',1.5,'markerfacecolor','w');
plot(unique([a.mid_val]),groupsummary([a.beta_arc_inter]',[a.mid_val]','median'),'ko-','markersize',10,'linewidth',1.5,'markerfacecolor','w');
% plot(unique([b.mid_val]),groupsummary([b.beta_arc_inter]',[b.mid_val]','median'),'ro-','markersize',10,'linewidth',1.5,'markerfacecolor','w');
fixPlot(gca,[-25 25],[-35 35],'arc mid value','intercept of saccade decoder',-20:20:40,-20:20:40,'arc time')

% histograms of the stuff in the previous plot
subplot(244); hold on;
histogram([a([a.mid_val]==-20).beta_arc_inter],linspace(-35,35,13),'DisplayStyle','stairs','linewidth',2,'edgecolor',[62 83 159]/255)
histogram([a([a.mid_val]==0).beta_arc_inter],linspace(-35,35,13),'DisplayStyle','stairs','linewidth',2,'edgecolor',[53 126 187]/255)
histogram([a([a.mid_val]==20).beta_arc_inter],linspace(-35,35,13),'DisplayStyle','stairs','linewidth',2,'edgecolor',[109 159 96]/255)
fixPlot(gca,[-35 35],[0 35],'intercept of saccade decoder','count',-30:15:30,0:10:30,'arc time')

subplot(245); cla; hold on;
a = intercept_len([intercept_len.mid_val]==0);
% b = stats_len([stats_len.mid_val]==-20);
c = intercept_len([intercept_len.mid_val]==20);
plot(-2+4*rand(size(a))+[a.len_val],[a.beta_stim_inter],'k.','markersize',6,'linewidth',1.5,'markerfacecolor','w'); hold on;
% plot(-2+4*rand(size(b))+[b.len_val],[b.beta_stim_inter],'r.','markersize',6,'linewidth',1.5,'markerfacecolor','w');
plot(-2+4*rand(size(c))+[c.len_val],[c.beta_stim_inter],'g.','markersize',6,'linewidth',1.5,'markerfacecolor','w');
plot(unique([a.len_val]),groupsummary([a.beta_stim_inter]',[a.len_val]','median'),'ko-','markersize',10,'linewidth',1.5,'markerfacecolor','w');
% plot(unique([b.len_val]),groupsummary([b.beta_stim_inter]',[b.len_val]','median'),'ro-','markersize',10,'linewidth',1.5,'markerfacecolor','w');
plot(unique([c.len_val]),groupsummary([c.beta_stim_inter]',[c.len_val]','median'),'go-','markersize',10,'linewidth',1.5,'markerfacecolor','w');
fixPlot(gca,[90 150],[-35 35],'arc len value','slope of saccade decoder',[100 140],-20:20:100,'stim time')

subplot(246); cla; hold on;
histogram([a([a.len_val]==100).beta_stim_inter],linspace(-35,35,13),'DisplayStyle','stairs','linewidth',2,'edgecolor',[115 192 235]/255)
histogram([a([a.len_val]==140).beta_stim_inter],linspace(-35,35,13),'DisplayStyle','stairs','linewidth',2,'edgecolor',[53 126 187]/255)
histogram([c([c.len_val]==100).beta_stim_inter],linspace(-35,35,13),'DisplayStyle','stairs','linewidth',2,'edgecolor',[187 217 174]/255)
histogram([c([c.len_val]==140).beta_stim_inter],linspace(-35,35,13),'DisplayStyle','stairs','linewidth',2,'edgecolor',[131 191 116]/255)
fixPlot(gca,[-35 35],[0 12],'intercept of saccade decoder','count',-30:15:30,0:10:30,'stim time')


subplot(247); cla; hold on;
plot(-2+4*rand(size(a))+[a.len_val],[a.beta_arc_inter],'k.','markersize',6,'linewidth',1.5,'markerfacecolor','w'); hold on;
% plot(-2+4*rand(size(b))+[b.len_val],[b.beta_arc_inter],'r.','markersize',6,'linewidth',1.5,'markerfacecolor','w');
plot(-2+4*rand(size(c))+[c.len_val],[c.beta_arc_inter],'g.','markersize',6,'linewidth',1.5,'markerfacecolor','w');
plot(unique([a.len_val]),groupsummary([a.beta_arc_inter]',[a.len_val]','median'),'ko-','markersize',10,'linewidth',1.5,'markerfacecolor','w');
% plot(unique([b.len_val]),groupsummary([b.beta_arc_inter]',[b.len_val]','median'),'ro-','markersize',10,'linewidth',1.5,'markerfacecolor','w');
plot(unique([c.len_val]),groupsummary([c.beta_arc_inter]',[c.len_val]','median'),'go-','markersize',10,'linewidth',1.5,'markerfacecolor','w');
fixPlot(gca,[90 150],[-35 35],'arc len value','slope of saccade decoder',[100 140],-20:20:100,'arc time')

subplot(248); cla; hold on;
histogram([a([a.len_val]==100).beta_arc_inter],linspace(-35,35,13),'DisplayStyle','stairs','linewidth',2,'edgecolor',[115 192 235]/255)
histogram([a([a.len_val]==140).beta_arc_inter],linspace(-35,35,13),'DisplayStyle','stairs','linewidth',2,'edgecolor',[53 126 187]/255)
histogram([c([c.len_val]==100).beta_arc_inter],linspace(-35,35,13),'DisplayStyle','stairs','linewidth',2,'edgecolor',[187 217 174]/255)
histogram([c([c.len_val]==140).beta_arc_inter],linspace(-35,35,13),'DisplayStyle','stairs','linewidth',2,'edgecolor',[131 191 116]/255)
fixPlot(gca,[-35 35],[0 12],'intercept of saccade decoder','count',-30:15:30,0:10:30,'arc time')


%% decoding accuracy histograms
figure('color','w','pos',[742,37,596,326])
subplot(121); hold on;
[~,uSess] = unique([[intercept_mid.sessId]' [intercept_mid.set]'],'rows');
histogram([intercept_mid(uSess).curv_dec_set_stim],linspace(0,1,15),'DisplayStyle','stairs','linewidth',1.5);
histogram([intercept_mid(uSess).curv_dec_set_arc],linspace(0,1,14),'DisplayStyle','stairs','linewidth',1.5)
histogram([intercept_mid(uSess).sacc_dec_set_arc],linspace(0,1,16),'DisplayStyle','stairs','linewidth',1.5)
fixPlot(gca,[-0.1 1.1],[0 30],'decoding acc','count',0:0.25:1,0:10:100,'mid sessions',{'stim time curv dec' 'arc time curv dec' 'arc time sacc dec'})
legend('Location','northwest')

subplot(122); hold on;
[~,uSess] = unique([[intercept_len.sessId]' [intercept_len.set]'],'rows');
histogram([intercept_len(uSess).curv_dec_set_stim],linspace(0,1,9),'DisplayStyle','stairs','linewidth',1.5);
histogram([intercept_len(uSess).curv_dec_set_arc],linspace(0,1,8),'DisplayStyle','stairs','linewidth',1.5)
histogram([intercept_len(uSess).sacc_dec_set_arc],linspace(0,1,10),'DisplayStyle','stairs','linewidth',1.5)
fixPlot(gca,[-0.1 1.1],[0 15],'decoding acc','count',0:0.25:1,0:5:100,'len sessions')