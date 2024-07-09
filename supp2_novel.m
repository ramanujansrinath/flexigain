clc; clear; close all;
load('data/supp2_novel.mat')

%%
cols = [0.9 0.5 0.2; 0.2 0.5 0.9];

figure('pos',[230,302,955,426],'color','w')
subplot(245)
plot(avg_error(monkey==1 & firstTrial==1 & trained==1),avg_error(monkey==1 & firstTrial==0 & trained==1),'.','color',cols(1,:)); hold on;
plot(avg_error(monkey==2 & firstTrial==1 & trained==1),avg_error(monkey==2 & firstTrial==0 & trained==1),'.','color',cols(2,:))
fixPlot(gca,[0 1],[0 1],'first trial error','11-end trial error',0:0.25:1,0:0.25:1,'trained shapes')

subplot(247)
plot(avg_error(monkey==1 & firstTrial==1 & trained==0),avg_error(monkey==1 & firstTrial==0 & trained==0),'.','color',cols(1,:)); hold on;
plot(avg_error(monkey==2 & firstTrial==1 & trained==0),avg_error(monkey==2 & firstTrial==0 & trained==0),'.','color',cols(2,:))
fixPlot(gca,[0 1],[0 1],'first trial error','11-end trial error',0:0.25:1,0:0.25:1,'novel shapes',{'monkey 1' 'monkey 2'})

subplot(241)
histogram(avg_error(monkey==1 & firstTrial==1 & trained==1),linspace(0,1,18),'DisplayStyle','stairs','EdgeColor',cols(1,:)); hold on
histogram(avg_error(monkey==2 & firstTrial==1 & trained==1),linspace(0,1,18),'DisplayStyle','stairs','EdgeColor',cols(2,:));
fixPlot(gca,[0 1],[0 200],'','')

subplot(246)
histogram(avg_error(monkey==1 & firstTrial==0 & trained==1),linspace(0,1,18),'DisplayStyle','stairs','EdgeColor',cols(1,:)); hold on
histogram(avg_error(monkey==2 & firstTrial==0 & trained==1),linspace(0,1,18),'DisplayStyle','stairs','EdgeColor',cols(2,:));
fixPlot(gca,[0 1],[0 200],'','')

subplot(243)
histogram(avg_error(monkey==1 & firstTrial==1 & trained==0),linspace(0,1,18),'DisplayStyle','stairs','EdgeColor',cols(1,:)); hold on
histogram(avg_error(monkey==2 & firstTrial==1 & trained==0),linspace(0,1,18),'DisplayStyle','stairs','EdgeColor',cols(2,:));
fixPlot(gca,[0 1],[0 200],'','')

subplot(248)
histogram(avg_error(monkey==1 & firstTrial==0 & trained==0),linspace(0,1,18),'DisplayStyle','stairs','EdgeColor',cols(1,:)); hold on
histogram(avg_error(monkey==2 & firstTrial==0 & trained==0),linspace(0,1,18),'DisplayStyle','stairs','EdgeColor',cols(2,:));
fixPlot(gca,[0 1],[0 200],'','')
