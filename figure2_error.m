%%
clc; clear; close all;
load('data/figure2_error.mat')

figure; plotGenError(rec_v1); sgtitle('V1')
figure; plotGenError(rec_v4); sgtitle('V4')

%%
function plotGenError(rec)
    x = [rec.corr_decErr_behErr_spe];
    y = [rec.corr_decErr_behErr_gen];
    t = [rec.exptType];

    cols = {[0 0 0] [0.5 0.25 0.5] [0.4 0.6 0.6]};
    
    set(gcf,'color','w','Position',[-1597 215 822 764])
    subplot(221); hold on;
    histogram(x(t==1),linspace(-0.2,1,20),'DisplayStyle','stairs','EdgeColor',cols{1},'LineWidth',2);
    histogram(x(t==2),linspace(-0.2,1,20),'DisplayStyle','stairs','EdgeColor',cols{2},'LineWidth',2);
    histogram(x(t==3),linspace(-0.2,1,20),'DisplayStyle','stairs','EdgeColor',cols{3},'LineWidth',2);
    line([mean(x(t==1)) mean(x(t==1))],[50 120],'linewidth',2,'color',cols{1},'linestyle','--')
    line([mean(x(t==2)) mean(x(t==2))],[50 120],'linewidth',2,'color',cols{2},'linestyle','--')
    line([mean(x(t==3)) mean(x(t==3))],[50 120],'linewidth',2,'color',cols{3},'linestyle','--')
    line([0 0],[0 120],'linewidth',2,'color','k','linestyle','--')
    fixPlot(gca,[-0.2 1],[0 60],'','',-1:0.5:1,0:30:200); set(gca,'XTickLabel',{})
    
    subplot(224); hold on;
    histogram(y(t==1),linspace(-0.2,1,20),'DisplayStyle','stairs','EdgeColor',cols{1},'LineWidth',2);
    histogram(y(t==2),linspace(-0.2,1,20),'DisplayStyle','stairs','EdgeColor',cols{2},'LineWidth',2);
    histogram(y(t==3),linspace(-0.2,1,20),'DisplayStyle','stairs','EdgeColor',cols{3},'LineWidth',2);
    line([mean(y(t==1)) mean(y(t==1))],[50 120],'linewidth',2,'color',cols{1},'linestyle','--')
    line([mean(y(t==2)) mean(y(t==2))],[50 120],'linewidth',2,'color',cols{2},'linestyle','--')
    line([mean(y(t==3)) mean(y(t==3))],[50 120],'linewidth',2,'color',cols{3},'linestyle','--')
    line([0 0],[0 120],'linewidth',2,'color','k','linestyle','--')
    view([90 -90])
    fixPlot(gca,[-0.2 1],[0 60],'','',-1:0.5:1,0:30:200); set(gca,'XTickLabel',{})
    
    subplot(223); hold on;
    plot(x(t==1),y(t==1),'.','MarkerSize',12,'MarkerFaceColor','w','LineWidth',2,'color',cols{1})
    plot(x(t==2),y(t==2),'.','MarkerSize',12,'MarkerFaceColor','w','LineWidth',2,'color',cols{2})
    plot(x(t==3),y(t==3),'.','MarkerSize',12,'MarkerFaceColor','w','LineWidth',2,'color',cols{3})
    line([-1 1],[-1 1],'color','k','linestyle','--','linewidth',2);
    fixPlot(gca,[-0.2 1],[-0.2 1],'shape specific','shape general',-1:0.5:1,-1:0.5:1)
    
    ht = sgtitle('correlation between the behavioral error and decoding error');
    ht.FontName = 'lato'; ht.FontSize = 20; ht.FontWeight = 'bold';
end