%%
clc; clear; close all;
load('data/figure2_cVals.mat')

% this code also generates supp6
figure; plotGenDec(cVals_v1); sgtitle('V1')
figure; plotGenDec(cVals_v4); sgtitle('V4')

figure; plotGenDec_marg(cVals_v1); sgtitle('V1')
figure; plotGenDec_marg(cVals_v4); sgtitle('V4')

%%
function plotGenDec(cVals)
    clf; set(gcf,'pos',[476,237,889,629],'color','w'); 
    subplot(131); hold on;
    line([-1 1],[-1 1],'color','k','linestyle','--','linewidth',2);
    plot(cVals(:,4),cVals(:,6),'.','markersize',11,'color',[0.5 0.5 0.5]);
    plot(mean(cVals(:,4)),mean(cVals(:,6)),'b+','markersize',20,'linewidth',4);
    fixPlot(gca,[-0.1 1],[-0.1 1],'shape specific','shape agnostic',0:0.25:1,0:0.25:1,'corr(decoded curv,sel)')
    
    subplot(132); hold on;
    plot(cVals(:,5),cVals(:,7),'.','markersize',11,'color',[0.5 0.5 0.5]);
    beta = regress(cVals(:,7),[ones(size(cVals,1),1) cVals(:,5)]);
    line([0 1],[beta(1) sum(beta)],'color','r','linestyle','--','linewidth',2);
    fixPlot(gca,[-0.1 1],[-0.04 0.4],'shape agnostic','shape error',0:0.25:1,0:0.2:1,{'agn corr(decoded curv,curv)' 'vs shape error'})
    
    exIds = unique(cVals(:,1));
    subplot(133);  hold on;
    line([-1 1],[-1 1],'color','k','linestyle','--','linewidth',2);
    pts_ii = [];
    for ii=1:length(exIds)
        dec = cVals(cVals(:,1) == exIds(ii),5); decPairs = nchoosek(dec,2);
        beh = cVals(cVals(:,1) == exIds(ii),7); behPairs = nchoosek(beh,2);
        pts_jj = nan(size(behPairs,1),2);
        for jj=1:size(behPairs,1)
            if decPairs(jj,1)>decPairs(jj,2)
                pts_jj(jj,:) = fliplr(behPairs(jj,:));
            else
                pts_jj(jj,:) = behPairs(jj,:);
            end
        end
        pts_ii = [pts_ii; pts_jj]; %#ok<AGROW>
    end
    plot(pts_ii(:,1),pts_ii(:,2),'.','markersize',11,'color',[0.5 0.5 0.5]);
    plot(mean(pts_ii(:,1)),mean(pts_ii(:,2)),'b+','markersize',20,'linewidth',4);
    fixPlot(gca,[-0.04 0.4],[-0.04 0.4],'avg beh error for worse dec','avg beh error for better dec',0:0.2:1,0:0.2:1,{'avg beh error for better' 'and worse decoded shapes'})
end

function plotGenDec_marg(cVals)
    x = cVals(:,4);
    y = cVals(:,6);
    t = cVals(:,8);

    cols = {[0 0 0] [0.5 0.25 0.5] [0.4 0.6 0.6]};
    
    set(gcf,'color','w','Position',[-1597 215 822 764])
    subplot(221); hold on;
    histogram(x(t==1),linspace(-0.1,1,20),'DisplayStyle','stairs','EdgeColor',cols{1},'LineWidth',2);
    histogram(x(t==2),linspace(-0.1,1,20),'DisplayStyle','stairs','EdgeColor',cols{2},'LineWidth',2);
    histogram(x(t==3),linspace(-0.1,1,20),'DisplayStyle','stairs','EdgeColor',cols{3},'LineWidth',2);
    line([mean(x(t==1)) mean(x(t==1))],[50 120],'linewidth',2,'color',cols{1},'linestyle','--')
    line([mean(x(t==2)) mean(x(t==2))],[50 120],'linewidth',2,'color',cols{2},'linestyle','--')
    line([mean(x(t==3)) mean(x(t==3))],[50 120],'linewidth',2,'color',cols{3},'linestyle','--')
    line([0 0],[0 120],'linewidth',2,'color','k','linestyle','--')
    fixPlot(gca,[-0.1 1],[0 40],'','',-1:0.5:1,0:30:200); set(gca,'XTickLabel',{})
    
    subplot(224); hold on;
    histogram(y(t==1),linspace(-0.1,1,20),'DisplayStyle','stairs','EdgeColor',cols{1},'LineWidth',2);
    histogram(y(t==2),linspace(-0.1,1,20),'DisplayStyle','stairs','EdgeColor',cols{2},'LineWidth',2);
    histogram(y(t==3),linspace(-0.1,1,20),'DisplayStyle','stairs','EdgeColor',cols{3},'LineWidth',2);
    line([mean(y(t==1)) mean(y(t==1))],[50 120],'linewidth',2,'color',cols{1},'linestyle','--')
    line([mean(y(t==2)) mean(y(t==2))],[50 120],'linewidth',2,'color',cols{2},'linestyle','--')
    line([mean(y(t==3)) mean(y(t==3))],[50 120],'linewidth',2,'color',cols{3},'linestyle','--')
    line([0 0],[0 120],'linewidth',2,'color','k','linestyle','--')
    view([90 -90])
    fixPlot(gca,[-0.1 1],[0 40],'','',-1:0.5:1,0:30:200); set(gca,'XTickLabel',{})
    
    subplot(223); hold on;
    plot(x(t==1),y(t==1),'o','MarkerSize',7,'MarkerFaceColor','w','LineWidth',2,'color',cols{1})
    plot(x(t==2),y(t==2),'o','MarkerSize',7,'MarkerFaceColor','w','LineWidth',2,'color',cols{2})
    plot(x(t==3),y(t==3),'o','MarkerSize',7,'MarkerFaceColor','w','LineWidth',2,'color',cols{3})
    line([-1 1],[-1 1],'color','k','linestyle','--','linewidth',2);
    fixPlot(gca,[-0.1 1],[-0.1 1],'shape specific','shape general',-1:0.25:1,-1:0.25:1)
end