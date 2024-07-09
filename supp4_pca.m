%%
clc; clear; close all;
plotSets = {4; [12 13 14]; [5 25]; [13 93 113]};
load('data/supp4_pca_data.mat')
figure('color','w','position',[119,117,637,694]);

%% select a session
sessId = 2;
trialMarkerAl = 0.2; % draw trial dots in low alpha
sets = plotSets{sessId};
resp_v4 = supp4_pca_data(sessId).resp_v4;
resp_v1 = supp4_pca_data(sessId).resp_v1;

red_v4 = pca(resp_v4','NumComponents',10);
red_v1 = pca(resp_v1','NumComponents',10);

cols = lines(length(sets));
legendStr_shape = cell(1,length(sets)*2);
legendStr_v4_dec = cell(1,length(sets)*2);
legendStr_v1_dec = cell(1,length(sets)*2);
clf
for ii=1:length(sets)
    idx = [supp4_pca_data(sessId).set] == sets(ii);
    cc = supp4_pca_data(sessId).curv(idx);
    ss = supp4_pca_data(sessId).sel(idx);

    red_v1_set = red_v1(idx,:);
    red_v4_set = red_v4(idx,:);

    % pca or tuning v4
    hold(subplot(221),'on'); 
    plotPCfn(subplot(221),red_v4_set,cc,cols(ii,:),'V4',trialMarkerAl);
    
    % pca v1
    plotPCfn(subplot(222),red_v1_set(:,1:3),cc,cols(ii,:),'V1',trialMarkerAl);

    % dec v4
    dec_cc_v4 = decodeFn(cc',red_v4_set);
    plotTrendPatchLine(subplot(223),cc,dec_cc_v4',cols(ii,:));

    % dec v1
    dec_cc_v1 = decodeFn(cc',red_v1_set);
    plotTrendPatchLine(subplot(224),cc,dec_cc_v1',cols(ii,:));
    
    % shape ids
    legendStr_shape{2*ii-1} = '';
    legendStr_shape{ii*2} = ['shape ' num2str(sets(ii))];

    % dec acc v4
    legendStr_v4_dec{2*ii-1} = '';
    legendStr_v4_dec{ii*2} = ['acc = ' num2str(round(corr(dec_cc_v4,cc'),2))];

    % dec acc v1
    legendStr_v1_dec{2*ii-1} = '';
    legendStr_v1_dec{ii*2} = ['acc = ' num2str(round(corr(dec_cc_v1,cc'),2))];
end

fixPlot(subplot(223),[-0.1 1.1],[-0.1 1.1],'curvature','decoded curvature',0:0.25:1,0:0.25:1,'v4 decoding',legendStr_v4_dec)
legend('location','southeast')
fixPlot(subplot(224),[-0.1 1.1],[-0.1 1.1],'curvature','decoded curvature',0:0.25:1,0:0.25:1,'v1 decoding',legendStr_v1_dec)
legend('location','southeast')

ht = sgtitle(supp4_pca_data(sessId).name);
ht.Interpreter = 'none';
ht.Color = 'k'; ht.FontSize = 26; ht.FontName = 'Lato';


%% functions
function plotTrendPatchLine(h,xVar,yVar,col)
    [xx_u,~,grp] = unique(xVar);
    yy_m = groupsummary(yVar',grp,'mean');
    yy_s = groupsummary(yVar',grp,'std')./sqrt(groupsummary(yVar',grp,'nnz'));

    hold(h,'on')
    patch([xx_u fliplr(xx_u)]',[yy_m-yy_s/2; flipud(yy_m+yy_s/2)],col,'edgecolor','none','facealpha',0.5,'parent',h)
    plot(h,xx_u,yy_m,'linewidth',2,'color',col);
end

function plotPCfn(h,redResp,param,colmap,titleStr,trialMarkerAl)
    [uParam,~,paramG] = unique(param);
    cols = [linspace(0.3,colmap(1),max(paramG))' linspace(0.3,colmap(2),max(paramG))' linspace(0.3,colmap(3),max(paramG))'];
    colsF = cols(paramG,:);

    mResp = groupsummary(redResp,paramG,'mean');
    respForFit = mResp;
    paramForFit = uParam;
    if trialMarkerAl > 0 % if you don't do this, then the axis gets scaled/zoomed out for no reason
        scatter3(h,redResp(:,1),redResp(:,2),redResp(:,3),50,colsF,'filled','MarkerFaceAlpha',trialMarkerAl); 
        respForFit = redResp;
        paramForFit = param;
    end
    hold(h,'on'); axis(h,'equal','square'); 
    
    scatter3(h,mResp(:,1),mResp(:,2),mResp(:,3),80,cols,'filled'); 

    % fit data parametrically, 3rd order polynomial
    t = paramForFit;
    polyOrder = [2 2 2];
    xcoefs=polyfit(t,respForFit(:,1),polyOrder(1)); x=polyval(xcoefs,t);
    ycoefs=polyfit(t,respForFit(:,2),polyOrder(2)); y=polyval(ycoefs,t);
    zcoefs=polyfit(t,respForFit(:,3),polyOrder(3)); z=polyval(zcoefs,t);
    plot3(h,x,y,z,'color',colmap,'linestyle','-','linewidth',2); 

    fix3dPlot(h,[],[],[],'PC1','PC2','PC3',[],[],[],titleStr)
    set(h,'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{})
    grid(h,'off')
end

