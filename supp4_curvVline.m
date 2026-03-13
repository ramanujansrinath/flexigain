%%
clc; clear; close all;
addpath('dep')
plotSets = {4; [12 13 14]; [5 25]; [13 93 113]};
load('data/supp4_pca_data.mat')
figure('color','w','position',[119,117,637,694]);

%% select a session
sessId = 2;
sets = plotSets{sessId};
resp_v4 = supp4_pca_data(sessId).resp_v4;
resp_v1 = supp4_pca_data(sessId).resp_v1;

red_v4 = pca(resp_v4','NumComponents',3);
red_v1 = pca(resp_v1','NumComponents',3);

for ii=1:length(sets)
    idx = [supp4_pca_data(sessId).set] == sets(ii);
    cc = supp4_pca_data(sessId).curv(idx);
    ss = supp4_pca_data(sessId).sel(idx);

    red_v1_set = red_v1(idx,:);
    red_v4_set = red_v4(idx,:);

    % pca or tuning v4
    hold(subplot(221),'on'); 
    [perf(ii).cc1,perf(ii).ss1] = get_spe_dec(red_v4_set,cc,ss,1);
    [perf(ii).cc2,perf(ii).ss2] = get_spe_dec(red_v4_set,cc,ss,2);
    [perf(ii).cc3,perf(ii).ss3] = get_spe_dec(red_v4_set,cc,ss,3);
end



%% functions
function [dec_curv,dec_sel] = get_spe_dec(resp_red,cc,ss,ord)
    nFold = 10;
    dec_curv = nan(1,nFold);
    dec_sel = nan(1,nFold);
    for ff=1:nFold
        trainIdx = 1:length(cc);
        testIdx = sort(randperm(length(cc),round(length(cc)/2)));
        trainIdx(ismember(trainIdx,testIdx)) = [];
        
        resp_train = resp_red(trainIdx,:);
        cc_train = cc(trainIdx);
    
        resp_test = resp_red(testIdx,:);
        cc_test = cc(testIdx);
        ss_test = ss(testIdx);
        
        % fit data parametrically, 3rd order polynomial
        xcoefs=polyfit(cc_train,resp_train(:,1),ord); x=polyval(xcoefs,cc_train);
        ycoefs=polyfit(cc_train,resp_train(:,2),ord); y=polyval(ycoefs,cc_train);
        zcoefs=polyfit(cc_train,resp_train(:,3),ord); z=polyval(zcoefs,cc_train);
        curv_fit = [x; y; z]';
    
        cc_pred = nan(1,size(resp_test,1));
        for ii=1:size(resp_test,1)
            cc_pred(ii) = cc_test(minind(sum((repmat(resp_test(ii,:),size(curv_fit,1),1) - curv_fit).^2,2)));
        end
    
        dec_curv(ff) = corr(cc_test',cc_pred');
        dec_sel(ff) = corr(ss_test',cc_pred');
    end
    dec_curv = mean(dec_curv);
    dec_sel = mean(dec_sel);
end

