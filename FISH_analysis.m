clear
close all

load("FISH_kr_hb.mat");
load("FISH_eve_hb.mat");
load('genotype_ID.mat')

reduce = 1;

%% smooth kr and kni signals
for i = 1:60
    krS(i,:) = smooth(kr(i,:));
    hbS(i,:) = smooth(hb(i,:));
end

for i = 1:60
    eveS(i,:) = eve(i,:);
end

for i = 1:length(wt_ind)
    wt_kr(i,:) = krS(wt_ind(i),:);
    wt_hb(i,:) = hbS(wt_ind(i),:);
    wt_eve(i,:) = eveS(wt_ind(i),:);
end
for i = 1:length(het_ind)
    het_kr(i,:) = krS(het_ind(i),:);
    het_hb(i,:) = hbS(het_ind(i),:);
    het_eve(i,:) = eveS(het_ind(i),:);
end
%% subtract minimun signals
if reduce == 1
    for i = 1:60
        krS(i,:) = krS(i,:)-min(krS(i,10:90));
        hbS(i,:) = hbS(i,:)-min(hbS(i,10:90));
        eveS(i,:) = eveS(i,:)-min(eveS(i,30:80));
    end
end

%% set wt and het arrays
for i = 1:length(wt_ind)
    wt_kr(i,:) = krS(wt_ind(i),:);
    wt_hb(i,:) = hbS(wt_ind(i),:);
    wt_eve(i,:) = eveS(wt_ind(i),:);
end
for i = 1:length(het_ind)
    het_kr(i,:) = krS(het_ind(i),:);
    het_hb(i,:) = hbS(het_ind(i),:);
    het_eve(i,:) = eveS(het_ind(i),:);
end

%% average trajectory LEM
[tmp wt_lem_ind] = intersect(wt_ind, wt_lem);
[tmp het_lem_ind] = intersect(het_ind, het_lem);

for i = 1:100
    krAveLEM_wt(i) = mean(wt_kr(wt_lem_ind,i)); %ave kr trajectory
    krAveLEM_het(i) = mean(het_kr(het_lem_ind,i)); %ave kr trajectory

    hbAveLEM_wt(i) = mean(wt_hb(wt_lem_ind,i)); %ave hb trajectory
    hbAveLEM_het(i) = mean(het_hb(het_lem_ind,i)); %ave hb trajectory

    eveAveLEM_wt(i) = mean(wt_eve(wt_lem_ind,i)); %ave eve trajectory
    eveAveLEM_het(i) = mean(het_eve(het_lem_ind,i)); %ave eve trajectory

    krAveLEMSE_wt(i) = std(wt_kr(wt_lem_ind,i))/sqrt(length(wt_lem_ind)); %ave kr trajectory
    krAveLEMSE_het(i) = std(het_kr(het_lem_ind,i))/sqrt(length(het_lem_ind)); %ave kr trajectory

    hbAveLEMSE_wt(i) = std(wt_hb(wt_lem_ind,i))/sqrt(length(wt_lem_ind)); %ave hb trajectory
    hbAveLEMSE_het(i) = std(het_hb(het_lem_ind,i))/sqrt(length(het_lem_ind)); %ave hb trajectory

    eveAveLEMSE_wt(i) = std(wt_eve(wt_lem_ind,i))/sqrt(length(wt_lem_ind)); %ave eve trajectory
    eveAveLEMSE_het(i) = std(het_eve(het_lem_ind,i))/sqrt(length(het_lem_ind)); %ave eve trajectory

end





%% make plots
figure(201); 
plot(11:90, krAveLEM_wt(11:90)/10000,'-k','LineWidth',3); hold on
plot(11:90, krAveLEM_het(11:90)/10000,'-r','LineWidth',3); 
shadedErrorBar2(11:90, krAveLEM_wt(11:90)/10000,krAveLEMSE_wt(11:90)/10000,'lineprops', {'-k','LineWidth',2});hold on;
shadedErrorBar2(11:90, krAveLEM_het(11:90)/10000,krAveLEMSE_het(11:90)/10000,'lineprops', {'-r','LineWidth',2});
legend('\itwt','\itKr^{1}/+');
% title('\itKr');
ylim([0 2]);
xlim([11 90]);


figure(202);
plot(11:70, hbAveLEM_wt(11:70)/1000,'-k','LineWidth',3); hold on
plot(11:70, hbAveLEM_het(11:70)/1000,'-r','LineWidth',3); 
shadedErrorBar2(11:70, hbAveLEM_wt(11:70)/1000,hbAveLEMSE_wt(11:70)/1000,'lineprops', {'-k','LineWidth',2});
shadedErrorBar2(11:70, hbAveLEM_het(11:70)/1000,hbAveLEMSE_het(11:70)/1000,'lineprops', {'-r','LineWidth',2});
legend('\itwt','\itKr^{1}/+','Location','northwest');
% title('\ithb');
ylim([0 6]);
xlim([11 70]);


%% measured mid-sagittal AP and DV lengths
TT = readtable('measured_EL_raw_image.xlsx');

lengthMid_wt = TT{:,3}';
widthMid_wt = TT{:,4}';

lengthMid_het = TT{1:22,9}';
widthMid_het = TT{1:22,10}';

[a,b,c,d] = ttest2(lengthMid_wt, lengthMid_het)
[d,e,g,g] = ttest2(widthMid_wt, widthMid_het)

figure(451);
    boxplot([lengthMid_wt lengthMid_het],...
        [ones(1,length(lengthMid_wt)), 2*ones(1,length(lengthMid_het))],...
        'Labels',{'\itwt','\itKr^{1}/+'},'symbol','');hold on
    set(gca, 'TickLabelInterpreter', 'tex');
    set(findobj(gca,'type','line'),'linew',2);
    s1 = swarmchart(ones(1,length(lengthMid_wt)),lengthMid_wt,10,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
    s1.XJitterWidth = 0.3;
    s2 = swarmchart(2*ones(1,length(lengthMid_het)),lengthMid_het,10,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
    s2.XJitterWidth = 0.3;
    % title('S2 mRNA output by nuclei (mature stripe nuclei)');
% ylabel('FISH AP length (\mum)');
    ylim([420 540]);
    hold off

    figure(452);
    boxplot([widthMid_wt widthMid_het],...
        [ones(1,length(widthMid_wt)), 2*ones(1,length(widthMid_het))],...
        'Labels',{'\itwt','\itKr^{1}/+'},'symbol','');hold on
    set(gca, 'TickLabelInterpreter', 'tex');
    set(findobj(gca,'type','line'),'linew',2);
    s1 = swarmchart(ones(1,length(widthMid_wt)),widthMid_wt,10,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
    s1.XJitterWidth = 0.3;
    s2 = swarmchart(2*ones(1,length(widthMid_het)),widthMid_het,10,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
    s2.XJitterWidth = 0.3;
    % title('S2 mRNA output by nuclei (mature stripe nuclei)');
%     ylabel('FISH DV length (\mum)');
    ylim([150 220]);
    hold off
