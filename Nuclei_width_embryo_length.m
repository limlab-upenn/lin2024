clear
close all

%% load mat
load('D:\LIM LAB\kr_project_code_cleanup\dataset_whole_kr_NC14_100f.mat');
load('D:\LIM LAB\kr_project_code_cleanup\dataset_whole_hb_NC14_100f.mat');

%% load embryo image and get average nuclei diameter
load('D:\LIM LAB\new-hb-whole embryo\wt01\segmentation_lineage.mat');
convexImageFixed_wt01 = convex_image(round(0.6*size(centroids_x,1)),AA);
for i = 1:length(convexImageFixed_wt01)
    xL_wt01(i) = size(convexImageFixed_wt01{i},1);%x axis length
end
xLAve_wt01 = mean(xL_wt01);

load('D:\LIM LAB\new-hb-whole embryo\wt02\segmentation_lineage.mat');
convexImageFixed_wt02 = convex_image(round(0.6*size(centroids_x,1)),AA);
for i = 1:length(convexImageFixed_wt02)
    xL_wt02(i) = size(convexImageFixed_wt02{i},1);%x axis length
end
xLAve_wt02 = mean(xL_wt02);

load('D:\LIM LAB\new-hb-whole embryo\wt03\segmentation_lineage.mat');
convexImageFixed_wt03 = convex_image(round(0.6*size(centroids_x,1)),AA);
for i = 1:length(convexImageFixed_wt03)
    xL_wt03(i) = size(convexImageFixed_wt03{i},1);%x axis length
end
xLAve_wt03 = mean(xL_wt03);

load('D:\LIM LAB\new-kr-whole embryo\wt01\segmentation_lineage.mat');
krconvexImageFixed_wt01 = convex_image(round(0.6*size(centroids_x,1)),AA);
for i = 1:length(krconvexImageFixed_wt01)
    krxL_wt01(i) = size(krconvexImageFixed_wt01{i},1);%x axis length
end
krxLAve_wt01 = mean(krxL_wt01);

load('D:\LIM LAB\new-kr-whole embryo\wt02\segmentation_lineage.mat');
krconvexImageFixed_wt02 = convex_image(round(0.6*size(centroids_x,1)),AA);
for i = 1:length(krconvexImageFixed_wt02)
    krxL_wt02(i) = size(krconvexImageFixed_wt02{i},1);%x axis length
end
krxLAve_wt02 = mean(krxL_wt02);

load('D:\LIM LAB\new-kr-whole embryo\wt03\segmentation_lineage.mat');
krconvexImageFixed_wt03 = convex_image(round(0.6*size(centroids_x,1)),AA);
for i = 1:length(krconvexImageFixed_wt03)
    krxL_wt03(i) = size(krconvexImageFixed_wt03{i},1);%x axis length
end
krxLAve_wt03 = mean(krxL_wt03);

load('D:\LIM LAB\new-kr-whole embryo\wt04\segmentation_lineage.mat');
krconvexImageFixed_wt04 = convex_image(round(0.6*size(centroids_x,1)),AA);
for i = 1:length(krconvexImageFixed_wt04)
    krxL_wt04(i) = size(krconvexImageFixed_wt04{i},1);%x axis length
end
krxLAve_wt04 = mean(krxL_wt04);

load('D:\LIM LAB\new-kr-whole embryo\het01\segmentation_lineage.mat');
krconvexImageFixed_het01 = convex_image(round(0.6*size(centroids_x,1)),AA);
for i = 1:length(krconvexImageFixed_het01)
    krxL_het01(i) = size(krconvexImageFixed_het01{i},1);%x axis length
end
krxLAve_het01 = mean(krxL_het01);

load('D:\LIM LAB\new-kr-whole embryo\het02\segmentation_lineage.mat');
krconvexImageFixed_het02 = convex_image(round(0.6*size(centroids_x,1)),AA);
for i = 1:length(krconvexImageFixed_het02)
    krxL_het02(i) = size(krconvexImageFixed_het02{i},1);%x axis length
end
krxLAve_het02 = mean(krxL_het02);

load('D:\LIM LAB\new-kr-whole embryo\het03\segmentation_lineage.mat');
krconvexImageFixed_het03 = convex_image(round(0.6*size(centroids_x,1)),AA);
for i = 1:length(krconvexImageFixed_het03)
    krxL_het03(i) = size(krconvexImageFixed_het03{i},1);%x axis length
end
krxLAve_het03 = mean(krxL_het03);


load('D:\LIM LAB\new-kr-whole embryo\het05\segmentation_lineage.mat');
krconvexImageFixed_het04 = convex_image(round(0.6*size(centroids_x,1)),AA);
for i = 1:length(krconvexImageFixed_het04)
    krxL_het04(i) = size(krconvexImageFixed_het04{i},1);%x axis length
end
krxLAve_het04 = mean(krxL_het04);

load('D:\LIM LAB\new-kr-whole embryo\het06\segmentation_lineage.mat');
krconvexImageFixed_het05 = convex_image(round(0.6*size(centroids_x,1)),AA);
for i = 1:length(krconvexImageFixed_het05)
    krxL_het05(i) = size(krconvexImageFixed_het05{i},1);%x axis length
end
krxLAve_het05 = mean(krxL_het05);

%% region prop
%wt
rpo_wt01 = regionprops(BW_wt01,'Extrema','Area','MajorAxisLength');
rpo_wt02 = regionprops(BW_wt02,'Extrema','Area','MajorAxisLength');
rpo_wt03 = regionprops(BW_wt03,'Extrema','Area','MajorAxisLength');
krrpo_wt01 = regionprops(krBW_wt01,'Extrema','Area','MajorAxisLength');
krrpo_wt02 = regionprops(krBW_wt02,'Extrema','Area','MajorAxisLength');
krrpo_wt03 = regionprops(krBW_wt03,'Extrema','Area','MajorAxisLength');
krrpo_wt04 = regionprops(krBW_wt04,'Extrema','Area','MajorAxisLength');

%het
krrpo_het01 = regionprops(krBW_het01,'Extrema','Area','MajorAxisLength');
krrpo_het02 = regionprops(krBW_het02,'Extrema','Area','MajorAxisLength');
krrpo_het03 = regionprops(krBW_het03,'Extrema','Area','MajorAxisLength');
krrpo_het04 = regionprops(krBW_het04,'Extrema','Area','MajorAxisLength');
krrpo_het05 = regionprops(krBW_het05,'Extrema','Area','MajorAxisLength');
%% get extrema 
extrema_wt01 = rpo_wt01.Extrema;
tip_wt01 = [min(extrema_wt01(:,1)) max(extrema_wt01(:,1))];
extrema_wt02 = rpo_wt02.Extrema;
tip_wt02 = [min(extrema_wt02(:,1)) max(extrema_wt02(:,1))];
extrema_wt03 = rpo_wt03.Extrema;
tip_wt03 = [min(extrema_wt03(:,1)) max(extrema_wt03(:,1))];
krextrema_wt01 = krrpo_wt01.Extrema;
krtip_wt01 = [min(krextrema_wt01(:,1)) max(krextrema_wt01(:,1))];
krextrema_wt02 = krrpo_wt02.Extrema;
krtip_wt02 = [min(krextrema_wt02(:,1)) max(krextrema_wt02(:,1))];
krextrema_wt03 = krrpo_wt03.Extrema;
krtip_wt03 = [min(krextrema_wt03(:,1)) max(krextrema_wt03(:,1))];
krextrema_wt04 = krrpo_wt04.Extrema;
krtip_wt04 = [min(krextrema_wt04(:,1)) max(krextrema_wt04(:,1))];

%het
krextrema_het01 = krrpo_het01.Extrema;
krtip_het01 = [min(krextrema_het01(:,1)) max(krextrema_het01(:,1))];
krextrema_het02 = krrpo_het02.Extrema;
krtip_het02 = [min(krextrema_het02(:,1)) max(krextrema_het02(:,1))];
krextrema_het03 = krrpo_het03.Extrema;
krtip_het03 = [min(krextrema_het03(:,1)) max(krextrema_het03(:,1))];
krextrema_het04 = krrpo_het04.Extrema;
krtip_het04 = [min(krextrema_het04(:,1)) max(krextrema_het04(:,1))];
krextrema_het05 = krrpo_het05.Extrema;
krtip_het05 = [min(krextrema_het05(:,1)) max(krextrema_het05(:,1))];

%% AP length
length_wt01 = abs(tip_wt01(1)-tip_wt01(2));
length_wt02 = abs(tip_wt02(1)-tip_wt02(2));
length_wt03 = abs(tip_wt03(1)-tip_wt03(2));

krlength_wt01 = abs(krtip_wt01(1)-krtip_wt01(2));
krlength_wt02 = abs(krtip_wt02(1)-krtip_wt02(2));
krlength_wt03 = abs(krtip_wt03(1)-krtip_wt03(2));
krlength_wt04 = abs(krtip_wt04(1)-krtip_wt04(2));

krlength_het01 = abs(krtip_het01(1)-krtip_het01(2));
krlength_het02 = abs(krtip_het02(1)-krtip_het02(2));
krlength_het03 = abs(krtip_het03(1)-krtip_het03(2));
krlength_het04 = abs(krtip_het04(1)-krtip_het04(2));
krlength_het05 = abs(krtip_het05(1)-krtip_het05(2));
%% calculate average width of nuclei (%EL) of each embryo
wEL_wt01 = xLAve_wt01/length_wt01*100;
wEL_wt02 = xLAve_wt02/length_wt02*100;
wEL_wt03 = xLAve_wt03/length_wt03*100;

krwEL_wt01 = krxLAve_wt01/krlength_wt01*100;
krwEL_wt02 = krxLAve_wt02/krlength_wt02*100;
krwEL_wt03 = krxLAve_wt03/krlength_wt03*100;
krwEL_wt04 = krxLAve_wt04/krlength_wt04*100;

krwEL_het01 = krxLAve_het01/krlength_het01*100;
krwEL_het02 = krxLAve_het02/krlength_het02*100;
krwEL_het03 = krxLAve_het03/krlength_het03*100;
krwEL_het04 = krxLAve_het04/krlength_het04*100;
krwEL_het05 = krxLAve_het05/krlength_het05*100;

aveEL = mean([wEL_wt01 wEL_wt02 wEL_wt03 ...
    krwEL_wt01 krwEL_wt02 krwEL_wt03 krwEL_wt04 ...
    krwEL_het01 krwEL_het02 krwEL_het03 krwEL_het04 krwEL_het05]);

avePx = mean([xLAve_wt01 xLAve_wt02 xLAve_wt03 ...
    krxLAve_wt01 krxLAve_wt02 krxLAve_wt03 krxLAve_wt04 ...
    krxLAve_het01 krxLAve_het02 krxLAve_het03 krxLAve_het04 krxLAve_het05]);%average length in pixel


%% measured by max projected image
load('embryo_length.mat');

% AP length
figure(2);
boxplot([lengthRaw_wt/1.8782 lengthRaw_het/1.8782],[ones(1,length(lengthRaw_wt)) 2*ones(1,length(lengthRaw_het))],'Labels', {'\itwt','\itKr^{1}/+'});hold on
set(gca, 'TickLabelInterpreter', 'tex');
    set(findobj(gca,'type','line'),'linew',2);
s1 = swarmchart(ones(1,length(lengthRaw_wt/1.8782)),lengthRaw_wt/1.8782,10,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);hold on
s1.XJitterWidth = 0.3;
s2 = swarmchart(2*ones(1,length(lengthRaw_het/1.8782)),lengthRaw_het/1.8782,10,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
s2.XJitterWidth = 0.3;
ylim([400 600]);

tmp1 = sort(lengthRaw_wt,'descend');
tmp2 = sort(lengthRaw_het,'descend');
[a,b,c,d] = ttest2(tmp1(1:end-1),tmp2(1:end-1));

% DV length
figure(3);
boxplot([widthRaw_wt/1.8782 widthRaw_het/1.8782],[ones(1,length(widthRaw_wt)) 2*ones(1,length(widthRaw_het))],'Labels', {'\itwt','\itKr^{1}/+'});hold on
set(gca, 'TickLabelInterpreter', 'tex');
    set(findobj(gca,'type','line'),'linew',2);
s1 = swarmchart(ones(1,length(widthRaw_wt/1.8782)),widthRaw_wt/1.8782,10,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);hold on
s1.XJitterWidth = 0.3;
s2 = swarmchart(2*ones(1,length(widthRaw_het/1.8782)),widthRaw_het/1.8782,10,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
s2.XJitterWidth = 0.3;
ylim([100 250]);

[a,b,c,d] = ttest2(widthRaw_wt,widthRaw_het);

%% get area of nuclei in each embyro
rpoSeg_wt01 = regionprops(seg_wt01,'Area');
rpoSeg_wt02 = regionprops(seg_wt02,'Area');
rpoSeg_wt03 = regionprops(seg_wt03,'Area');
krrpoSeg_wt01 = regionprops(krseg_wt01,'Area');
krrpoSeg_wt02 = regionprops(krseg_wt02,'Area');
krrpoSeg_wt03 = regionprops(krseg_wt03,'Area');
krrpoSeg_wt04 = regionprops(krseg_wt04,'Area');
krrpoSeg_het01 = regionprops(krseg_het01,'Area');
krrpoSeg_het02 = regionprops(krseg_het02,'Area');
krrpoSeg_het03 = regionprops(krseg_het03,'Area');
krrpoSeg_het04 = regionprops(krseg_het04,'Area');
krrpoSeg_het05 = regionprops(krseg_het05,'Area');

nucleiSize_wt = [mean([rpoSeg_wt01.Area]) mean([rpoSeg_wt02.Area]) mean([rpoSeg_wt03.Area]) ...
    mean([krrpoSeg_wt01.Area]) mean([krrpoSeg_wt02.Area]) mean([krrpoSeg_wt03.Area]) mean([krrpoSeg_wt04.Area])];
nucleiSize_het = [mean([krrpoSeg_het01.Area]) mean([krrpoSeg_het02.Area]) mean([krrpoSeg_het03.Area]) mean([krrpoSeg_het04.Area]) mean([krrpoSeg_het05.Area])];


%% average length of NC14 in min
NC14_het = [size(krnuc_het01,2) size(krnuc_het02,2) size(krnuc_het03,2) size(krnuc_het04,2) size(krnuc_het05,2)];
NC14_wt = [size(nuc_wt01,2) size(nuc_wt02,2) size(nuc_wt03,2)... 
    size(krnuc_wt01,2) size(krnuc_wt02,2) size(krnuc_wt03,2) size(krnuc_wt04,2)];

NC14_het = NC14_het - 1;
NC14_wt = NC14_wt - 1;

[a,b,c,d] = ttest2(NC14_wt, NC14_het);

NC14AveLength = mean([NC14_wt NC14_het])*1.02;
