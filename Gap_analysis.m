clear
close all

%% load files and set parameters
% trajectories used in this script is smoothed
load('dataset_whole_kr_kniMS2_NC14_100f.mat');

R = 75;
thr = 1000;
binsize = 0.02;
bincount = 1/binsize;
checkFlip = 0;

%% smooth trajectory


cNaN_wt01 = length(find(isnan(M_wt01(:,1))==1));
cNaN_wt02 = length(find(isnan(M_wt02(:,1))==1));
cNaN_wt03 = length(find(isnan(M_wt03(:,1))==1));
cNaN_het01 = length(find(isnan(M_het01(:,1))==1));
cNaN_het02 = length(find(isnan(M_het02(:,1))==1));
cNaN_het03 = length(find(isnan(M_het03(:,1))==1));



for i = 1:size(M_wt01,2)
    M_wt01(cNaN_wt01+1:end,i) = smooth(M_wt01(cNaN_wt01+1:end,i));
end
for i = 1:size(M_wt02,2)
    M_wt02(cNaN_wt02+1:end,i) = smooth(M_wt02(cNaN_wt02+1:end,i));
end
for i = 1:size(M_wt03,2)
    M_wt03(cNaN_wt03+1:end,i) = smooth(M_wt03(cNaN_wt03+1:end,i));
end

for i = 1:size(M_het01,2)
    M_het01(cNaN_het01+1:end,i) = smooth(M_het01(cNaN_het01+1:end,i));
end
for i = 1:size(M_het02,2)
    M_het02(cNaN_het02+1:end,i) = smooth(M_het02(cNaN_het02+1:end,i));
end
for i = 1:size(M_het03,2)
    M_het03(cNaN_het03+1:end,i) = smooth(M_het03(cNaN_het03+1:end,i));
end


%% extract max intensity

for i=1:size(M_wt01,2) % all nuclei
  maxM_wt01(i) = max(M_wt01(:,i)); % obtain max intensity for each nuclues
  maxM2_wt01(i) = max(M_wt01(round(0.75*size(M_wt01,1)):end,i));%max intensity for last 1/4 frames
end 

for i=1:size(M_wt02,2) % all nuclei
  maxM_wt02(i) = max(M_wt02(:,i)); % obtain max intensity for each nuclues
  maxM2_wt02(i) = max(M_wt02(round(0.75*size(M_wt02,1)):end,i));%max intensity for last 1/4 frames
end 

for i=1:size(M_wt03,2) % all nuclei
  maxM_wt03(i) = max(M_wt03(:,i)); % obtain max intensity for each nuclues
  maxM2_wt03(i) = max(M_wt03(round(0.75*size(M_wt03,1)):end,i));%max intensity for last 1/4 frames
end 

for i=1:size(M_het01,2) % all nuclei
  maxM_het01(i) = max(M_het01(:,i)); % obtain max intensity for each nuclues
  maxM2_het01(i) = max(M_het01(round(0.75*size(M_het01,1)):end,i));%max intensity for last 1/4 frames
end 

for i=1:size(M_het02,2) % all nuclei
  maxM_het02(i) = max(M_het02(:,i)); % obtain max intensity for each nuclues
  maxM2_het02(i) = max(M_het02(round(0.75*size(M_het02,1)):end,i));%max intensity for last 1/4 frames
end 

for i=1:size(M_het03,2) % all nuclei
  maxM_het03(i) = max(M_het03(:,i)); % obtain max intensity for each nuclues
  maxM2_het03(i) = max(M_het03(round(0.75*size(M_het03,1)):end,i));%max intensity for last 1/4 frames
end 

%% active nuclei

act_wt01 = find(maxM_wt01 > thr);
actLate_wt01 = find(maxM2_wt01 > thr);
act_wt02 = find(maxM_wt02 > thr);
actLate_wt02 = find(maxM2_wt02 > thr);
act_wt03 = find(maxM_wt03 > thr);
actLate_wt03 = find(maxM2_wt03 > thr);

act_het01 = find(maxM_het01 > thr);
actLate_het01 = find(maxM2_het01 > thr);
act_het02 = find(maxM_het02 > thr);
actLate_het02 = find(maxM2_het02 > thr);
act_het03 = find(maxM_het03 > thr);
actLate_het03 = find(maxM2_het03 > thr);


%% get regionprops - use Extrema to calibrate EL
rpo_wt01 = regionprops(BW_wt01,'Extrema');
rpo_wt02 = regionprops(BW_wt02,'Extrema');
rpo_wt03 = regionprops(BW_wt03,'Extrema');

rpo_het01 = regionprops(BW_het01,'Extrema');
rpo_het02 = regionprops(BW_het02,'Extrema');
rpo_het03 = regionprops(BW_het03,'Extrema');

%% extract Extrema get extreme values

extrema_wt01 = rpo_wt01.Extrema;
tip_wt01 = [min(extrema_wt01(:,1)) max(extrema_wt01(:,1))];
extrema_wt02 = rpo_wt02.Extrema;
tip_wt02 = [min(extrema_wt02(:,1)) max(extrema_wt02(:,1))];
extrema_wt03 = rpo_wt03.Extrema;
tip_wt03 = [min(extrema_wt03(:,1)) max(extrema_wt03(:,1))];

extrema_het01 = rpo_het01.Extrema;
tip_het01 = [min(extrema_het01(:,1)) max(extrema_het01(:,1))];
extrema_het02 = rpo_het02.Extrema;
tip_het02 = [min(extrema_het02(:,1)) max(extrema_het02(:,1))];
extrema_het03 = rpo_het03.Extrema;
tip_het03 = [min(extrema_het03(:,1)) max(extrema_het03(:,1))];

%% flip embryos - note: tips need to need min and max for anterior and posterior

%flip nuclei
lineageCx_wt01 = abs(lineageCx_wt01 - 950);
lineageCx_het02 = abs(lineageCx_het02 - 950);
lineageCx_het03 = abs(lineageCx_het03 - 950);

%flip nuclei fixed cx
cx_wt01 = abs(cx_wt01 - 950); %flip x axis
cx_het02 = abs(cx_het02 - 950); %flip x axis
cx_het03 = abs(cx_het03 - 950); %flip x axis

%flip tip points
tip_wt01 = abs(tip_wt01 - 1400); %flip x axis
tip_het02 = abs(tip_het02 - 1400); %flip x axis
tip_het03 = abs(tip_het03 - 1400); %flip x axis

%flip nuc image
BW_wt01 = flip(BW_wt01,2);
BW_het02 = flip(BW_het02,2);
BW_het03 = flip(BW_het03,2);

%flip seg image
seg_wt01 = flip(seg_wt01,2);
seg_het02 = flip(seg_het02,2);
seg_het03 = flip(seg_het03,2);


%% calibrate x position by EL
cxEL_wt01 = (cx_wt01+225 - min(tip_wt01))/(max(tip_wt01)- min(tip_wt01));
cxEL_wt02 = (cx_wt02+225 - min(tip_wt02))/(max(tip_wt02)- min(tip_wt02));
cxEL_wt03 = (cx_wt03+225 - min(tip_wt03))/(max(tip_wt03)- min(tip_wt03));

cxEL_het01 = (cx_het01+225 - min(tip_het01))/(max(tip_het01)- min(tip_het01));
cxEL_het02 = (cx_het02+225 - min(tip_het02))/(max(tip_het02)- min(tip_het02));
cxEL_het03 = (cx_het03+225 - min(tip_het03))/(max(tip_het03)- min(tip_het03));
%% get the midline of each embryo
midline_wt01 = (max(cy_wt01) + min(cy_wt01))/2;
midline_wt02 = (max(cy_wt02) + min(cy_wt02))/2;
midline_wt03 = (max(cy_wt03) + min(cy_wt03))/2;
midline_het01 = (max(cy_het01) + min(cy_het01))/2;
midline_het02 = (max(cy_het02) + min(cy_het02))/2;
midline_het03 = (max(cy_het03) + min(cy_het03))/2;
%% bin by 2% EL

%wt01
for i = 1:50
    binEL_wt01{i} = find(cxEL_wt01(:) > (i-1)*binsize & cxEL_wt01(:) <= i*binsize & cy_wt01(:) >= midline_wt01-R & cy_wt01(:) <= midline_wt01+R);
end
%wt02
for i = 1:50
    binEL_wt02{i} = find(cxEL_wt02(:) > (i-1)*binsize & cxEL_wt02(:) <= i*binsize & cy_wt02(:) >= midline_wt02-R & cy_wt02(:) <= midline_wt02+R);
end
%wt03
for i = 1:50
    binEL_wt03{i} = find(cxEL_wt03(:) > (i-1)*binsize & cxEL_wt03(:) <= i*binsize & cy_wt03(:) >= midline_wt03-R & cy_wt03(:) <= midline_wt03+R);
end

%het01
for i = 1:50
    binEL_het01{i} = find(cxEL_het01(:) > (i-1)*binsize & cxEL_het01(:) <= i*binsize & cy_het01(:) >= midline_het01-R & cy_het01(:) <= midline_het01+R);
end
%het02
for i = 1:50
    binEL_het02{i} = find(cxEL_het02(:) > (i-1)*binsize & cxEL_het02(:) <= i*binsize & cy_het02(:) >= midline_het02-R & cy_het02(:) <= midline_het02+R);
end
%het03
for i = 1:50
    binEL_het03{i} = find(cxEL_het03(:) > (i-1)*binsize & cxEL_het03(:) <= i*binsize & cy_het03(:) >= midline_het03-R & cy_het03(:) <= midline_het03+R);
end


%% mean bin trajectory

%get bin amp for each embryo
for i = 1:100
    for j = 1:bincount
        binAveAmp_wt01{i,j} = mean(M_wt01(i,binEL_wt01{j}));
        binAveAmp_wt02{i,j} = mean(M_wt02(i,binEL_wt02{j}));
        binAveAmp_wt03{i,j} = mean(M_wt03(i,binEL_wt03{j}));

        binAveAmp_het01{i,j} = mean(M_het01(i,binEL_het01{j}));
        binAveAmp_het02{i,j} = mean(M_het02(i,binEL_het02{j}));
        binAveAmp_het03{i,j} = mean(M_het03(i,binEL_het03{j}));

    end
end

%bin amp for average embryo + errorbar
for i = 1:100
    for j = 1:bincount
        binAveAmp_wt{i,j} = mean([binAveAmp_wt01{i,j}, binAveAmp_wt02{i,j}, binAveAmp_wt03{i,j}]);
        binAveAmp_het{i,j} = mean([binAveAmp_het01{i,j}, binAveAmp_het02{i,j}, binAveAmp_het03{i,j}]);

        binAveAmpSE_wt{i,j} = std([binAveAmp_wt01{i,j}, binAveAmp_wt02{i,j}, binAveAmp_wt03{i,j}])/sqrt(3);
        binAveAmpSE_het{i,j} = std([binAveAmp_het01{i,j}, binAveAmp_het02{i,j}, binAveAmp_het03{i,j}])/sqrt(3);
    end
end

%plot trajectory
figure(1064);
plot(2:2:100,[binAveAmp_wt{52,:}]/1000,'-k','LineWidth',3);hold on
plot(2:2:100,[binAveAmp_het{52,:}]/1000,'-r','LineWidth',3);
shadedErrorBar2(2:2:100,[binAveAmp_wt{52,:}]/1000, [binAveAmpSE_wt{52,:}]/1000, 'lineprops', {'-k','LineWidth',2});hold on
shadedErrorBar2(2:2:100,[binAveAmp_het{52,:}]/1000,[binAveAmpSE_het{52,:}]/1000, 'lineprops', {'-r','LineWidth',2});
legend('\itwt','\itKr^{1}/+');
ylim([0 4.5]); xlim([40 80]); hold off
