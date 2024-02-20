clear 
close all


%% load and set


load('D:\LIM LAB\kr_project_code_cleanup\dataset_whole_hb_NC14_100f.mat');
load('D:\LIM LAB\kr_project_code_cleanup\dataset_whole_kr_NC14_100f.mat');
load('all_stripes_nuclei_index_NC14_100f.mat');

bincount = 50;

thr = 500;
thr2 = 250;
%% smooth trajectory

cNaN_wt01 = length(find(isnan(M_wt01(:,1))==1));
cNaN_wt02 = length(find(isnan(M_wt02(:,1))==1));
cNaN_wt03 = length(find(isnan(M_wt03(:,1))==1));
krcNaN_wt01 = length(find(isnan(krM_wt01(:,1))==1));
krcNaN_wt02 = length(find(isnan(krM_wt02(:,1))==1));
krcNaN_wt03 = length(find(isnan(krM_wt03(:,1))==1));
krcNaN_wt04 = length(find(isnan(krM_wt04(:,1))==1));
krcNaN_het01 = length(find(isnan(krM_het01(:,1))==1));
krcNaN_het02 = length(find(isnan(krM_het02(:,1))==1));
krcNaN_het03 = length(find(isnan(krM_het03(:,1))==1));
krcNaN_het04 = length(find(isnan(krM_het04(:,1))==1));
krcNaN_het05 = length(find(isnan(krM_het05(:,1))==1));


for i = 1:size(M_wt01,2)
    M_wt01(cNaN_wt01+1:end,i) = smooth(M_wt01(cNaN_wt01+1:end,i));
end
for i = 1:size(M_wt02,2)
    M_wt02(cNaN_wt02+1:end,i) = smooth(M_wt02(cNaN_wt02+1:end,i));
end
for i = 1:size(M_wt03,2)
    M_wt03(cNaN_wt03+1:end,i) = smooth(M_wt03(cNaN_wt03+1:end,i));
end
for i = 1:size(krM_wt01,2)
    krM_wt01(krcNaN_wt01+1:end,i) = smooth(krM_wt01(krcNaN_wt01+1:end,i));
end
for i = 1:size(krM_wt02,2)
    krM_wt02(krcNaN_wt02+1:end,i) = smooth(krM_wt02(krcNaN_wt02+1:end,i));
end
for i = 1:size(krM_wt03,2)
    krM_wt03(krcNaN_wt03+1:end,i) = smooth(krM_wt03(krcNaN_wt03+1:end,i));
end
for i = 1:size(krM_wt01,2)
    krM_wt01(krcNaN_wt01+1:end,i) = smooth(krM_wt01(krcNaN_wt01+1:end,i));
end
for i = 1:size(krM_wt04,2)
    krM_wt04(krcNaN_wt04+1:end,i) = smooth(krM_wt04(krcNaN_wt04+1:end,i));
end
for i = 1:size(krM_het01,2)
    krM_het01(krcNaN_het01+1:end,i) = smooth(krM_het01(krcNaN_het01+1:end,i));
end
for i = 1:size(krM_het02,2)
    krM_het02(krcNaN_het02+1:end,i) = smooth(krM_het02(krcNaN_het02+1:end,i));
end
for i = 1:size(krM_het03,2)
    krM_het03(krcNaN_het03+1:end,i) = smooth(krM_het03(krcNaN_het03+1:end,i));
end
for i = 1:size(krM_het01,2)
    krM_het01(krcNaN_het01+1:end,i) = smooth(krM_het01(krcNaN_het01+1:end,i));
end
for i = 1:size(krM_het04,2)
    krM_het04(krcNaN_het04+1:end,i) = smooth(krM_het04(krcNaN_het04+1:end,i));
end
for i = 1:size(krM_het05,2)
    krM_het05(krcNaN_het05+1:end,i) = smooth(krM_het05(krcNaN_het05+1:end,i));
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

for i=1:size(krM_wt01,2) % all nuclei
  krmaxM_wt01(i) = max(krM_wt01(:,i)); % obtain max intensity for each nuclues
  krmaxM2_wt01(i) = max(krM_wt01(round(0.75*size(krM_wt01,1)):end,i));%max intensity for last 1/4 frames
end 

for i=1:size(krM_wt02,2) % all nuclei
  krmaxM_wt02(i) = max(krM_wt02(:,i)); % obtain max intensity for each nuclues
  krmaxM2_wt02(i) = max(krM_wt02(round(0.75*size(krM_wt02,1)):end,i));%max intensity for last 1/4 frames
end 

for i=1:size(krM_wt03,2) % all nuclei
  krmaxM_wt03(i) = max(krM_wt03(:,i)); % obtain max intensity for each nuclues
  krmaxM2_wt03(i) = max(krM_wt03(round(0.75*size(krM_wt03,1)):end,i));%max intensity for last 1/4 frames
end 

for i=1:size(krM_wt04,2) % all nuclei
  krmaxM_wt04(i) = max(krM_wt04(:,i)); % obtain max intensity for each nuclues
  krmaxM2_wt04(i) = max(krM_wt04(round(0.75*size(krM_wt04,1)):end,i));%max intensity for last 1/4 frames
end 


for i=1:size(krM_het01,2) % all nuclei
  krmaxM_het01(i) = max(krM_het01(:,i)); % obtain max intensity for each nuclues
  krmaxM2_het01(i) = max(krM_het01(round(0.75*size(krM_het01,1)):end,i));%max intensity for last 1/4 frames
end 

for i=1:size(krM_het02,2) % all nuclei
  krmaxM_het02(i) = max(krM_het02(:,i)); % obtain max intensity for each nuclues
  krmaxM2_het02(i) = max(krM_het02(round(0.75*size(krM_het02,1)):end,i));%max intensity for last 1/4 frames
end 

for i=1:size(krM_het03,2) % all nuclei
  krmaxM_het03(i) = max(krM_het03(:,i)); % obtain max intensity for each nuclues
  krmaxM2_het03(i) = max(krM_het03(round(0.75*size(krM_het03,1)):end,i));%max intensity for last 1/4 frames
end 

for i=1:size(krM_het04,2) % all nuclei
  krmaxM_het04(i) = max(krM_het04(:,i)); % obtain max intensity for each nuclues
  krmaxM2_het04(i) = max(krM_het04(round(0.75*size(krM_het04,1)):end,i));%max intensity for last 1/4 frames
end 

for i=1:size(krM_het05,2) % all nuclei
  krmaxM_het05(i) = max(krM_het05(:,i)); % obtain max intensity for each nuclues
  krmaxM2_het05(i) = max(krM_het05(round(0.75*size(krM_het05,1)):end,i));%max intensity for last 1/4 frames
end 

%% get regionprops

rpo_wt01 = regionprops(BW_wt01,'Extrema','Area','MajorAxisLength');
rpo_wt02 = regionprops(BW_wt02,'Extrema','Area','MajorAxisLength');
rpo_wt03 = regionprops(BW_wt03,'Extrema','Area','MajorAxisLength');
krrpo_wt01 = regionprops(krBW_wt01,'Extrema','Area','MajorAxisLength');
krrpo_wt02 = regionprops(krBW_wt02,'Extrema','Area','MajorAxisLength');
krrpo_wt03 = regionprops(krBW_wt03,'Extrema','Area','MajorAxisLength');
krrpo_wt04 = regionprops(krBW_wt04,'Extrema','Area','MajorAxisLength');

krrpo_het01 = regionprops(krBW_het01,'Extrema','Area','MajorAxisLength');
krrpo_het02 = regionprops(krBW_het02,'Extrema','Area','MajorAxisLength');
krrpo_het03 = regionprops(krBW_het03,'Extrema','Area','MajorAxisLength');
krrpo_het04 = regionprops(krBW_het04,'Extrema','Area','MajorAxisLength');
krrpo_het05 = regionprops(krBW_het05,'Extrema','Area','MajorAxisLength');

%% extract Extrema get extreme values

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

fuse_wt01 = imfuse(BW_wt01,seg_wt01);
fuse_wt02 = imfuse(BW_wt02,seg_wt02);
fuse_wt03 = imfuse(BW_wt03,seg_wt03);
krfuse_wt01 = imfuse(krBW_wt01,krseg_wt01);
krfuse_wt02 = imfuse(krBW_wt02,krseg_wt02);
krfuse_wt03 = imfuse(krBW_wt03,krseg_wt03);
krfuse_wt04 = imfuse(krBW_wt04,krseg_wt04);
krfuse_het01 = imfuse(krBW_het01,krseg_het01);
krfuse_het02 = imfuse(krBW_het02,krseg_het02);
krfuse_het03 = imfuse(krBW_het03,krseg_het03);
krfuse_het04 = imfuse(krBW_het04,krseg_het04);
krfuse_het05 = imfuse(krBW_het05,krseg_het05);
%% flip embryos - note: tips need to need min and max for anterior and posterior
%flip nuclei
lineageCx_wt02 = abs(lineageCx_wt02 - 950); %flip x axis
krlineageCx_wt04 = abs(krlineageCx_wt04 - 950); %flip x axis
krlineageCx_het02 = abs(krlineageCx_het02 - 950); %flip x axis

%flip nuclei fixed cx
cx_wt02 = abs(cx_wt02 - 950); %flip x axis
krcx_wt04 = abs(krcx_wt04 - 950); %flip x axis
krcx_het02 = abs(krcx_het02 - 950); %flip x axis

%flip tip points
tip_wt02 = abs(tip_wt02 - 950); %flip x axis
krtip_wt04 = abs(krtip_wt04 - 950); %flip x axis
krtip_het02 = abs(krtip_het02 - 950); %flip x axis

%flip fused nuc image
fuseFlip_wt02 = flip(fuse_wt02,2);
krfuseFlip_wt04 = flip(krfuse_wt04,2);
krfuseFlip_het02 = flip(krfuse_het02,2);

%asign new name to embryos that don't need flip
fuseFlip_wt01 = fuse_wt01;
fuseFlip_wt03 = fuse_wt03;
krfuseFlip_wt01 = krfuse_wt01;
krfuseFlip_wt02 = krfuse_wt02;
krfuseFlip_wt03 = krfuse_wt03;
krfuseFlip_het01 = krfuse_het01;
krfuseFlip_het03 = krfuse_het03;
krfuseFlip_het04 = krfuse_het04;
krfuseFlip_het05 = krfuse_het05;

%% surface major axis length
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

%% calibrate x position by EL
cxEL_wt01 = (cx_wt01 - min(tip_wt01))/(max(tip_wt01)- min(tip_wt01));
cxEL_wt02 = (cx_wt02 - min(tip_wt02))/(max(tip_wt02)- min(tip_wt02));
cxEL_wt03 = (cx_wt03 - min(tip_wt03))/(max(tip_wt03)- min(tip_wt03));
krcxEL_wt01 = (krcx_wt01 - min(krtip_wt01))/(max(krtip_wt01)- min(krtip_wt01));
krcxEL_wt02 = (krcx_wt02 - min(krtip_wt02))/(max(krtip_wt02)- min(krtip_wt02));
krcxEL_wt03 = (krcx_wt03 - min(krtip_wt03))/(max(krtip_wt03)- min(krtip_wt03));
krcxEL_wt04 = (krcx_wt04 - min(krtip_wt04))/(max(krtip_wt04)- min(krtip_wt04));
krcxEL_het01 = (krcx_het01 - min(krtip_het01))/(max(krtip_het01)- min(krtip_het01));
krcxEL_het02 = (krcx_het02 - min(krtip_het02))/(max(krtip_het02)- min(krtip_het02));
krcxEL_het03 = (krcx_het03 - min(krtip_het03))/(max(krtip_het03)- min(krtip_het03));
krcxEL_het04 = (krcx_het04 - min(krtip_het04))/(max(krtip_het04)- min(krtip_het04));
krcxEL_het05 = (krcx_het05 - min(krtip_het05))/(max(krtip_het05)- min(krtip_het05));

%% bin y axis to 5 bins
Ybincount = 5;
cybinW_wt01 = (max(cy_wt01) - min(cy_wt01))/Ybincount;%width of cy bin
cybinW_wt02 = (max(cy_wt02) - min(cy_wt02))/Ybincount;%width of cy bin
cybinW_wt03 = (max(cy_wt03) - min(cy_wt03))/Ybincount;%width of cy bin
krcybinW_wt01 = (max(krcy_wt01) - min(krcy_wt01))/Ybincount;%width of krcy bin
krcybinW_wt02 = (max(krcy_wt02) - min(krcy_wt02))/Ybincount;%width of krcy bin
krcybinW_wt03 = (max(krcy_wt03) - min(krcy_wt03))/Ybincount;%width of krcy bin
krcybinW_wt04 = (max(krcy_wt04) - min(krcy_wt04))/Ybincount;%width of krcy bin
krcybinW_het01 = (max(krcy_het01) - min(krcy_het01))/Ybincount;%width of krcy bin
krcybinW_het02 = (max(krcy_het02) - min(krcy_het02))/Ybincount;%width of krcy bin
krcybinW_het03 = (max(krcy_het03) - min(krcy_het03))/Ybincount;%width of krcy bin
krcybinW_het04 = (max(krcy_het04) - min(krcy_het04))/Ybincount;%width of krcy bin
krcybinW_het05 = (max(krcy_het05) - min(krcy_het05))/Ybincount;%width of krcy bin

for i = 1:Ybincount
    cybin_wt01{i} = find(cy_wt01 >= min(cy_wt01) + (i-1)*cybinW_wt01 & cy_wt01 <= min(cy_wt01) + i*cybinW_wt01);
    cybin_wt02{i} = find(cy_wt02 >= min(cy_wt02) + (i-1)*cybinW_wt02 & cy_wt02 <= min(cy_wt02) + i*cybinW_wt02);
    cybin_wt03{i} = find(cy_wt03 >= min(cy_wt03) + (i-1)*cybinW_wt03 & cy_wt03 <= min(cy_wt03) + i*cybinW_wt03);
    krcybin_wt01{i} = find(krcy_wt01 >= min(krcy_wt01) + (i-1)*krcybinW_wt01 & krcy_wt01 <= min(krcy_wt01) + i*krcybinW_wt01);
    krcybin_wt02{i} = find(krcy_wt02 >= min(krcy_wt02) + (i-1)*krcybinW_wt02 & krcy_wt02 <= min(krcy_wt02) + i*krcybinW_wt02);
    krcybin_wt03{i} = find(krcy_wt03 >= min(krcy_wt03) + (i-1)*krcybinW_wt03 & krcy_wt03 <= min(krcy_wt03) + i*krcybinW_wt03);
    krcybin_wt04{i} = find(krcy_wt04 >= min(krcy_wt04) + (i-1)*krcybinW_wt04 & krcy_wt04 <= min(krcy_wt04) + i*krcybinW_wt04);
    krcybin_het01{i} = find(krcy_het01 >= min(krcy_het01) + (i-1)*krcybinW_het01 & krcy_het01 <= min(krcy_het01) + i*krcybinW_het01);
    krcybin_het02{i} = find(krcy_het02 >= min(krcy_het02) + (i-1)*krcybinW_het02 & krcy_het02 <= min(krcy_het02) + i*krcybinW_het02);
    krcybin_het03{i} = find(krcy_het03 >= min(krcy_het03) + (i-1)*krcybinW_het03 & krcy_het03 <= min(krcy_het03) + i*krcybinW_het03);
    krcybin_het04{i} = find(krcy_het04 >= min(krcy_het04) + (i-1)*krcybinW_het04 & krcy_het04 <= min(krcy_het04) + i*krcybinW_het04);
    krcybin_het05{i} = find(krcy_het05 >= min(krcy_het05) + (i-1)*krcybinW_het05 & krcy_het05 <= min(krcy_het05) + i*krcybinW_het05);
end


%% identify max and min cxEL per bin per frame in 2nd half NC14
D = 30;
for i = 1:100
    tmp1 = find(M_wt01(i,:) >= thr);
    S5actOT_wt01{i} = intersect(tmp1, stripe5_wt01); 

    for j = 1:5
        tmp3 = intersect(cybin_wt01{j}, S5actOT_wt01{i});
        if length(tmp3) > 2
            tmp4 = sort(cx_wt01(tmp3),'descend');
            if abs(tmp4(1) - tmp4(2))>D
                S5A_wt01(i,j) = min(cx_wt01(tmp3));%S5 anterior
                S5P_wt01(i,j) = tmp4(2);%S5 posterior
            else
                S5A_wt01(i,j) = min(cx_wt01(tmp3));%S5 anterior
                S5P_wt01(i,j) = max(cx_wt01(tmp3));%S5 posterior
            end
        elseif length(tmp3) == 2
            S5A_wt01(i,j) = min(cx_wt01(tmp3));%S5 anterior
            S5P_wt01(i,j) = max(cx_wt01(tmp3));%S5 posterior
        end 
    end
end
for i = 1:100
    tmp1 = find(M_wt02(i,:) >= thr);
    S5actOT_wt02{i} = intersect(tmp1, stripe5_wt02); 

    for j = 1:5
        tmp3 = intersect(cybin_wt02{j}, S5actOT_wt02{i});
        if length(tmp3) > 2
            tmp4 = sort(cx_wt02(tmp3),'descend');
            if abs(tmp4(1) - tmp4(2))>D
                S5A_wt02(i,j) = min(cx_wt02(tmp3));%S5 anterior
                S5P_wt02(i,j) = tmp4(2);%S5 posterior
            else
                S5A_wt02(i,j) = min(cx_wt02(tmp3));%S5 anterior
                S5P_wt02(i,j) = max(cx_wt02(tmp3));%S5 posterior
            end
        elseif length(tmp3) == 2
            S5A_wt02(i,j) = min(cx_wt02(tmp3));%S5 anterior
            S5P_wt02(i,j) = max(cx_wt02(tmp3));%S5 posterior
        end 
    end
end
for i = 1:100
    tmp1 = find(M_wt03(i,:) >= thr);
    S5actOT_wt03{i} = intersect(tmp1, stripe5_wt03); 

    for j = 1:5
        tmp3 = intersect(cybin_wt03{j}, S5actOT_wt03{i});
        if length(tmp3) > 2
            tmp4 = sort(cx_wt03(tmp3),'descend');
            if abs(tmp4(1) - tmp4(2))>D
                S5A_wt03(i,j) = min(cx_wt03(tmp3));%S5 anterior
                S5P_wt03(i,j) = tmp4(2);%S5 posterior
            else
                S5A_wt03(i,j) = min(cx_wt03(tmp3));%S5 anterior
                S5P_wt03(i,j) = max(cx_wt03(tmp3));%S5 posterior
            end
        elseif length(tmp3) == 2
            S5A_wt03(i,j) = min(cx_wt03(tmp3));%S5 anterior
            S5P_wt03(i,j) = max(cx_wt03(tmp3));%S5 posterior
        end 
    
    end
end
for i = 1:100
    tmp1 = find(krM_wt01(i,:) >= thr);
    krS5actOT_wt01{i} = intersect(tmp1, krstripe5_wt01); 
    for j = 1:5
        tmp3 = intersect(krcybin_wt01{j}, krS5actOT_wt01{i});
        if length(tmp3) > 2
            tmp4 = sort(krcx_wt01(tmp3),'descend');
            if abs(tmp4(1) - tmp4(2))>D
                krS5A_wt01(i,j) = min(krcx_wt01(tmp3));%S5 anterior
                krS5P_wt01(i,j) = tmp4(2);%S5 posterior
            else
                krS5A_wt01(i,j) = min(krcx_wt01(tmp3));%krS5 anterior
                krS5P_wt01(i,j) = max(krcx_wt01(tmp3));%krS5 posterior
            end
        elseif length(tmp3) == 2
            krS5A_wt01(i,j) = min(krcx_wt01(tmp3));%S5 anterior
            krS5P_wt01(i,j) = max(krcx_wt01(tmp3));%S5 posterior
        end 
    
    end
end
for i = 1:100
    tmp1 = find(krM_wt02(i,:) >= thr);
    krS5actOT_wt02{i} = intersect(tmp1, krstripe5_wt02); 
    for j = 1:5
        tmp3 = intersect(krcybin_wt02{j}, krS5actOT_wt02{i});
        if length(tmp3) > 2
            tmp4 = sort(krcx_wt02(tmp3),'descend');
            if abs(tmp4(1) - tmp4(2))>D
                krS5A_wt02(i,j) = min(krcx_wt02(tmp3));%S5 anterior
                krS5P_wt02(i,j) = tmp4(2);%S5 posterior
            else
                krS5A_wt02(i,j) = min(krcx_wt02(tmp3));%krS5 anterior
                krS5P_wt02(i,j) = max(krcx_wt02(tmp3));%krS5 posterior
            end
        elseif length(tmp3) == 2
            krS5A_wt02(i,j) = min(krcx_wt02(tmp3));%S5 anterior
            krS5P_wt02(i,j) = max(krcx_wt02(tmp3));%S5 posterior
        end 
    end
end
for i = 1:100
    tmp1 = find(krM_wt03(i,:) >= thr);
    krS5actOT_wt03{i} = intersect(tmp1, krstripe5_wt03); 
    for j = 1:5
        tmp3 = intersect(krcybin_wt03{j}, krS5actOT_wt03{i});
        if length(tmp3) > 2
            tmp4 = sort(krcx_wt03(tmp3),'descend');
            if abs(tmp4(1) - tmp4(2))>D
                krS5A_wt03(i,j) = min(krcx_wt03(tmp3));%S5 anterior
                krS5P_wt03(i,j) = tmp4(2);%S5 posterior
            else
                krS5A_wt03(i,j) = min(krcx_wt03(tmp3));%krS5 anterior
                krS5P_wt03(i,j) = max(krcx_wt03(tmp3));%krS5 posterior
            end
        elseif length(tmp3) == 2
            krS5A_wt03(i,j) = min(krcx_wt03(tmp3));%S5 anterior
            krS5P_wt03(i,j) = max(krcx_wt03(tmp3));%S5 posterior
        end 
    end
end
for i = 1:100
    tmp1 = find(krM_wt04(i,:) >= thr);
    krS5actOT_wt04{i} = intersect(tmp1, krstripe5_wt04); 
    for j = 1:5
        tmp3 = intersect(krcybin_wt04{j}, krS5actOT_wt04{i});
        if length(tmp3) > 2
            tmp4 = sort(krcx_wt04(tmp3),'descend');
            if abs(tmp4(1) - tmp4(2))>D
                krS5A_wt04(i,j) = min(krcx_wt04(tmp3));%S5 anterior
                krS5P_wt04(i,j) = tmp4(2);%S5 posterior
            else
                krS5A_wt04(i,j) = min(krcx_wt04(tmp3));%krS5 anterior
                krS5P_wt04(i,j) = max(krcx_wt04(tmp3));%krS5 posterior
            end
        elseif length(tmp3) == 2
            krS5A_wt04(i,j) = min(krcx_wt04(tmp3));%S5 anterior
            krS5P_wt04(i,j) = max(krcx_wt04(tmp3));%S5 posterior
        end 
    end
end

for i = 1:100
    tmp1 = find(krM_het01(i,:) >= thr);
    krS5actOT_het01{i} = intersect(tmp1, krstripe5_het01); 
    for j = 1:5
        tmp3 = intersect(krcybin_het01{j}, krS5actOT_het01{i});
        if length(tmp3) > 2
            tmp4 = sort(krcx_het01(tmp3),'descend');
            if abs(tmp4(1) - tmp4(2))>D
                krS5A_het01(i,j) = min(krcx_het01(tmp3));%S5 anterior
                krS5P_het01(i,j) = tmp4(2);%S5 posterior
            else
                krS5A_het01(i,j) = min(krcx_het01(tmp3));%krS5 anterior
                krS5P_het01(i,j) = max(krcx_het01(tmp3));%krS5 posterior
            end
        elseif length(tmp3) == 2
            krS5A_het01(i,j) = min(krcx_het01(tmp3));%krS5 anterior
            krS5P_het01(i,j) = max(krcx_het01(tmp3));%krS5 posterior           
        end
        
    end
end
for i = 1:100
    tmp1 = find(krM_het02(i,:) >= thr);
    krS5actOT_het02{i} = intersect(tmp1, krstripe5_het02); 
    for j = 1:5
        tmp3 = intersect(krcybin_het02{j}, krS5actOT_het02{i});
        if length(tmp3) > 2
            tmp4 = sort(krcx_het02(tmp3),'descend');
            if abs(tmp4(1) - tmp4(2))>D
                krS5A_het02(i,j) = min(krcx_het02(tmp3));%S5 anterior
                krS5P_het02(i,j) = tmp4(2);%S5 posterior
            else
                krS5A_het02(i,j) = min(krcx_het02(tmp3));%krS5 anterior
                krS5P_het02(i,j) = max(krcx_het02(tmp3));%krS5 posterior
            end
        elseif length(tmp3) == 2
                krS5A_het02(i,j) = min(krcx_het02(tmp3));%krS5 anterior
                krS5P_het02(i,j) = max(krcx_het02(tmp3));%krS5 posterior
        end 
    end

end
for i = 1:100
    tmp1 = find(krM_het03(i,:) >= thr);
    krS5actOT_het03{i} = intersect(tmp1, krstripe5_het03); 
    for j = 1:5
        tmp3 = intersect(krcybin_het03{j}, krS5actOT_het03{i});
        if length(tmp3) >=2
            tmp4 = sort(krcx_het03(tmp3),'descend');
            if abs(tmp4(1) - tmp4(2))>D
                krS5A_het03(i,j) = min(krcx_het03(tmp3));%S5 anterior
                krS5P_het03(i,j) = tmp4(2);%S5 posterior
            else
                krS5A_het03(i,j) = min(krcx_het03(tmp3));%krS5 anterior
                krS5P_het03(i,j) = max(krcx_het03(tmp3));%krS5 posterior
            end
        elseif length(tmp3) == 2
                krS5A_het03(i,j) = min(krcx_het03(tmp3));%krS5 anterior
                krS5P_het03(i,j) = max(krcx_het03(tmp3));%krS5 posterior
        end 
    end
end
for i = 1:100
    tmp1 = find(krM_het04(i,:) >= thr);
    krS5actOT_het04{i} = intersect(tmp1, krstripe5_het04); 
    for j = 1:5
        tmp3 = intersect(krcybin_het04{j}, krS5actOT_het04{i});
        if length(tmp3) >=2
            tmp4 = sort(krcx_het04(tmp3),'descend');
            if abs(tmp4(1) - tmp4(2))>D
                krS5A_het04(i,j) = min(krcx_het04(tmp3));%S5 anterior
                krS5P_het04(i,j) = tmp4(2);%S5 posterior
            else
                krS5A_het04(i,j) = min(krcx_het04(tmp3));%krS5 anterior
                krS5P_het04(i,j) = max(krcx_het04(tmp3));%krS5 posterior
            end
        elseif length(tmp3) == 2
                krS5A_het04(i,j) = min(krcx_het04(tmp3));%krS5 anterior
                krS5P_het04(i,j) = max(krcx_het04(tmp3));%krS5 posterior
        end 
    end
end
for i = 1:100
    tmp1 = find(krM_het05(i,:) >= thr);
    krS5actOT_het05{i} = intersect(tmp1, krstripe5_het05);
    for j = 1:5
        tmp3 = intersect(krcybin_het05{j}, krS5actOT_het05{i});
        
        if length(tmp3) >=2
            tmp4 = sort(krcx_het05(tmp3),'descend');
            if ismember(i,[70 71 72 75]) == 1 && j == 5
                krS5A_het05(i,j) = min(krcx_het05(tmp3));%krS5 anterior
                krS5P_het05(i,j) = tmp4(2);%krS5 posterior
            else
                if abs(tmp4(1) - tmp4(2))>D
                    krS5A_het05(i,j) = min(krcx_het05(tmp3));%S5 anterior
                    krS5P_het05(i,j) = tmp4(2);%S5 posterior
                else
                    krS5A_het05(i,j) = min(krcx_het05(tmp3));%krS5 anterior
                    krS5P_het05(i,j) = max(krcx_het05(tmp3));%krS5 posterior
                end 
            end
        elseif length(tmp3) == 2
                krS5A_het05(i,j) = min(krcx_het05(tmp3));%krS5 anterior
                krS5P_het05(i,j) = max(krcx_het05(tmp3));%krS5 posterior
        end
    end
end

S5A_wt01(S5A_wt01 == 0) = NaN;
S5A_wt02(S5A_wt02 == 0) = NaN;
S5A_wt03(S5A_wt03 == 0) = NaN;
S5P_wt01(S5P_wt01 == 0) = NaN;
S5P_wt02(S5P_wt02 == 0) = NaN;
S5P_wt03(S5P_wt03 == 0) = NaN;
krS5A_wt01(krS5A_wt01 == 0) = NaN;
krS5A_wt02(krS5A_wt02 == 0) = NaN;
krS5A_wt03(krS5A_wt03 == 0) = NaN;
krS5A_wt04(krS5A_wt04 == 0) = NaN;
krS5P_wt01(krS5P_wt01 == 0) = NaN;
krS5P_wt02(krS5P_wt02 == 0) = NaN;
krS5P_wt03(krS5P_wt03 == 0) = NaN;
krS5P_wt04(krS5P_wt04 == 0) = NaN;
krS5A_het01(krS5A_het01 == 0) = NaN;
krS5A_het02(krS5A_het02 == 0) = NaN;
krS5A_het03(krS5A_het03 == 0) = NaN;
krS5A_het04(krS5A_het04 == 0) = NaN;
krS5A_het05(krS5A_het05 == 0) = NaN;
krS5P_het01(krS5P_het01 == 0) = NaN;
krS5P_het02(krS5P_het02 == 0) = NaN;
krS5P_het03(krS5P_het03 == 0) = NaN;
krS5P_het04(krS5P_het04 == 0) = NaN;
krS5P_het05(krS5P_het05 == 0) = NaN;

%% calculate mean distance
%ignore A&P that shares the same nuclei
for i = 1:100
    tmp = find(S5A_wt01(i,:) == S5P_wt01(i,:));
    if length(tmp) > 0 
        S5A_wt01(i,tmp) = NaN;
        S5P_wt01(i,tmp) = NaN;
    end
end
for i = 1:100
    tmp = find(S5A_wt02(i,:) == S5P_wt02(i,:));
    if length(tmp) > 0 
        S5A_wt02(i,tmp) = NaN;
        S5P_wt02(i,tmp) = NaN;
    end
end
for i = 1:100
    tmp = find(S5A_wt03(i,:) == S5P_wt03(i,:));
    if length(tmp) > 0 
        S5A_wt03(i,tmp) = NaN;
        S5P_wt03(i,tmp) = NaN;
    end
end
for i = 1:100
    tmp = find(krS5A_wt01(i,:) == krS5P_wt01(i,:));
    if length(tmp) > 0 
        krS5A_wt01(i,tmp) = NaN;
        krS5P_wt01(i,tmp) = NaN;
    end
end
for i = 1:100
    tmp = find(krS5A_wt02(i,:) == krS5P_wt02(i,:));
    if length(tmp) > 0 
        krS5A_wt02(i,tmp) = NaN;
        krS5P_wt02(i,tmp) = NaN;
    end
end
for i = 1:100
    tmp = find(krS5A_wt03(i,:) == krS5P_wt03(i,:));
    if length(tmp) > 0 
        krS5A_wt03(i,tmp) = NaN;
        krS5P_wt03(i,tmp) = NaN;
    end
end
for i = 1:100
    tmp = find(krS5A_wt04(i,:) == krS5P_wt04(i,:));
    if length(tmp) > 0 
        krS5A_wt04(i,tmp) = NaN;
        krS5P_wt04(i,tmp) = NaN;
    end
end
for i = 1:100
    tmp = find(krS5A_het01(i,:) == krS5P_het01(i,:));
    if length(tmp) > 0 
        krS5A_het01(i,tmp) = NaN;
        krS5P_het01(i,tmp) = NaN;
    end
end
for i = 1:100
    tmp = find(krS5A_het02(i,:) == krS5P_het02(i,:));
    if length(tmp) > 0 
        krS5A_het02(i,tmp) = NaN;
        krS5P_het02(i,tmp) = NaN;
    end
end
for i = 1:100
    tmp = find(krS5A_het03(i,:) == krS5P_het03(i,:));
    if length(tmp) > 0 
        krS5A_het03(i,tmp) = NaN;
        krS5P_het03(i,tmp) = NaN;
    end
end
for i = 1:100
    tmp = find(krS5A_het04(i,:) == krS5P_het04(i,:));
    if length(tmp) > 0 
        krS5A_het04(i,tmp) = NaN;
        krS5P_het04(i,tmp) = NaN;
    end
end
for i = 1:100
    tmp = find(krS5A_het05(i,:) == krS5P_het05(i,:));
    if length(tmp) > 0 
        krS5A_het05(i,tmp) = NaN;
        krS5P_het05(i,tmp) = NaN;
    end
end

%% calculate stripe width when >=3 bins have signals
for i = 1:100
    if length(find(isnan(S5A_wt01(i,:)) == 0))>=3
       sigS5A_wt01(i) = 1;
    else
       sigS5A_wt01(i) = 0;
    end

    if length(find(isnan(S5A_wt02(i,:)) == 0))>=3
       sigS5A_wt02(i) = 1;
       else
       sigS5A_wt02(i) = 0;
    end
    if length(find(isnan(S5A_wt03(i,:)) == 0))>=3
       sigS5A_wt03(i) = 1;
       else
       sigS5A_wt03(i) = 0;
    end

    if length(find(isnan(krS5A_wt01(i,:)) == 0))>=3
       krsigS5A_wt01(i) = 1;
    else
       krsigS5A_wt01(i) = 0;
    end
    if length(find(isnan(krS5A_wt02(i,:)) == 0))>=3
       krsigS5A_wt02(i) = 1;
       else
       krsigS5A_wt02(i) = 0;
    end
    if length(find(isnan(krS5A_wt03(i,:)) == 0))>=3
       krsigS5A_wt03(i) = 1;
       else
       krsigS5A_wt03(i) = 0;
    end
    if length(find(isnan(krS5A_wt04(i,:)) == 0))>=3
       krsigS5A_wt04(i) = 1;
       else
       krsigS5A_wt04(i) = 0;
    end
   
   
   if length(find(isnan(krS5A_het01(i,:)) == 0))>=3
       krsigS5A_het01(i) = 1;
    else
       krsigS5A_het01(i) = 0;
    end
    if length(find(isnan(krS5A_het02(i,:)) == 0))>=3
       krsigS5A_het02(i) = 1;
       else
       krsigS5A_het02(i) = 0;
    end
    if length(find(isnan(krS5A_het03(i,:)) == 0))>=3
       krsigS5A_het03(i) = 1;
       else
       krsigS5A_het03(i) = 0;
    end
    if length(find(isnan(krS5A_het04(i,:)) == 0))>=3
       krsigS5A_het04(i) = 1;
       else
       krsigS5A_het04(i) = 0;
    end
    if length(find(isnan(krS5A_het05(i,:)) == 0))>=3
       krsigS5A_het05(i) = 1;
       else
       krsigS5A_het05(i) = 0;
    end
   
end

%calculate width
for i = 1:100
    if sigS5A_wt01(i) == 1
        S5W_wt01(i) = nanmean(S5P_wt01(i,:) - S5A_wt01(i,:))/length_wt01;
    else
        S5W_wt01(i) = 0;
    end
    if sigS5A_wt02(i) == 1
        S5W_wt02(i) = nanmean(S5P_wt02(i,:) - S5A_wt02(i,:))/length_wt02;
    else
        S5W_wt02(i) = 0;
    end
    if sigS5A_wt03(i) == 1
        S5W_wt03(i) = nanmean(S5P_wt03(i,:) - S5A_wt03(i,:))/length_wt03;
    else
        S5W_wt03(i) = 0;
    end

    if krsigS5A_wt01(i) == 1
        krS5W_wt01(i) = nanmean(krS5P_wt01(i,:) - krS5A_wt01(i,:))/krlength_wt01;
    else
        krS5W_wt01(i) = 0;
    end
    if krsigS5A_wt02(i) == 1
        krS5W_wt02(i) = nanmean(krS5P_wt02(i,:) - krS5A_wt02(i,:))/krlength_wt02;
    else
        krS5W_wt02(i) = 0;
    end
    if krsigS5A_wt03(i) == 1
        krS5W_wt03(i) = nanmean(krS5P_wt03(i,:) - krS5A_wt03(i,:))/krlength_wt03;
    else
        krS5W_wt03(i) = 0;
    end
    if krsigS5A_wt04(i) == 1
        krS5W_wt04(i) = nanmean(krS5P_wt04(i,:) - krS5A_wt04(i,:))/krlength_wt04;
    else
        krS5W_wt04(i) = 0;
    end

    if krsigS5A_het01(i) == 1
        krS5W_het01(i) = nanmean(krS5P_het01(i,:) - krS5A_het01(i,:))/krlength_het01;
    else
        krS5W_het01(i) = 0;
    end
    if krsigS5A_het02(i) == 1
        krS5W_het02(i) = nanmean(krS5P_het02(i,:) - krS5A_het02(i,:))/krlength_het02;
    else
        krS5W_het02(i) = 0;
    end
    if krsigS5A_het03(i) == 1
        krS5W_het03(i) = nanmean(krS5P_het03(i,:) - krS5A_het03(i,:))/krlength_het03;
    else
        krS5W_het03(i) = 0;
    end
    if krsigS5A_het04(i) == 1
        krS5W_het04(i) = nanmean(krS5P_het04(i,:) - krS5A_het04(i,:))/krlength_het04;
    else
        krS5W_het04(i) = 0;
    end
    if krsigS5A_het05(i) == 1
        krS5W_het05(i) = nanmean(krS5P_het05(i,:) - krS5A_het05(i,:))/krlength_het05;
    else
        krS5W_het05(i) = 0;
    end
    
    tmp_wt = [S5W_wt01(i) S5W_wt02(i) S5W_wt03(i) krS5W_wt01(i) krS5W_wt02(i) krS5W_wt03(i) krS5W_wt04(i)];
    tmp_het = [krS5W_het01(i) krS5W_het02(i) krS5W_het03(i) krS5W_het04(i) krS5W_het05(i)];
    S5W_wt(i) = nanmean(tmp_wt);
    S5W_het(i) = nanmean(tmp_het);

    S5Wse_wt(i) = nanstd(tmp_wt)/sqrt(length(find(tmp_wt ~= NaN)));
    S5Wse_het(i) = nanstd(tmp_het)/sqrt(length(find(tmp_het ~= NaN)));
end


%% anterior and posterior border
for i = 1:100
    S5AA_wt01(i) = nanmean(S5A_wt01(i,:)); %average anterior border position
    S5AA_wt02(i) = nanmean(S5A_wt02(i,:));
    S5AA_wt03(i) = nanmean(S5A_wt03(i,:));
    krS5AA_wt01(i) = nanmean(krS5A_wt01(i,:)); 
    krS5AA_wt02(i) = nanmean(krS5A_wt02(i,:));
    krS5AA_wt03(i) = nanmean(krS5A_wt03(i,:));
    krS5AA_wt04(i) = nanmean(krS5A_wt04(i,:));
    krS5AA_het01(i) = nanmean(krS5A_het01(i,:)); 
    krS5AA_het02(i) = nanmean(krS5A_het02(i,:));
    krS5AA_het03(i) = nanmean(krS5A_het03(i,:));
    krS5AA_het04(i) = nanmean(krS5A_het04(i,:));
    krS5AA_het05(i) = nanmean(krS5A_het05(i,:));

    S5PA_wt01(i) = nanmean(S5P_wt01(i,:)); %posterior average border position
    S5PA_wt02(i) = nanmean(S5P_wt02(i,:));
    S5PA_wt03(i) = nanmean(S5P_wt03(i,:));
    krS5PA_wt01(i) = nanmean(krS5P_wt01(i,:)); 
    krS5PA_wt02(i) = nanmean(krS5P_wt02(i,:));
    krS5PA_wt03(i) = nanmean(krS5P_wt03(i,:));
    krS5PA_wt04(i) = nanmean(krS5P_wt04(i,:));
    krS5PA_het01(i) = nanmean(krS5P_het01(i,:)); 
    krS5PA_het02(i) = nanmean(krS5P_het02(i,:));
    krS5PA_het03(i) = nanmean(krS5P_het03(i,:));
    krS5PA_het04(i) = nanmean(krS5P_het04(i,:));
    krS5PA_het05(i) = nanmean(krS5P_het05(i,:));
end

%calibrate to %EL
S5AAEL_wt01 = (S5AA_wt01 - min(tip_wt01))/(max(tip_wt01)- min(tip_wt01));
S5AAEL_wt02 = (S5AA_wt02 - min(tip_wt02))/(max(tip_wt02)- min(tip_wt02));
S5AAEL_wt03 = (S5AA_wt03 - min(tip_wt03))/(max(tip_wt03)- min(tip_wt03));
krS5AAEL_wt01 = (krS5AA_wt01 - min(krtip_wt01))/(max(krtip_wt01)- min(krtip_wt01));
krS5AAEL_wt02 = (krS5AA_wt02 - min(krtip_wt02))/(max(krtip_wt02)- min(krtip_wt02));
krS5AAEL_wt03 = (krS5AA_wt03 - min(krtip_wt03))/(max(krtip_wt03)- min(krtip_wt03));
krS5AAEL_wt04 = (krS5AA_wt04 - min(krtip_wt04))/(max(krtip_wt04)- min(krtip_wt04));
krS5AAEL_het01 = (krS5AA_het01 - min(krtip_het01))/(max(krtip_het01)- min(krtip_het01));
krS5AAEL_het02 = (krS5AA_het02 - min(krtip_het02))/(max(krtip_het02)- min(krtip_het02));
krS5AAEL_het03 = (krS5AA_het03 - min(krtip_het03))/(max(krtip_het03)- min(krtip_het03));
krS5AAEL_het04 = (krS5AA_het04 - min(krtip_het04))/(max(krtip_het04)- min(krtip_het04));
krS5AAEL_het05 = (krS5AA_het05 - min(krtip_het05))/(max(krtip_het05)- min(krtip_het05));

S5PAEL_wt01 = (S5PA_wt01 - min(tip_wt01))/(max(tip_wt01)- min(tip_wt01));
S5PAEL_wt02 = (S5PA_wt02 - min(tip_wt02))/(max(tip_wt02)- min(tip_wt02));
S5PAEL_wt03 = (S5PA_wt03 - min(tip_wt03))/(max(tip_wt03)- min(tip_wt03));
krS5PAEL_wt01 = (krS5PA_wt01 - min(krtip_wt01))/(max(krtip_wt01)- min(krtip_wt01));
krS5PAEL_wt02 = (krS5PA_wt02 - min(krtip_wt02))/(max(krtip_wt02)- min(krtip_wt02));
krS5PAEL_wt03 = (krS5PA_wt03 - min(krtip_wt03))/(max(krtip_wt03)- min(krtip_wt03));
krS5PAEL_wt04 = (krS5PA_wt04 - min(krtip_wt04))/(max(krtip_wt04)- min(krtip_wt04));
krS5PAEL_het01 = (krS5PA_het01 - min(krtip_het01))/(max(krtip_het01)- min(krtip_het01));
krS5PAEL_het02 = (krS5PA_het02 - min(krtip_het02))/(max(krtip_het02)- min(krtip_het02));
krS5PAEL_het03 = (krS5PA_het03 - min(krtip_het03))/(max(krtip_het03)- min(krtip_het03));
krS5PAEL_het04 = (krS5PA_het04 - min(krtip_het04))/(max(krtip_het04)- min(krtip_het04));
krS5PAEL_het05 = (krS5PA_het05 - min(krtip_het05))/(max(krtip_het05)- min(krtip_het05));

%average S5 border position
for i = 1:100
    tmp_wt = [S5AAEL_wt01(i) S5AAEL_wt02(i) S5AAEL_wt03(i) krS5AAEL_wt01(i) krS5AAEL_wt02(i) krS5AAEL_wt03(i) krS5AAEL_wt04(i)];
    tmp_het = [krS5AAEL_het01(i) krS5AAEL_het02(i) krS5AAEL_het03(i) krS5AAEL_het04(i) krS5AAEL_het05(i)];

    S5AAEL_wt(i) = nanmean(tmp_wt);
    S5AAEL_het(i) = nanmean(tmp_het);
    S5AAELse_wt(i) = nanstd(tmp_wt)/sqrt(length(find(tmp_wt ~= NaN)));
    S5AAELse_het(i) = nanstd(tmp_het)/sqrt(length(find(tmp_het ~= NaN)));
end

for i = 1:100
    tmp_wt = [S5PAEL_wt01(i) S5PAEL_wt02(i) S5PAEL_wt03(i) krS5PAEL_wt01(i) krS5PAEL_wt02(i) krS5PAEL_wt03(i) krS5PAEL_wt04(i)];
    tmp_het = [krS5PAEL_het01(i) krS5PAEL_het02(i) krS5PAEL_het03(i) krS5PAEL_het04(i) krS5PAEL_het05(i)];

    S5PAEL_wt(i) = nanmean(tmp_wt);
    S5PAEL_het(i) = nanmean(tmp_het);
    S5PAELse_wt(i) = nanstd(tmp_wt)/sqrt(length(find(tmp_wt ~= NaN)));
    S5PAELse_het(i) = nanstd(tmp_het)/sqrt(length(find(tmp_het ~= NaN)));
end

S5AAELseL_wt = [zeros(1,length(S5AAELse_wt));S5AAELse_wt];
S5PAELseP_wt = [S5PAELse_wt;zeros(1,length(S5PAELse_wt))];
S5AAELseL_het = [zeros(1,length(S5AAELse_het));S5AAELse_het];
S5PAELseP_het = [S5PAELse_het;zeros(1,length(S5PAELse_het))];

cNaN_wt = length(find(isnan(S5AAEL_wt)==1));
cNaN_het = length(find(isnan(S5AAEL_het)==1));

%% plot 2nd half NC14
figure(232);
plot(0.5:0.01:1, S5W_wt(50:100)*100,'LineWidth',3,'Color','k');hold on
plot(0.5:0.01:1, S5W_het(50:100)*100,'LineWidth',3,'Color','r');
shadedErrorBar2(0.50:0.01:1, S5W_wt(50:100)*100, S5Wse_wt(50:100)*100,'lineprops', {'-k','LineWidth',3}); hold on
shadedErrorBar2(0.50:0.01:1, S5W_het(50:100)*100, S5Wse_het(50:100)*100,'lineprops', {'-r','LineWidth',3}); 
legend({'wt', 'Kr^{1}/+'},'Location','northeast');
xlim([0.5 1]);
ylim([0 5]);

figure(238);
patch([0.50:0.01:1 fliplr(0.50:0.01:1)],[S5AAEL_wt(50:100) fliplr(S5PAEL_wt(50:100))]*100,'k','FaceAlpha',.5,'EdgeColor','none');hold on
patch([0.50:0.01:1 fliplr(0.50:0.01:1)],[S5AAEL_het(50:100) fliplr(S5PAEL_het(50:100))]*100,'r','FaceAlpha',.5,'EdgeColor','none');
shadedErrorBar2(0.50:0.01:1, S5AAEL_wt(50:100)*100, S5AAELseL_wt(:,50:100)*100,'lineprops', {'-k','LineWidth',1}); hold on
shadedErrorBar2(0.50:0.01:1, S5AAEL_het(50:100)*100, S5AAELseL_het(:,50:100)*100,'lineprops', {'-r','LineWidth',1}); 
shadedErrorBar2(0.50:0.01:1, S5PAEL_wt(50:100)*100, S5PAELseP_wt(:,50:100)*100,'lineprops', {'-k','LineWidth',1});
shadedErrorBar2(0.50:0.01:1, S5PAEL_het(50:100)*100, S5PAELseP_het(:,50:100)*100,'lineprops', {'-r','LineWidth',1});
legend({'wt', 'Kr^{1}/+'},'Location','southeast');
% title('Stripe 2 location');
set(gca,'View',[90 90]);
% xlabel('nc 14 (normalized)');
% ylabel('anterior \rightarrow posterior (%EL)');
xlim([0.5 1]);
ylim([54 67]);
%% load S3 data and compare S3-5 width
load('S3_boundary_tracing.mat');

for i = 1:100
    S35W_wt01(i) = nanmean(S5P_wt01(i,:) - S3A_wt01(i,:))/length_wt01;
    S35W_wt02(i) = nanmean(S5P_wt02(i,:) - S3A_wt02(i,:))/length_wt02;
    S35W_wt03(i) = nanmean(S5P_wt03(i,:) - S3A_wt03(i,:))/length_wt03;

    krS35W_wt01(i) = nanmean(krS5P_wt01(i,:) - krS3A_wt01(i,:))/krlength_wt01;
    krS35W_wt02(i) = nanmean(krS5P_wt02(i,:) - krS3A_wt02(i,:))/krlength_wt02;
    krS35W_wt03(i) = nanmean(krS5P_wt03(i,:) - krS3A_wt03(i,:))/krlength_wt03;
    krS35W_wt04(i) = nanmean(krS5P_wt04(i,:) - krS3A_wt04(i,:))/krlength_wt04;

    krS35W_het01(i) = nanmean(krS5P_het01(i,:) - krS3A_het01(i,:))/krlength_het01;
    krS35W_het02(i) = nanmean(krS5P_het02(i,:) - krS3A_het02(i,:))/krlength_het02;
    krS35W_het03(i) = nanmean(krS5P_het03(i,:) - krS3A_het03(i,:))/krlength_het03;
    krS35W_het04(i) = nanmean(krS5P_het04(i,:) - krS3A_het04(i,:))/krlength_het04;
    krS35W_het05(i) = nanmean(krS5P_het05(i,:) - krS3A_het05(i,:))/krlength_het05;
    
    tmp_wt = [S35W_wt01(i) S35W_wt02(i) S35W_wt03(i) krS35W_wt01(i) krS35W_wt02(i) krS35W_wt03(i) krS35W_wt04(i)];
    tmp_het = [krS35W_het01(i) krS35W_het02(i) krS35W_het03(i) krS35W_het04(i) krS35W_het05(i)];
    S35W_wt(i) = nanmean(tmp_wt);
    S35W_het(i) = nanmean(tmp_het);

    S35Wse_wt(i) = nanstd(tmp_wt)/sqrt(length(find(tmp_wt ~= NaN)));
    S35Wse_het(i) = nanstd(tmp_het)/sqrt(length(find(tmp_het ~= NaN)));
end

figure(312);
plot(0.5:0.01:1, S35W_wt(50:100)*100,'LineWidth',3,'Color','k');hold on
plot(0.5:0.01:1, S35W_het(50:100)*100,'LineWidth',3,'Color','r');
shadedErrorBar2(0.5:0.01:1, S35W_wt(50:100)*100, S35Wse_wt(50:100)*100,'lineprops', {'-k','LineWidth',2}); hold on
shadedErrorBar2(0.5:0.01:1, S35W_het(50:100)*100, S35Wse_het(50:100)*100,'lineprops', {'-r','LineWidth',2}); 
legend({'{\itwt}', ['\itKr^{1}/+']},'Location','northeast');
xlim([0.5 1]);
ax = gca;
ax.YAxis.FontSize = 24;
ax.XAxis.FontSize = 24;
ax.Legend.FontSize = 24;
ax.Title.FontSize = 24;
hold off

%% measure position difference of stripe - getting the middle line between anterior and posterior boundary
for i = 1:100
    if S5AAEL_wt(i) ~= NaN
        S5EL_wt(i) = S5AAEL_wt(i) + 0.5*(S5PAEL_wt(i) - S5AAEL_wt(i));
    else
        S5EL_wt(i) = NaN;
    end

    if S5AAEL_het(i) ~= NaN
        S5EL_het(i) = S5AAEL_het(i) + 0.5*(S5PAEL_het(i) - S5AAEL_het(i));
    else
        S5EL_wt(i) = NaN;
    end
end

for i = 1:100
    S5EL_D(i) = (S5EL_wt(i) - S5EL_het(i))*100;
end

mean(S5EL_D(80:90))
