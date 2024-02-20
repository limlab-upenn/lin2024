clear 
close all

%% load and set
load('D:\LIM LAB\kr_project_code_cleanup\dataset_whole_hb_NC14_100f_fixedY-150.mat');
load('D:\LIM LAB\kr_project_code_cleanup\dataset_whole_kr_NC14_100f_fixedY-150.mat');
load('D:\LIM LAB\kr_project_code_cleanup\Stripe_three_to_six_anterior_borders.mat');
load('D:\LIM LAB\kr_project_code_cleanup\Stripe_one_two_seven_anterior_borders.mat');
load('myCustomColormap.mat');

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
  maxM3_wt01(i) = max(M_wt01(round(0.51*size(M_wt01,1)):end,i));%max intensity for last 1/2 frames
end 
for i=1:size(M_wt02,2) % all nuclei
  maxM_wt02(i) = max(M_wt02(:,i)); % obtain max intensity for each nuclues
  maxM2_wt02(i) = max(M_wt02(round(0.75*size(M_wt02,1)):end,i));%max intensity for last 1/4 frames
  maxM3_wt02(i) = max(M_wt02(round(0.51*size(M_wt02,1)):end,i));%max intensity for last 1/2 frames
end 
for i=1:size(M_wt03,2) % all nuclei
  maxM_wt03(i) = max(M_wt03(:,i)); % obtain max intensity for each nuclues
  maxM2_wt03(i) = max(M_wt03(round(0.75*size(M_wt03,1)):end,i));%max intensity for last 1/4 frames
  maxM3_wt03(i) = max(M_wt03(round(0.51*size(M_wt03,1)):end,i));%max intensity for last 1/2 frames
end 
for i=1:size(krM_wt01,2) % all nuclei
  krmaxM_wt01(i) = max(krM_wt01(:,i)); % obtain max intensity for each nuclues
  krmaxM2_wt01(i) = max(krM_wt01(round(0.75*size(krM_wt01,1)):end,i));%max intensity for last 1/4 frames
  krmaxM3_wt01(i) = max(krM_wt01(round(0.51*size(krM_wt01,1)):end,i));%max intensity for last 1/2 frames
end 
for i=1:size(krM_wt02,2) % all nuclei
  krmaxM_wt02(i) = max(krM_wt02(:,i)); % obtain max intensity for each nuclues
  krmaxM2_wt02(i) = max(krM_wt02(round(0.75*size(krM_wt02,1)):end,i));%max intensity for last 1/4 frames
  krmaxM3_wt02(i) = max(krM_wt02(round(0.51*size(krM_wt02,1)):end,i));%max intensity for last 1/2 frames
end 
for i=1:size(krM_wt03,2) % all nuclei
  krmaxM_wt03(i) = max(krM_wt03(:,i)); % obtain max intensity for each nuclues
  krmaxM2_wt03(i) = max(krM_wt03(round(0.75*size(krM_wt03,1)):end,i));%max intensity for last 1/4 frames
  krmaxM3_wt03(i) = max(krM_wt03(round(0.51*size(krM_wt03,1)):end,i));%max intensity for last 1/2 frames
end 
for i=1:size(krM_wt04,2) % all nuclei
  krmaxM_wt04(i) = max(krM_wt04(:,i)); % obtain max intensity for each nuclues
  krmaxM2_wt04(i) = max(krM_wt04(round(0.75*size(krM_wt04,1)):end,i));%max intensity for last 1/4 frames
  krmaxM3_wt04(i) = max(krM_wt04(round(0.51*size(krM_wt04,1)):end,i));%max intensity for last 1/2 frames
end 

%het
for i=1:size(krM_het01,2) % all nuclei
  krmaxM_het01(i) = max(krM_het01(:,i)); % obtain max intensity for each nuclues
  krmaxM2_het01(i) = max(krM_het01(round(0.75*size(krM_het01,1)):end,i));%max intensity for last 1/4 frames
  krmaxM3_het01(i) = max(krM_het01(round(0.51*size(krM_het01,1)):end,i));%max intensity for last 1/2 frames
end 
for i=1:size(krM_het02,2) % all nuclei
  krmaxM_het02(i) = max(krM_het02(:,i)); % obtain max intensity for each nuclues
  krmaxM2_het02(i) = max(krM_het02(round(0.75*size(krM_het02,1)):end,i));%max intensity for last 1/4 frames
  krmaxM3_het02(i) = max(krM_het02(round(0.51*size(krM_het02,1)):end,i));%max intensity for last 1/2 frames
end 
for i=1:size(krM_het03,2) % all nuclei
  krmaxM_het03(i) = max(krM_het03(:,i)); % obtain max intensity for each nuclues
  krmaxM2_het03(i) = max(krM_het03(round(0.75*size(krM_het03,1)):end,i));%max intensity for last 1/4 frames
  krmaxM3_het03(i) = max(krM_het03(round(0.51*size(krM_het03,1)):end,i));%max intensity for last 1/2 frames
end 
for i=1:size(krM_het04,2) % all nuclei
  krmaxM_het04(i) = max(krM_het04(:,i)); % obtain max intensity for each nuclues
  krmaxM2_het04(i) = max(krM_het04(round(0.75*size(krM_het04,1)):end,i));%max intensity for last 1/4 frames
  krmaxM3_het04(i) = max(krM_het04(round(0.51*size(krM_het04,1)):end,i));%max intensity for last 1/2 frames
end 
for i=1:size(krM_het05,2) % all nuclei
  krmaxM_het05(i) = max(krM_het05(:,i)); % obtain max intensity for each nuclues
  krmaxM2_het05(i) = max(krM_het05(round(0.75*size(krM_het05,1)):end,i));%max intensity for last 1/4 frames
  krmaxM3_het05(i) = max(krM_het05(round(0.51*size(krM_het05,1)):end,i));%max intensity for last 1/2 frames
end 

%% active nuclei
act_wt01 = find(maxM_wt01 > thr);
actLate_wt01 = find(maxM2_wt01 > thr);
act_wt02 = find(maxM_wt02 > thr);
actLate_wt02 = find(maxM2_wt02 > thr);
act_wt03 = find(maxM_wt03 > thr);
actLate_wt03 = find(maxM2_wt03 > thr);

kract_wt01 = find(krmaxM_wt01 > thr);
kractLate_wt01 = find(krmaxM2_wt01 > thr);
kract_wt02 = find(krmaxM_wt02 > thr);
kractLate_wt02 = find(krmaxM2_wt02 > thr);
kract_wt03 = find(krmaxM_wt03 > thr);
kractLate_wt03 = find(krmaxM2_wt03 > thr);
kract_wt04 = find(krmaxM_wt04 > thr);
kractLate_wt04 = find(krmaxM2_wt04 > thr);

kract_het01 = find(krmaxM_het01 > thr);
kractLate_het01 = find(krmaxM2_het01 > thr);
kract_het02 = find(krmaxM_het02 > thr);
kractLate_het02 = find(krmaxM2_het02 > thr);
kract_het03 = find(krmaxM_het03 > thr);
kractLate_het03 = find(krmaxM2_het03 > thr);
kract_het04 = find(krmaxM_het04 > thr);
kractLate_het04 = find(krmaxM2_het04 > thr);
kract_het05 = find(krmaxM_het05 > thr);
kractLate_het05 = find(krmaxM2_het05 > thr);

act2nd_wt01 = find(maxM3_wt01 > thr);
act2nd_wt02 = find(maxM3_wt02 > thr);
act2nd_wt03 = find(maxM3_wt03 > thr);

kract2nd_wt01 = find(krmaxM3_wt01 > thr);
kract2nd_wt02 = find(krmaxM3_wt02 > thr);
kract2nd_wt03 = find(krmaxM3_wt03 > thr);
kract2nd_wt04 = find(krmaxM3_wt04 > thr);

kract2nd_het01 = find(krmaxM3_het01 > thr);
kract2nd_het02 = find(krmaxM3_het02 > thr);
kract2nd_het03 = find(krmaxM3_het03 > thr);
kract2nd_het04 = find(krmaxM3_het04 > thr);
kract2nd_het05 = find(krmaxM3_het05 > thr);

%wt01
%S1
S1Polygonx_wt01 = cat(1,S1Anterior_wt01(:,1), flip(S2Anterior_wt01(:,1)))';
S1Polygony_wt01 = cat(1,S1Anterior_wt01(:,2), flip(S2Anterior_wt01(:,2)))';
S1Pgon_wt01 = polyshape(S1Polygonx_wt01,S1Polygony_wt01);
stripe1_wt01 = find(inpolygon(cx_wt01,cy_wt01,S1Polygonx_wt01,S1Polygony_wt01)== 1);
stripe1Act_wt01 = intersect(stripe1_wt01, act_wt01);
%S2
S2Polygonx_wt01 = cat(1,S2Anterior_wt01(:,1), flip(S3Anterior_wt01(:,1)))';
S2Polygony_wt01 = cat(1,S2Anterior_wt01(:,2), flip(S3Anterior_wt01(:,2)))';
S2Pgon_wt01 = polyshape(S2Polygonx_wt01,S2Polygony_wt01);
stripe2_wt01 = find(inpolygon(cx_wt01,cy_wt01,S2Polygonx_wt01,S2Polygony_wt01)== 1);
stripe2Act_wt01 = intersect(stripe2_wt01, act_wt01);
%S6
S6Polygonx_wt01 = cat(1,S6Anterior_wt01(:,1), flip(S7Anterior_wt01(:,1)))';
S6Polygony_wt01 = cat(1,S6Anterior_wt01(:,2), flip(S7Anterior_wt01(:,2)))';
S6Pgon_wt01 = polyshape(S6Polygonx_wt01,S6Polygony_wt01);
stripe6_wt01 = find(inpolygon(cx_wt01,cy_wt01,S6Polygonx_wt01,S6Polygony_wt01)== 1);
stripe6Act_wt01 = intersect(stripe6_wt01, act_wt01);

%wt02
%S1
S1Polygonx_wt02 = cat(1,S1Anterior_wt02(:,1), flip(S2Anterior_wt02(:,1)))';
S1Polygony_wt02 = cat(1,S1Anterior_wt02(:,2), flip(S2Anterior_wt02(:,2)))';
S1Pgon_wt02 = polyshape(S1Polygonx_wt02,S1Polygony_wt02);
stripe1_wt02 = find(inpolygon(cx_wt02,cy_wt02,S1Polygonx_wt02,S1Polygony_wt02)== 1);
stripe1Act_wt02 = intersect(stripe1_wt02, act_wt02);
%S2
S2Polygonx_wt02 = cat(1,S2Anterior_wt02(:,1), flip(S3Anterior_wt02(:,1)))';
S2Polygony_wt02 = cat(1,S2Anterior_wt02(:,2), flip(S3Anterior_wt02(:,2)))';
S2Pgon_wt02 = polyshape(S2Polygonx_wt02,S2Polygony_wt02);
stripe2_wt02 = find(inpolygon(cx_wt02,cy_wt02,S2Polygonx_wt02,S2Polygony_wt02)== 1);
stripe2Act_wt02 = intersect(stripe2_wt02, act_wt02);
%S6
S6Polygonx_wt02 = cat(1,S6Anterior_wt02(:,1), flip(S7Anterior_wt02(:,1)))';
S6Polygony_wt02 = cat(1,S6Anterior_wt02(:,2), flip(S7Anterior_wt02(:,2)))';
S6Pgon_wt02 = polyshape(S6Polygonx_wt02,S6Polygony_wt02);
stripe6_wt02 = find(inpolygon(cx_wt02,cy_wt02,S6Polygonx_wt02,S6Polygony_wt02)== 1);
stripe6Act_wt02 = intersect(stripe6_wt02, act_wt02);

%wt03
%S1
S1Polygonx_wt03 = cat(1,S1Anterior_wt03(:,1), flip(S2Anterior_wt03(:,1)))';
S1Polygony_wt03 = cat(1,S1Anterior_wt03(:,2), flip(S2Anterior_wt03(:,2)))';
S1Pgon_wt03 = polyshape(S1Polygonx_wt03,S1Polygony_wt03);
stripe1_wt03 = find(inpolygon(cx_wt03,cy_wt03,S1Polygonx_wt03,S1Polygony_wt03)== 1);
stripe1Act_wt03 = intersect(stripe1_wt03, act_wt03);
%S2
S2Polygonx_wt03 = cat(1,S2Anterior_wt03(:,1), flip(S3Anterior_wt03(:,1)))';
S2Polygony_wt03 = cat(1,S2Anterior_wt03(:,2), flip(S3Anterior_wt03(:,2)))';
S2Pgon_wt03 = polyshape(S2Polygonx_wt03,S2Polygony_wt03);
stripe2_wt03 = find(inpolygon(cx_wt03,cy_wt03,S2Polygonx_wt03,S2Polygony_wt03)== 1);
stripe2Act_wt03 = intersect(stripe2_wt03, act_wt03);
%S6
S6Polygonx_wt03 = cat(1,S6Anterior_wt03(:,1), flip(S7Anterior_wt03(:,1)))';
S6Polygony_wt03 = cat(1,S6Anterior_wt03(:,2), flip(S7Anterior_wt03(:,2)))';
S6Pgon_wt03 = polyshape(S6Polygonx_wt03,S6Polygony_wt03);
stripe6_wt03 = find(inpolygon(cx_wt03,cy_wt03,S6Polygonx_wt03,S6Polygony_wt03)== 1);
stripe6Act_wt03 = intersect(stripe6_wt03, act_wt03);

%krwt01
%krS1
krS1Polygonx_wt01 = cat(1,krS1Anterior_wt01(:,1), flip(krS2Anterior_wt01(:,1)))';
krS1Polygony_wt01 = cat(1,krS1Anterior_wt01(:,2), flip(krS2Anterior_wt01(:,2)))';
krS1Pgon_wt01 = polyshape(krS1Polygonx_wt01,krS1Polygony_wt01);
krstripe1_wt01 = find(inpolygon(krcx_wt01,krcy_wt01,krS1Polygonx_wt01,krS1Polygony_wt01)== 1);
krstripe1Act_wt01 = intersect(krstripe1_wt01, kract_wt01);
%krS2
krS2Polygonx_wt01 = cat(1,krS2Anterior_wt01(:,1), flip(krS3Anterior_wt01(:,1)))';
krS2Polygony_wt01 = cat(1,krS2Anterior_wt01(:,2), flip(krS3Anterior_wt01(:,2)))';
krS2Pgon_wt01 = polyshape(krS2Polygonx_wt01,krS2Polygony_wt01);
krstripe2_wt01 = find(inpolygon(krcx_wt01,krcy_wt01,krS2Polygonx_wt01,krS2Polygony_wt01)== 1);
krstripe2Act_wt01 = intersect(krstripe2_wt01, kract_wt01);
%krS6
krS6Polygonx_wt01 = cat(1,krS6Anterior_wt01(:,1), flip(krS7Anterior_wt01(:,1)))';
krS6Polygony_wt01 = cat(1,krS6Anterior_wt01(:,2), flip(krS7Anterior_wt01(:,2)))';
krS6Pgon_wt01 = polyshape(krS6Polygonx_wt01,krS6Polygony_wt01);
krstripe6_wt01 = find(inpolygon(krcx_wt01,krcy_wt01,krS6Polygonx_wt01,krS6Polygony_wt01)== 1);
krstripe6Act_wt01 = intersect(krstripe6_wt01, kract_wt01);

%krwt02
%krS1
krS1Polygonx_wt02 = cat(1,krS1Anterior_wt02(:,1), flip(krS2Anterior_wt02(:,1)))';
krS1Polygony_wt02 = cat(1,krS1Anterior_wt02(:,2), flip(krS2Anterior_wt02(:,2)))';
krS1Pgon_wt02 = polyshape(krS1Polygonx_wt02,krS1Polygony_wt02);
krstripe1_wt02 = find(inpolygon(krcx_wt02,krcy_wt02,krS1Polygonx_wt02,krS1Polygony_wt02)== 1);
krstripe1Act_wt02 = intersect(krstripe1_wt02, kract_wt02);
%krS2
krS2Polygonx_wt02 = cat(1,krS2Anterior_wt02(:,1), flip(krS3Anterior_wt02(:,1)))';
krS2Polygony_wt02 = cat(1,krS2Anterior_wt02(:,2), flip(krS3Anterior_wt02(:,2)))';
krS2Pgon_wt02 = polyshape(krS2Polygonx_wt02,krS2Polygony_wt02);
krstripe2_wt02 = find(inpolygon(krcx_wt02,krcy_wt02,krS2Polygonx_wt02,krS2Polygony_wt02)== 1);
krstripe2Act_wt02 = intersect(krstripe2_wt02, kract_wt02);
%krS6
krS6Polygonx_wt02 = cat(1,krS6Anterior_wt02(:,1), flip(krS7Anterior_wt02(:,1)))';
krS6Polygony_wt02 = cat(1,krS6Anterior_wt02(:,2), flip(krS7Anterior_wt02(:,2)))';
krS6Pgon_wt02 = polyshape(krS6Polygonx_wt02,krS6Polygony_wt02);
krstripe6_wt02 = find(inpolygon(krcx_wt02,krcy_wt02,krS6Polygonx_wt02,krS6Polygony_wt02)== 1);
krstripe6Act_wt02 = intersect(krstripe6_wt02, kract_wt02);

%krwt03
%krS1
krS1Polygonx_wt03 = cat(1,krS1Anterior_wt03(:,1), flip(krS2Anterior_wt03(:,1)))';
krS1Polygony_wt03 = cat(1,krS1Anterior_wt03(:,2), flip(krS2Anterior_wt03(:,2)))';
krS1Pgon_wt03 = polyshape(krS1Polygonx_wt03,krS1Polygony_wt03);
krstripe1_wt03 = find(inpolygon(krcx_wt03,krcy_wt03,krS1Polygonx_wt03,krS1Polygony_wt03)== 1);
krstripe1Act_wt03 = intersect(krstripe1_wt03, kract_wt03);
%krS2
krS2Polygonx_wt03 = cat(1,krS2Anterior_wt03(:,1), flip(krS3Anterior_wt03(:,1)))';
krS2Polygony_wt03 = cat(1,krS2Anterior_wt03(:,2), flip(krS3Anterior_wt03(:,2)))';
krS2Pgon_wt03 = polyshape(krS2Polygonx_wt03,krS2Polygony_wt03);
krstripe2_wt03 = find(inpolygon(krcx_wt03,krcy_wt03,krS2Polygonx_wt03,krS2Polygony_wt03)== 1);
krstripe2Act_wt03 = intersect(krstripe2_wt03, kract_wt03);
%krS6
krS6Polygonx_wt03 = cat(1,krS6Anterior_wt03(:,1), flip(krS7Anterior_wt03(:,1)))';
krS6Polygony_wt03 = cat(1,krS6Anterior_wt03(:,2), flip(krS7Anterior_wt03(:,2)))';
krS6Pgon_wt03 = polyshape(krS6Polygonx_wt03,krS6Polygony_wt03);
krstripe6_wt03 = find(inpolygon(krcx_wt03,krcy_wt03,krS6Polygonx_wt03,krS6Polygony_wt03)== 1);
krstripe6Act_wt03 = intersect(krstripe6_wt03, kract_wt03);

%krwt04
%krS1
krS1Polygonx_wt04 = cat(1,krS1Anterior_wt04(:,1), flip(krS2Anterior_wt04(:,1)))';
krS1Polygony_wt04 = cat(1,krS1Anterior_wt04(:,2), flip(krS2Anterior_wt04(:,2)))';
krS1Pgon_wt04 = polyshape(krS1Polygonx_wt04,krS1Polygony_wt04);
krstripe1_wt04 = find(inpolygon(krcx_wt04,krcy_wt04,krS1Polygonx_wt04,krS1Polygony_wt04)== 1);
krstripe1Act_wt04 = intersect(krstripe1_wt04, kract_wt04);
%krS2
krS2Polygonx_wt04 = cat(1,krS2Anterior_wt04(:,1), flip(krS3Anterior_wt04(:,1)))';
krS2Polygony_wt04 = cat(1,krS2Anterior_wt04(:,2), flip(krS3Anterior_wt04(:,2)))';
krS2Pgon_wt04 = polyshape(krS2Polygonx_wt04,krS2Polygony_wt04);
krstripe2_wt04 = find(inpolygon(krcx_wt04,krcy_wt04,krS2Polygonx_wt04,krS2Polygony_wt04)== 1);
krstripe2Act_wt04 = intersect(krstripe2_wt04, kract_wt04);
%krS6
krS6Polygonx_wt04 = cat(1,krS6Anterior_wt04(:,1), flip(krS7Anterior_wt04(:,1)))';
krS6Polygony_wt04 = cat(1,krS6Anterior_wt04(:,2), flip(krS7Anterior_wt04(:,2)))';
krS6Pgon_wt04 = polyshape(krS6Polygonx_wt04,krS6Polygony_wt04);
krstripe6_wt04 = find(inpolygon(krcx_wt04,krcy_wt04,krS6Polygonx_wt04,krS6Polygony_wt04)== 1);
krstripe6Act_wt04 = intersect(krstripe6_wt04, kract_wt04);

%krhet01
%krS1
krS1Polygonx_het01 = cat(1,krS1Anterior_het01(:,1), flip(krS2Anterior_het01(:,1)))';
krS1Polygony_het01 = cat(1,krS1Anterior_het01(:,2), flip(krS2Anterior_het01(:,2)))';
krS1Pgon_het01 = polyshape(krS1Polygonx_het01,krS1Polygony_het01);
krstripe1_het01 = find(inpolygon(krcx_het01,krcy_het01,krS1Polygonx_het01,krS1Polygony_het01)== 1);
krstripe1Act_het01 = intersect(krstripe1_het01, kract_het01);
%krS2
krS2Polygonx_het01 = cat(1,krS2Anterior_het01(:,1), flip(krS3Anterior_het01(:,1)))';
krS2Polygony_het01 = cat(1,krS2Anterior_het01(:,2), flip(krS3Anterior_het01(:,2)))';
krS2Pgon_het01 = polyshape(krS2Polygonx_het01,krS2Polygony_het01);
krstripe2_het01 = find(inpolygon(krcx_het01,krcy_het01,krS2Polygonx_het01,krS2Polygony_het01)== 1);
krstripe2Act_het01 = intersect(krstripe2_het01, kract_het01);
%krS6
krS6Polygonx_het01 = cat(1,krS6Anterior_het01(:,1), flip(krS7Anterior_het01(:,1)))';
krS6Polygony_het01 = cat(1,krS6Anterior_het01(:,2), flip(krS7Anterior_het01(:,2)))';
krS6Pgon_het01 = polyshape(krS6Polygonx_het01,krS6Polygony_het01);
krstripe6_het01 = find(inpolygon(krcx_het01,krcy_het01,krS6Polygonx_het01,krS6Polygony_het01)== 1);
krstripe6Act_het01 = intersect(krstripe6_het01, kract_het01);

%krhet02
%krS1
krS1Polygonx_het02 = cat(1,krS1Anterior_het02(:,1), flip(krS2Anterior_het02(:,1)))';
krS1Polygony_het02 = cat(1,krS1Anterior_het02(:,2), flip(krS2Anterior_het02(:,2)))';
krS1Pgon_het02 = polyshape(krS1Polygonx_het02,krS1Polygony_het02);
krstripe1_het02 = find(inpolygon(krcx_het02,krcy_het02,krS1Polygonx_het02,krS1Polygony_het02)== 1);
krstripe1Act_het02 = intersect(krstripe1_het02, kract_het02);
%krS2
krS2Polygonx_het02 = cat(1,krS2Anterior_het02(:,1), flip(krS3Anterior_het02(:,1)))';
krS2Polygony_het02 = cat(1,krS2Anterior_het02(:,2), flip(krS3Anterior_het02(:,2)))';
krS2Pgon_het02 = polyshape(krS2Polygonx_het02,krS2Polygony_het02);
krstripe2_het02 = find(inpolygon(krcx_het02,krcy_het02,krS2Polygonx_het02,krS2Polygony_het02)== 1);
krstripe2Act_het02 = intersect(krstripe2_het02, kract_het02);
%krS6
krS6Polygonx_het02 = cat(1,krS6Anterior_het02(:,1), flip(krS7Anterior_het02(:,1)))';
krS6Polygony_het02 = cat(1,krS6Anterior_het02(:,2), flip(krS7Anterior_het02(:,2)))';
krS6Pgon_het02 = polyshape(krS6Polygonx_het02,krS6Polygony_het02);
krstripe6_het02 = find(inpolygon(krcx_het02,krcy_het02,krS6Polygonx_het02,krS6Polygony_het02)== 1);
krstripe6Act_het02 = intersect(krstripe6_het02, kract_het02);

%krhet03
%krS1
krS1Polygonx_het03 = cat(1,krS1Anterior_het03(:,1), flip(krS2Anterior_het03(:,1)))';
krS1Polygony_het03 = cat(1,krS1Anterior_het03(:,2), flip(krS2Anterior_het03(:,2)))';
krS1Pgon_het03 = polyshape(krS1Polygonx_het03,krS1Polygony_het03);
krstripe1_het03 = find(inpolygon(krcx_het03,krcy_het03,krS1Polygonx_het03,krS1Polygony_het03)== 1);
krstripe1Act_het03 = intersect(krstripe1_het03, kract_het03);
%krS2
krS2Polygonx_het03 = cat(1,krS2Anterior_het03(:,1), flip(krS3Anterior_het03(:,1)))';
krS2Polygony_het03 = cat(1,krS2Anterior_het03(:,2), flip(krS3Anterior_het03(:,2)))';
krS2Pgon_het03 = polyshape(krS2Polygonx_het03,krS2Polygony_het03);
krstripe2_het03 = find(inpolygon(krcx_het03,krcy_het03,krS2Polygonx_het03,krS2Polygony_het03)== 1);
krstripe2Act_het03 = intersect(krstripe2_het03, kract_het03);
%krS6
krS6Polygonx_het03 = cat(1,krS6Anterior_het03(:,1), flip(krS7Anterior_het03(:,1)))';
krS6Polygony_het03 = cat(1,krS6Anterior_het03(:,2), flip(krS7Anterior_het03(:,2)))';
krS6Pgon_het03 = polyshape(krS6Polygonx_het03,krS6Polygony_het03);
krstripe6_het03 = find(inpolygon(krcx_het03,krcy_het03,krS6Polygonx_het03,krS6Polygony_het03)== 1);
krstripe6Act_het03 = intersect(krstripe6_het03, kract_het03);

%krhet04
%krS1
krS1Polygonx_het04 = cat(1,krS1Anterior_het04(:,1), flip(krS2Anterior_het04(:,1)))';
krS1Polygony_het04 = cat(1,krS1Anterior_het04(:,2), flip(krS2Anterior_het04(:,2)))';
krS1Pgon_het04 = polyshape(krS1Polygonx_het04,krS1Polygony_het04);
krstripe1_het04 = find(inpolygon(krcx_het04,krcy_het04,krS1Polygonx_het04,krS1Polygony_het04)== 1);
krstripe1Act_het04 = intersect(krstripe1_het04, kract_het04);
%krS2
krS2Polygonx_het04 = cat(1,krS2Anterior_het04(:,1), flip(krS3Anterior_het04(:,1)))';
krS2Polygony_het04 = cat(1,krS2Anterior_het04(:,2), flip(krS3Anterior_het04(:,2)))';
krS2Pgon_het04 = polyshape(krS2Polygonx_het04,krS2Polygony_het04);
krstripe2_het04 = find(inpolygon(krcx_het04,krcy_het04,krS2Polygonx_het04,krS2Polygony_het04)== 1);
krstripe2Act_het04 = intersect(krstripe2_het04, kract_het04);
%krS6
krS6Polygonx_het04 = cat(1,krS6Anterior_het04(:,1), flip(krS7Anterior_het04(:,1)))';
krS6Polygony_het04 = cat(1,krS6Anterior_het04(:,2), flip(krS7Anterior_het04(:,2)))';
krS6Pgon_het04 = polyshape(krS6Polygonx_het04,krS6Polygony_het04);
krstripe6_het04 = find(inpolygon(krcx_het04,krcy_het04,krS6Polygonx_het04,krS6Polygony_het04)== 1);
krstripe6Act_het04 = intersect(krstripe6_het04, kract_het04);

%krhet05
%krS1
krS1Polygonx_het05 = cat(1,krS1Anterior_het05(:,1), flip(krS2Anterior_het05(:,1)))';
krS1Polygony_het05 = cat(1,krS1Anterior_het05(:,2), flip(krS2Anterior_het05(:,2)))';
krS1Pgon_het05 = polyshape(krS1Polygonx_het05,krS1Polygony_het05);
krstripe1_het05 = find(inpolygon(krcx_het05,krcy_het05,krS1Polygonx_het05,krS1Polygony_het05)== 1);
krstripe1Act_het05 = intersect(krstripe1_het05, kract_het05);
%krS2
krS2Polygonx_het05 = cat(1,krS2Anterior_het05(:,1), flip(krS3Anterior_het05(:,1)))';
krS2Polygony_het05 = cat(1,krS2Anterior_het05(:,2), flip(krS3Anterior_het05(:,2)))';
krS2Pgon_het05 = polyshape(krS2Polygonx_het05,krS2Polygony_het05);
krstripe2_het05 = find(inpolygon(krcx_het05,krcy_het05,krS2Polygonx_het05,krS2Polygony_het05)== 1);
krstripe2Act_het05 = intersect(krstripe2_het05, kract_het05);
%krS6
krS6Polygonx_het05 = cat(1,krS6Anterior_het05(:,1), flip(krS7Anterior_het05(:,1)))';
krS6Polygony_het05 = cat(1,krS6Anterior_het05(:,2), flip(krS7Anterior_het05(:,2)))';
krS6Pgon_het05 = polyshape(krS6Polygonx_het05,krS6Polygony_het05);
krstripe6_het05 = find(inpolygon(krcx_het05,krcy_het05,krS6Polygonx_het05,krS6Polygony_het05)== 1);
krstripe6Act_het05 = intersect(krstripe6_het05, kract_het05);

%wt01
%S3
S3Polygonx_wt01 = cat(1,S3Anterior_wt01(:,1), flip(S4Anterior_wt01(:,1)))';
S3Polygony_wt01 = cat(1,S3Anterior_wt01(:,2), flip(S4Anterior_wt01(:,2)))';
S3Pgon_wt01 = polyshape(S3Polygonx_wt01,S3Polygony_wt01);
stripe3_wt01 = find(inpolygon(cx_wt01,cy_wt01,S3Polygonx_wt01,S3Polygony_wt01)== 1);
stripe3Act_wt01 = intersect(stripe3_wt01, act_wt01);
%S4
S4Polygonx_wt01 = cat(1,S4Anterior_wt01(:,1), flip(S5Anterior_wt01(:,1)))';
S4Polygony_wt01 = cat(1,S4Anterior_wt01(:,2), flip(S5Anterior_wt01(:,2)))';
S4Pgon_wt01 = polyshape(S4Polygonx_wt01,S4Polygony_wt01);
stripe4_wt01 = find(inpolygon(cx_wt01,cy_wt01,S4Polygonx_wt01,S4Polygony_wt01)== 1);
stripe4Act_wt01 = intersect(stripe4_wt01, act_wt01);
%S5
S5Polygonx_wt01 = cat(1,S5Anterior_wt01(:,1), flip(S6Anterior_wt01(:,1)))';
S5Polygony_wt01 = cat(1,S5Anterior_wt01(:,2), flip(S6Anterior_wt01(:,2)))';
S5Pgon_wt01 = polyshape(S5Polygonx_wt01,S5Polygony_wt01);
stripe5_wt01 = find(inpolygon(cx_wt01,cy_wt01,S5Polygonx_wt01,S5Polygony_wt01)== 1);
stripe5Act_wt01 = intersect(stripe5_wt01, act_wt01);

%wt02
%S3
S3Polygonx_wt02 = cat(1,S3Anterior_wt02(:,1), flip(S4Anterior_wt02(:,1)))';
S3Polygony_wt02 = cat(1,S3Anterior_wt02(:,2), flip(S4Anterior_wt02(:,2)))';
S3Pgon_wt02 = polyshape(S3Polygonx_wt02,S3Polygony_wt02);
stripe3_wt02 = find(inpolygon(cx_wt02,cy_wt02,S3Polygonx_wt02,S3Polygony_wt02)== 1);
stripe3Act_wt02 = intersect(stripe3_wt02, act_wt02);
%S4
S4Polygonx_wt02 = cat(1,S4Anterior_wt02(:,1), flip(S5Anterior_wt02(:,1)))';
S4Polygony_wt02 = cat(1,S4Anterior_wt02(:,2), flip(S5Anterior_wt02(:,2)))';
S4Pgon_wt02 = polyshape(S4Polygonx_wt02,S4Polygony_wt02);
stripe4_wt02 = find(inpolygon(cx_wt02,cy_wt02,S4Polygonx_wt02,S4Polygony_wt02)== 1);
stripe4Act_wt02 = intersect(stripe4_wt02, act_wt02);
%S5
S5Polygonx_wt02 = cat(1,S5Anterior_wt02(:,1), flip(S6Anterior_wt02(:,1)))';
S5Polygony_wt02 = cat(1,S5Anterior_wt02(:,2), flip(S6Anterior_wt02(:,2)))';
S5Pgon_wt02 = polyshape(S5Polygonx_wt02,S5Polygony_wt02);
stripe5_wt02 = find(inpolygon(cx_wt02,cy_wt02,S5Polygonx_wt02,S5Polygony_wt02)== 1);
stripe5Act_wt02 = intersect(stripe5_wt02, act_wt02);

%wt03
%S3
S3Polygonx_wt03 = cat(1,S3Anterior_wt03(:,1), flip(S4Anterior_wt03(:,1)))';
S3Polygony_wt03 = cat(1,S3Anterior_wt03(:,2), flip(S4Anterior_wt03(:,2)))';
S3Pgon_wt03 = polyshape(S3Polygonx_wt03,S3Polygony_wt03);
stripe3_wt03 = find(inpolygon(cx_wt03,cy_wt03,S3Polygonx_wt03,S3Polygony_wt03)== 1);
stripe3Act_wt03 = intersect(stripe3_wt03, act_wt03);
%S4
S4Polygonx_wt03 = cat(1,S4Anterior_wt03(:,1), flip(S5Anterior_wt03(:,1)))';
S4Polygony_wt03 = cat(1,S4Anterior_wt03(:,2), flip(S5Anterior_wt03(:,2)))';
S4Pgon_wt03 = polyshape(S4Polygonx_wt03,S4Polygony_wt03);
stripe4_wt03 = find(inpolygon(cx_wt03,cy_wt03,S4Polygonx_wt03,S4Polygony_wt03)== 1);
stripe4Act_wt03 = intersect(stripe4_wt03, act_wt03);
%S5
S5Polygonx_wt03 = cat(1,S5Anterior_wt03(:,1), flip(S6Anterior_wt03(:,1)))';
S5Polygony_wt03 = cat(1,S5Anterior_wt03(:,2), flip(S6Anterior_wt03(:,2)))';
S5Pgon_wt03 = polyshape(S5Polygonx_wt03,S5Polygony_wt03);
stripe5_wt03 = find(inpolygon(cx_wt03,cy_wt03,S5Polygonx_wt03,S5Polygony_wt03)== 1);
stripe5Act_wt03 = intersect(stripe5_wt03, act_wt03);

%krwt01
%krS3
krS3Polygonx_wt01 = cat(1,krS3Anterior_wt01(:,1), flip(krS4Anterior_wt01(:,1)))';
krS3Polygony_wt01 = cat(1,krS3Anterior_wt01(:,2), flip(krS4Anterior_wt01(:,2)))';
krS3Pgon_wt01 = polyshape(krS3Polygonx_wt01,krS3Polygony_wt01);
krstripe3_wt01 = find(inpolygon(krcx_wt01,krcy_wt01,krS3Polygonx_wt01,krS3Polygony_wt01)== 1);
krstripe3Act_wt01 = intersect(krstripe3_wt01, kract_wt01);
%krS4
krS4Polygonx_wt01 = cat(1,krS4Anterior_wt01(:,1), flip(krS5Anterior_wt01(:,1)))';
krS4Polygony_wt01 = cat(1,krS4Anterior_wt01(:,2), flip(krS5Anterior_wt01(:,2)))';
krS4Pgon_wt01 = polyshape(krS4Polygonx_wt01,krS4Polygony_wt01);
krstripe4_wt01 = find(inpolygon(krcx_wt01,krcy_wt01,krS4Polygonx_wt01,krS4Polygony_wt01)== 1);
krstripe4Act_wt01 = intersect(krstripe4_wt01, kract_wt01);
%krS5
krS5Polygonx_wt01 = cat(1,krS5Anterior_wt01(:,1), flip(krS6Anterior_wt01(:,1)))';
krS5Polygony_wt01 = cat(1,krS5Anterior_wt01(:,2), flip(krS6Anterior_wt01(:,2)))';
krS5Pgon_wt01 = polyshape(krS5Polygonx_wt01,krS5Polygony_wt01);
krstripe5_wt01 = find(inpolygon(krcx_wt01,krcy_wt01,krS5Polygonx_wt01,krS5Polygony_wt01)== 1);
krstripe5Act_wt01 = intersect(krstripe5_wt01, kract_wt01);

%krwt02
%krS3
krS3Polygonx_wt02 = cat(1,krS3Anterior_wt02(:,1), flip(krS4Anterior_wt02(:,1)))';
krS3Polygony_wt02 = cat(1,krS3Anterior_wt02(:,2), flip(krS4Anterior_wt02(:,2)))';
krS3Pgon_wt02 = polyshape(krS3Polygonx_wt02,krS3Polygony_wt02);
krstripe3_wt02 = find(inpolygon(krcx_wt02,krcy_wt02,krS3Polygonx_wt02,krS3Polygony_wt02)== 1);
krstripe3Act_wt02 = intersect(krstripe3_wt02, kract_wt02);
%krS4
krS4Polygonx_wt02 = cat(1,krS4Anterior_wt02(:,1), flip(krS5Anterior_wt02(:,1)))';
krS4Polygony_wt02 = cat(1,krS4Anterior_wt02(:,2), flip(krS5Anterior_wt02(:,2)))';
krS4Pgon_wt02 = polyshape(krS4Polygonx_wt02,krS4Polygony_wt02);
krstripe4_wt02 = find(inpolygon(krcx_wt02,krcy_wt02,krS4Polygonx_wt02,krS4Polygony_wt02)== 1);
krstripe4Act_wt02 = intersect(krstripe4_wt02, kract_wt02);
%krS5
krS5Polygonx_wt02 = cat(1,krS5Anterior_wt02(:,1), flip(krS6Anterior_wt02(:,1)))';
krS5Polygony_wt02 = cat(1,krS5Anterior_wt02(:,2), flip(krS6Anterior_wt02(:,2)))';
krS5Pgon_wt02 = polyshape(krS5Polygonx_wt02,krS5Polygony_wt02);
krstripe5_wt02 = find(inpolygon(krcx_wt02,krcy_wt02,krS5Polygonx_wt02,krS5Polygony_wt02)== 1);
krstripe5Act_wt02 = intersect(krstripe5_wt02, kract_wt02);

%krwt03
%krS3
krS3Polygonx_wt03 = cat(1,krS3Anterior_wt03(:,1), flip(krS4Anterior_wt03(:,1)))';
krS3Polygony_wt03 = cat(1,krS3Anterior_wt03(:,2), flip(krS4Anterior_wt03(:,2)))';
krS3Pgon_wt03 = polyshape(krS3Polygonx_wt03,krS3Polygony_wt03);
krstripe3_wt03 = find(inpolygon(krcx_wt03,krcy_wt03,krS3Polygonx_wt03,krS3Polygony_wt03)== 1);
krstripe3Act_wt03 = intersect(krstripe3_wt03, kract_wt03);
%krS4
krS4Polygonx_wt03 = cat(1,krS4Anterior_wt03(:,1), flip(krS5Anterior_wt03(:,1)))';
krS4Polygony_wt03 = cat(1,krS4Anterior_wt03(:,2), flip(krS5Anterior_wt03(:,2)))';
krS4Pgon_wt03 = polyshape(krS4Polygonx_wt03,krS4Polygony_wt03);
krstripe4_wt03 = find(inpolygon(krcx_wt03,krcy_wt03,krS4Polygonx_wt03,krS4Polygony_wt03)== 1);
krstripe4Act_wt03 = intersect(krstripe4_wt03, kract_wt03);
%krS5
krS5Polygonx_wt03 = cat(1,krS5Anterior_wt03(:,1), flip(krS6Anterior_wt03(:,1)))';
krS5Polygony_wt03 = cat(1,krS5Anterior_wt03(:,2), flip(krS6Anterior_wt03(:,2)))';
krS5Pgon_wt03 = polyshape(krS5Polygonx_wt03,krS5Polygony_wt03);
krstripe5_wt03 = find(inpolygon(krcx_wt03,krcy_wt03,krS5Polygonx_wt03,krS5Polygony_wt03)== 1);
krstripe5Act_wt03 = intersect(krstripe5_wt03, kract_wt03);

%krwt04
%krS3
krS3Polygonx_wt04 = cat(1,krS3Anterior_wt04(:,1), flip(krS4Anterior_wt04(:,1)))';
krS3Polygony_wt04 = cat(1,krS3Anterior_wt04(:,2), flip(krS4Anterior_wt04(:,2)))';
krS3Pgon_wt04 = polyshape(krS3Polygonx_wt04,krS3Polygony_wt04);
krstripe3_wt04 = find(inpolygon(krcx_wt04,krcy_wt04,krS3Polygonx_wt04,krS3Polygony_wt04)== 1);
krstripe3Act_wt04 = intersect(krstripe3_wt04, kract_wt04);
%krS4
krS4Polygonx_wt04 = cat(1,krS4Anterior_wt04(:,1), flip(krS5Anterior_wt04(:,1)))';
krS4Polygony_wt04 = cat(1,krS4Anterior_wt04(:,2), flip(krS5Anterior_wt04(:,2)))';
krS4Pgon_wt04 = polyshape(krS4Polygonx_wt04,krS4Polygony_wt04);
krstripe4_wt04 = find(inpolygon(krcx_wt04,krcy_wt04,krS4Polygonx_wt04,krS4Polygony_wt04)== 1);
krstripe4Act_wt04 = intersect(krstripe4_wt04, kract_wt04);
%krS5
krS5Polygonx_wt04 = cat(1,krS5Anterior_wt04(:,1), flip(krS6Anterior_wt04(:,1)))';
krS5Polygony_wt04 = cat(1,krS5Anterior_wt04(:,2), flip(krS6Anterior_wt04(:,2)))';
krS5Pgon_wt04 = polyshape(krS5Polygonx_wt04,krS5Polygony_wt04);
krstripe5_wt04 = find(inpolygon(krcx_wt04,krcy_wt04,krS5Polygonx_wt04,krS5Polygony_wt04)== 1);
krstripe5Act_wt04 = intersect(krstripe5_wt04, kract_wt04);

%krhet01
%krS3
krS3Polygonx_het01 = cat(1,krS3Anterior_het01(:,1), flip(krS4Anterior_het01(:,1)))';
krS3Polygony_het01 = cat(1,krS3Anterior_het01(:,2), flip(krS4Anterior_het01(:,2)))';
krS3Pgon_het01 = polyshape(krS3Polygonx_het01,krS3Polygony_het01);
krstripe3_het01 = find(inpolygon(krcx_het01,krcy_het01,krS3Polygonx_het01,krS3Polygony_het01)== 1);
krstripe3Act_het01 = intersect(krstripe3_het01, kract_het01);
%krS4
krS4Polygonx_het01 = cat(1,krS4Anterior_het01(:,1), flip(krS5Anterior_het01(:,1)))';
krS4Polygony_het01 = cat(1,krS4Anterior_het01(:,2), flip(krS5Anterior_het01(:,2)))';
krS4Pgon_het01 = polyshape(krS4Polygonx_het01,krS4Polygony_het01);
krstripe4_het01 = find(inpolygon(krcx_het01,krcy_het01,krS4Polygonx_het01,krS4Polygony_het01)== 1);
krstripe4Act_het01 = intersect(krstripe4_het01, kract_het01);
%krS5
krS5Polygonx_het01 = cat(1,krS5Anterior_het01(:,1), flip(krS6Anterior_het01(:,1)))';
krS5Polygony_het01 = cat(1,krS5Anterior_het01(:,2), flip(krS6Anterior_het01(:,2)))';
krS5Pgon_het01 = polyshape(krS5Polygonx_het01,krS5Polygony_het01);
krstripe5_het01 = find(inpolygon(krcx_het01,krcy_het01,krS5Polygonx_het01,krS5Polygony_het01)== 1);
krstripe5Act_het01 = intersect(krstripe5_het01, kract_het01);

%krhet02
%krS3
krS3Polygonx_het02 = cat(1,krS3Anterior_het02(:,1), flip(krS4Anterior_het02(:,1)))';
krS3Polygony_het02 = cat(1,krS3Anterior_het02(:,2), flip(krS4Anterior_het02(:,2)))';
krS3Pgon_het02 = polyshape(krS3Polygonx_het02,krS3Polygony_het02);
krstripe3_het02 = find(inpolygon(krcx_het02,krcy_het02,krS3Polygonx_het02,krS3Polygony_het02)== 1);
krstripe3Act_het02 = intersect(krstripe3_het02, kract_het02);
%krS4
krS4Polygonx_het02 = cat(1,krS4Anterior_het02(:,1), flip(krS5Anterior_het02(:,1)))';
krS4Polygony_het02 = cat(1,krS4Anterior_het02(:,2), flip(krS5Anterior_het02(:,2)))';
krS4Pgon_het02 = polyshape(krS4Polygonx_het02,krS4Polygony_het02);
krstripe4_het02 = find(inpolygon(krcx_het02,krcy_het02,krS4Polygonx_het02,krS4Polygony_het02)== 1);
krstripe4Act_het02 = intersect(krstripe4_het02, kract_het02);
%krS5
krS5Polygonx_het02 = cat(1,krS5Anterior_het02(:,1), flip(krS6Anterior_het02(:,1)))';
krS5Polygony_het02 = cat(1,krS5Anterior_het02(:,2), flip(krS6Anterior_het02(:,2)))';
krS5Pgon_het02 = polyshape(krS5Polygonx_het02,krS5Polygony_het02);
krstripe5_het02 = find(inpolygon(krcx_het02,krcy_het02,krS5Polygonx_het02,krS5Polygony_het02)== 1);
krstripe5Act_het02 = intersect(krstripe5_het02, kract_het02);

%krhet03
%krS3
krS3Polygonx_het03 = cat(1,krS3Anterior_het03(:,1), flip(krS4Anterior_het03(:,1)))';
krS3Polygony_het03 = cat(1,krS3Anterior_het03(:,2), flip(krS4Anterior_het03(:,2)))';
krS3Pgon_het03 = polyshape(krS3Polygonx_het03,krS3Polygony_het03);
krstripe3_het03 = find(inpolygon(krcx_het03,krcy_het03,krS3Polygonx_het03,krS3Polygony_het03)== 1);
krstripe3Act_het03 = intersect(krstripe3_het03, kract_het03);
%krS4
krS4Polygonx_het03 = cat(1,krS4Anterior_het03(:,1), flip(krS5Anterior_het03(:,1)))';
krS4Polygony_het03 = cat(1,krS4Anterior_het03(:,2), flip(krS5Anterior_het03(:,2)))';
krS4Pgon_het03 = polyshape(krS4Polygonx_het03,krS4Polygony_het03);
krstripe4_het03 = find(inpolygon(krcx_het03,krcy_het03,krS4Polygonx_het03,krS4Polygony_het03)== 1);
krstripe4Act_het03 = intersect(krstripe4_het03, kract_het03);
%krS5
krS5Polygonx_het03 = cat(1,krS5Anterior_het03(:,1), flip(krS6Anterior_het03(:,1)))';
krS5Polygony_het03 = cat(1,krS5Anterior_het03(:,2), flip(krS6Anterior_het03(:,2)))';
krS5Pgon_het03 = polyshape(krS5Polygonx_het03,krS5Polygony_het03);
krstripe5_het03 = find(inpolygon(krcx_het03,krcy_het03,krS5Polygonx_het03,krS5Polygony_het03)== 1);
krstripe5Act_het03 = intersect(krstripe5_het03, kract_het03);

%krhet04
%krS3
krS3Polygonx_het04 = cat(1,krS3Anterior_het04(:,1), flip(krS4Anterior_het04(:,1)))';
krS3Polygony_het04 = cat(1,krS3Anterior_het04(:,2), flip(krS4Anterior_het04(:,2)))';
krS3Pgon_het04 = polyshape(krS3Polygonx_het04,krS3Polygony_het04);
krstripe3_het04 = find(inpolygon(krcx_het04,krcy_het04,krS3Polygonx_het04,krS3Polygony_het04)== 1);
krstripe3Act_het04 = intersect(krstripe3_het04, kract_het04);
%krS4
krS4Polygonx_het04 = cat(1,krS4Anterior_het04(:,1), flip(krS5Anterior_het04(:,1)))';
krS4Polygony_het04 = cat(1,krS4Anterior_het04(:,2), flip(krS5Anterior_het04(:,2)))';
krS4Pgon_het04 = polyshape(krS4Polygonx_het04,krS4Polygony_het04);
krstripe4_het04 = find(inpolygon(krcx_het04,krcy_het04,krS4Polygonx_het04,krS4Polygony_het04)== 1);
krstripe4Act_het04 = intersect(krstripe4_het04, kract_het04);
%krS5
krS5Polygonx_het04 = cat(1,krS5Anterior_het04(:,1), flip(krS6Anterior_het04(:,1)))';
krS5Polygony_het04 = cat(1,krS5Anterior_het04(:,2), flip(krS6Anterior_het04(:,2)))';
krS5Pgon_het04 = polyshape(krS5Polygonx_het04,krS5Polygony_het04);
krstripe5_het04 = find(inpolygon(krcx_het04,krcy_het04,krS5Polygonx_het04,krS5Polygony_het04)== 1);
krstripe5Act_het04 = intersect(krstripe5_het04, kract_het04);

%krhet05
%krS3
krS3Polygonx_het05 = cat(1,krS3Anterior_het05(:,1), flip(krS4Anterior_het05(:,1)))';
krS3Polygony_het05 = cat(1,krS3Anterior_het05(:,2), flip(krS4Anterior_het05(:,2)))';
krS3Pgon_het05 = polyshape(krS3Polygonx_het05,krS3Polygony_het05);
krstripe3_het05 = find(inpolygon(krcx_het05,krcy_het05,krS3Polygonx_het05,krS3Polygony_het05)== 1);
krstripe3Act_het05 = intersect(krstripe3_het05, kract_het05);
%krS4
krS4Polygonx_het05 = cat(1,krS4Anterior_het05(:,1), flip(krS5Anterior_het05(:,1)))';
krS4Polygony_het05 = cat(1,krS4Anterior_het05(:,2), flip(krS5Anterior_het05(:,2)))';
krS4Pgon_het05 = polyshape(krS4Polygonx_het05,krS4Polygony_het05);
krstripe4_het05 = find(inpolygon(krcx_het05,krcy_het05,krS4Polygonx_het05,krS4Polygony_het05)== 1);
krstripe4Act_het05 = intersect(krstripe4_het05, kract_het05);
%krS5
krS5Polygonx_het05 = cat(1,krS5Anterior_het05(:,1), flip(krS6Anterior_het05(:,1)))';
krS5Polygony_het05 = cat(1,krS5Anterior_het05(:,2), flip(krS6Anterior_het05(:,2)))';
krS5Pgon_het05 = polyshape(krS5Polygonx_het05,krS5Polygony_het05);
krstripe5_het05 = find(inpolygon(krcx_het05,krcy_het05,krS5Polygonx_het05,krS5Polygony_het05)== 1);
krstripe5Act_het05 = intersect(krstripe5_het05, kract_het05);


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

%% width of embryos

extrema_wt01 = rpo_wt01.Extrema;
wtip_wt01 = [min(extrema_wt01(:,2)) max(extrema_wt01(:,2))];
extrema_wt02 = rpo_wt02.Extrema;
wtip_wt02 = [min(extrema_wt02(:,2)) max(extrema_wt02(:,2))];
extrema_wt03 = rpo_wt03.Extrema;
wtip_wt03 = [min(extrema_wt03(:,2)) max(extrema_wt03(:,2))];
krextrema_wt01 = krrpo_wt01.Extrema;
krwtip_wt01 = [min(krextrema_wt01(:,2)) max(krextrema_wt01(:,2))];
krextrema_wt02 = krrpo_wt02.Extrema;
krwtip_wt02 = [min(krextrema_wt02(:,2)) max(krextrema_wt02(:,2))];
krextrema_wt03 = krrpo_wt03.Extrema;
krwtip_wt03 = [min(krextrema_wt03(:,2)) max(krextrema_wt03(:,2))];
krextrema_wt04 = krrpo_wt04.Extrema;
krwtip_wt04 = [min(krextrema_wt04(:,2)) max(krextrema_wt04(:,2))];

krextrema_het01 = krrpo_het01.Extrema;
krwtip_het01 = [min(krextrema_het01(:,2)) max(krextrema_het01(:,2))];
krextrema_het02 = krrpo_het02.Extrema;
krwtip_het02 = [min(krextrema_het02(:,2)) max(krextrema_het02(:,2))];
krextrema_het03 = krrpo_het03.Extrema;
krwtip_het03 = [min(krextrema_het03(:,2)) max(krextrema_het03(:,2))];
krextrema_het04 = krrpo_het04.Extrema;
krwtip_het04 = [min(krextrema_het04(:,2)) max(krextrema_het04(:,2))];
krextrema_het05 = krrpo_het05.Extrema;
krwtip_het05 = [min(krextrema_het05(:,2)) max(krextrema_het05(:,2))];

%% width of embryos
width_wt01 = abs(wtip_wt01(1)-wtip_wt01(2));
width_wt02 = abs(wtip_wt02(1)-wtip_wt02(2));
width_wt03 = abs(wtip_wt03(1)-wtip_wt03(2));

krwidth_wt01 = abs(krwtip_wt01(1)-krwtip_wt01(2));
krwidth_wt02 = abs(krwtip_wt02(1)-krwtip_wt02(2));
krwidth_wt03 = abs(krwtip_wt03(1)-krwtip_wt03(2));
krwidth_wt04 = abs(krwtip_wt04(1)-krwtip_wt04(2));

krwidth_het01 = abs(krwtip_het01(1)-krwtip_het01(2));
krwidth_het02 = abs(krwtip_het02(1)-krwtip_het02(2));
krwidth_het03 = abs(krwtip_het03(1)-krwtip_het03(2));
krwidth_het04 = abs(krwtip_het04(1)-krwtip_het04(2));
krwidth_het05 = abs(krwtip_het05(1)-krwtip_het05(2));

%% check extreme points
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

%% flip embryos - notD: tips need to need min and max for anterior and posterior
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

%% act at 2nd half NC14

%stripe1
actS1second_wt01 = intersect(act2nd_wt01, stripe1_wt01);
actS1second_wt02 = intersect(act2nd_wt02, stripe1_wt02);
actS1second_wt03 = intersect(act2nd_wt03, stripe1_wt03);
kractS1second_wt01 = intersect(kract2nd_wt01, krstripe1_wt01);
kractS1second_wt02 = intersect(kract2nd_wt02, krstripe1_wt02);
kractS1second_wt03 = intersect(kract2nd_wt03, krstripe1_wt03);
kractS1second_wt04 = intersect(kract2nd_wt04, krstripe1_wt04);
kractS1second_het01 = intersect(kract2nd_het01, krstripe1_het01);
kractS1second_het02 = intersect(kract2nd_het02, krstripe1_het02);
kractS1second_het03 = intersect(kract2nd_het03, krstripe1_het03);
kractS1second_het04 = intersect(kract2nd_het04, krstripe1_het04);
kractS1second_het05 = intersect(kract2nd_het05, krstripe1_het05);

%stripe2
actS2second_wt01 = intersect(act2nd_wt01, stripe2_wt01);
actS2second_wt02 = intersect(act2nd_wt02, stripe2_wt02);
actS2second_wt03 = intersect(act2nd_wt03, stripe2_wt03);
kractS2second_wt01 = intersect(kract2nd_wt01, krstripe2_wt01);
kractS2second_wt02 = intersect(kract2nd_wt02, krstripe2_wt02);
kractS2second_wt03 = intersect(kract2nd_wt03, krstripe2_wt03);
kractS2second_wt04 = intersect(kract2nd_wt04, krstripe2_wt04);
kractS2second_het01 = intersect(kract2nd_het01, krstripe2_het01);
kractS2second_het02 = intersect(kract2nd_het02, krstripe2_het02);
kractS2second_het03 = intersect(kract2nd_het03, krstripe2_het03);
kractS2second_het04 = intersect(kract2nd_het04, krstripe2_het04);
kractS2second_het05 = intersect(kract2nd_het05, krstripe2_het05);

%stripe3
actS3second_wt01 = intersect(act2nd_wt01, stripe3_wt01);
actS3second_wt02 = intersect(act2nd_wt02, stripe3_wt02);
actS3second_wt03 = intersect(act2nd_wt03, stripe3_wt03);
kractS3second_wt01 = intersect(kract2nd_wt01, krstripe3_wt01);
kractS3second_wt02 = intersect(kract2nd_wt02, krstripe3_wt02);
kractS3second_wt03 = intersect(kract2nd_wt03, krstripe3_wt03);
kractS3second_wt04 = intersect(kract2nd_wt04, krstripe3_wt04);
kractS3second_het01 = intersect(kract2nd_het01, krstripe3_het01);
kractS3second_het02 = intersect(kract2nd_het02, krstripe3_het02);
kractS3second_het03 = intersect(kract2nd_het03, krstripe3_het03);
kractS3second_het04 = intersect(kract2nd_het04, krstripe3_het04);
kractS3second_het05 = intersect(kract2nd_het05, krstripe3_het05);

%stripe4
actS4second_wt01 = intersect(act2nd_wt01, stripe4_wt01);
actS4second_wt02 = intersect(act2nd_wt02, stripe4_wt02);
actS4second_wt03 = intersect(act2nd_wt03, stripe4_wt03);
kractS4second_wt01 = intersect(kract2nd_wt01, krstripe4_wt01);
kractS4second_wt02 = intersect(kract2nd_wt02, krstripe4_wt02);
kractS4second_wt03 = intersect(kract2nd_wt03, krstripe4_wt03);
kractS4second_wt04 = intersect(kract2nd_wt04, krstripe4_wt04);
kractS4second_het01 = intersect(kract2nd_het01, krstripe4_het01);
kractS4second_het02 = intersect(kract2nd_het02, krstripe4_het02);
kractS4second_het03 = intersect(kract2nd_het03, krstripe4_het03);
kractS4second_het04 = intersect(kract2nd_het04, krstripe4_het04);
kractS4second_het05 = intersect(kract2nd_het05, krstripe4_het05);

%stripe5
actS5second_wt01 = intersect(act2nd_wt01, stripe5_wt01);
actS5second_wt02 = intersect(act2nd_wt02, stripe5_wt02);
actS5second_wt03 = intersect(act2nd_wt03, stripe5_wt03);
kractS5second_wt01 = intersect(kract2nd_wt01, krstripe5_wt01);
kractS5second_wt02 = intersect(kract2nd_wt02, krstripe5_wt02);
kractS5second_wt03 = intersect(kract2nd_wt03, krstripe5_wt03);
kractS5second_wt04 = intersect(kract2nd_wt04, krstripe5_wt04);
kractS5second_het01 = intersect(kract2nd_het01, krstripe5_het01);
kractS5second_het02 = intersect(kract2nd_het02, krstripe5_het02);
kractS5second_het03 = intersect(kract2nd_het03, krstripe5_het03);
kractS5second_het04 = intersect(kract2nd_het04, krstripe5_het04);
kractS5second_het05 = intersect(kract2nd_het05, krstripe5_het05);

%stripe6
actS6second_wt01 = intersect(act2nd_wt01, stripe6_wt01);
actS6second_wt02 = intersect(act2nd_wt02, stripe6_wt02);
actS6second_wt03 = intersect(act2nd_wt03, stripe6_wt03);
kractS6second_wt01 = intersect(kract2nd_wt01, krstripe6_wt01);
kractS6second_wt02 = intersect(kract2nd_wt02, krstripe6_wt02);
kractS6second_wt03 = intersect(kract2nd_wt03, krstripe6_wt03);
kractS6second_wt04 = intersect(kract2nd_wt04, krstripe6_wt04);
kractS6second_het01 = intersect(kract2nd_het01, krstripe6_het01);
kractS6second_het02 = intersect(kract2nd_het02, krstripe6_het02);
kractS6second_het03 = intersect(kract2nd_het03, krstripe6_het03);
kractS6second_het04 = intersect(kract2nd_het04, krstripe6_het04);
kractS6second_het05 = intersect(kract2nd_het05, krstripe6_het05);

%% bin by 2% EL
binsize = 0.02;

%wt01
for i = 1:50
    binEL_wt01{i} = find(cxEL_wt01(:) > (i-1)*binsize & cxEL_wt01(:) <= i*binsize);
end
%wt02
for i = 1:50
    binEL_wt02{i} = find(cxEL_wt02(:) > (i-1)*binsize & cxEL_wt02(:) <= i*binsize);
end
%wt03
for i = 1:50
    binEL_wt03{i} = find(cxEL_wt03(:) > (i-1)*binsize & cxEL_wt03(:) <= i*binsize);
end
%krwt01
for i = 1:50
    krbinEL_wt01{i} = find(krcxEL_wt01(:) > (i-1)*binsize & krcxEL_wt01(:) <= i*binsize);
end
%krwt02
for i = 1:50
    krbinEL_wt02{i} = find(krcxEL_wt02(:) > (i-1)*binsize & krcxEL_wt02(:) <= i*binsize);
end
%krwt03
for i = 1:50
    krbinEL_wt03{i} = find(krcxEL_wt03(:) > (i-1)*binsize & krcxEL_wt03(:) <= i*binsize);
end
%krwt04
for i = 1:50
    krbinEL_wt04{i} = find(krcxEL_wt04(:) > (i-1)*binsize & krcxEL_wt04(:) <= i*binsize);
end
%krhet01
for i = 1:50
    krbinEL_het01{i} = find(krcxEL_het01(:) > (i-1)*binsize & krcxEL_het01(:) <= i*binsize);
end
%krhet02
for i = 1:50
    krbinEL_het02{i} = find(krcxEL_het02(:) > (i-1)*binsize & krcxEL_het02(:) <= i*binsize);
end
%krhet03
for i = 1:50
    krbinEL_het03{i} = find(krcxEL_het03(:) > (i-1)*binsize & krcxEL_het03(:) <= i*binsize);
end
%krhet04
for i = 1:50
    krbinEL_het04{i} = find(krcxEL_het04(:) > (i-1)*binsize & krcxEL_het04(:) <= i*binsize);
end
%krhet05
for i = 1:50
    krbinEL_het05{i} = find(krcxEL_het05(:) > (i-1)*binsize & krcxEL_het05(:) <= i*binsize);
end

%% 2nd half NC14 output
for i = 1:size(M_wt01,2)
    Output_wt01(i) = trapz(M_wt01(51:end,i));
end
for i = 1:size(M_wt02,2)
    Output_wt02(i) = trapz(M_wt02(51:end,i));
end
for i = 1:size(M_wt03,2)
    Output_wt03(i) = trapz(M_wt03(51:end,i));
end

for i = 1:size(krM_wt01,2)
    krOutput_wt01(i) = trapz(krM_wt01(51:end,i));
end
for i = 1:size(krM_wt02,2)
    krOutput_wt02(i) = trapz(krM_wt02(51:end,i));
end
for i = 1:size(krM_wt03,2)
    krOutput_wt03(i) = trapz(krM_wt03(51:end,i));
end
for i = 1:size(krM_wt04,2)
    krOutput_wt04(i) = trapz(krM_wt04(51:end,i));
end

for i = 1:size(krM_het01,2)
    krOutput_het01(i) = trapz(krM_het01(51:end,i));
end
for i = 1:size(krM_het02,2)
    krOutput_het02(i) = trapz(krM_het02(51:end,i));
end
for i = 1:size(krM_het03,2)
    krOutput_het03(i) = trapz(krM_het03(51:end,i));
end
for i = 1:size(krM_het04,2)
    krOutput_het04(i) = trapz(krM_het04(51:end,i));
end
for i = 1:size(krM_het05,2)
    krOutput_het05(i) = trapz(krM_het05(51:end,i));
end


%% individual nuclei boxplot for S3
OutputS3_wt = [Output_wt01(actS3second_wt01) Output_wt02(actS3second_wt02) Output_wt03(actS3second_wt03) ...
    krOutput_wt01(kractS3second_wt01) krOutput_wt02(kractS3second_wt02) krOutput_wt03(kractS3second_wt03) krOutput_wt04(kractS3second_wt04)];
OutputS3_het = [krOutput_het01(kractS3second_het01) krOutput_het02(kractS3second_het02) krOutput_het03(kractS3second_het03)...
    krOutput_het04(kractS3second_het04) krOutput_het05(kractS3second_het05)];

OutputS3P_wt = OutputS3_wt/10000;%for plot
OutputS3P_het = OutputS3_het/10000;%for plot

figure(171);
boxplot([OutputS3P_wt OutputS3P_het],...
    [ones(1,length(OutputS3_wt)), 2*ones(1,length(OutputS3_het))],...
    'Labels',{'\it\bfwt','\it\bfKr^{1}/+'},'symbol','');hold on
set(gca, 'TickLabelInterpreter', 'tex');
set(findobj(gca,'type','line'),'linew',2);
s1 = swarmchart(ones(1,length(OutputS3P_wt)),OutputS3P_wt,10,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
s1.XJitterWidth = 0.3;
s2 = swarmchart(2*ones(1,length(OutputS3P_het)),OutputS3P_het,10,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
s2.XJitterWidth = 0.3;
hold off
%% individual nuclei boxplot for S4
OutputS4_wt = [Output_wt01(actS4second_wt01) Output_wt02(actS4second_wt02) Output_wt03(actS4second_wt03) ...
    krOutput_wt01(kractS4second_wt01) krOutput_wt02(kractS4second_wt02) krOutput_wt03(kractS4second_wt03) krOutput_wt04(kractS4second_wt04)];
OutputS4_het = [krOutput_het01(kractS4second_het01) krOutput_het02(kractS4second_het02) krOutput_het03(kractS4second_het03)...
    krOutput_het04(kractS4second_het04) krOutput_het05(kractS4second_het05)];

OutputS4P_wt = OutputS4_wt/10000;%for plot
OutputS4P_het = OutputS4_het/10000;%for plot

figure(172);
boxplot([OutputS4P_wt OutputS4P_het],...
    [ones(1,length(OutputS4_wt)), 2*ones(1,length(OutputS4_het))],...
    'Labels',{'\it\bfwt','\it\bfKr^{1}/+'},'symbol','');hold on
set(gca, 'TickLabelInterpreter', 'tex');
set(findobj(gca,'type','line'),'linew',2);
s1 = swarmchart(ones(1,length(OutputS4P_wt)),OutputS4P_wt,10,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
s1.XJitterWidth = 0.3;
s2 = swarmchart(2*ones(1,length(OutputS4P_het)),OutputS4P_het,10,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
s2.XJitterWidth = 0.3;
hold off

%% individual nuclei boxplot for S5
OutputS5_wt = [Output_wt01(actS5second_wt01) Output_wt02(actS5second_wt02) Output_wt03(actS5second_wt03) ...
    krOutput_wt01(kractS5second_wt01) krOutput_wt02(kractS5second_wt02) krOutput_wt03(kractS5second_wt03) krOutput_wt04(kractS5second_wt04)];
OutputS5_het = [krOutput_het01(kractS5second_het01) krOutput_het02(kractS5second_het02) krOutput_het03(kractS5second_het03)...
    krOutput_het04(kractS5second_het04) krOutput_het05(kractS5second_het05)];

OutputS5P_wt = OutputS5_wt/10000;%for plot
OutputS5P_het = OutputS5_het/10000;%for plot

figure(173);
boxplot([OutputS5P_wt OutputS5P_het],...
    [ones(1,length(OutputS5_wt)), 2*ones(1,length(OutputS5_het))],...
    'Labels',{'\it\bfwt','\it\bfKr^{1}/+'},'symbol','');hold on
set(gca, 'TickLabelInterpreter', 'tex');
% set(gca, 'ActivePositionProperty', 'position', 'FontWeight','bold', 'FontSize',15);
set(findobj(gca,'type','line'),'linew',2);
s1 = swarmchart(ones(1,length(OutputS5P_wt)),OutputS5P_wt,10,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
s1.XJitterWidth = 0.3;
s2 = swarmchart(2*ones(1,length(OutputS5P_het)),OutputS5P_het,10,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
s2.XJitterWidth = 0.3;
hold off

%% individual nuclei boxplot for S2
OutputS2_wt = [Output_wt01(actS2second_wt01) Output_wt02(actS2second_wt02) Output_wt03(actS2second_wt03) ...
    krOutput_wt01(kractS2second_wt01) krOutput_wt02(kractS2second_wt02) krOutput_wt03(kractS2second_wt03) krOutput_wt04(kractS2second_wt04)];
OutputS2_het = [krOutput_het01(kractS2second_het01) krOutput_het02(kractS2second_het02) krOutput_het03(kractS2second_het03)...
    krOutput_het04(kractS2second_het04) krOutput_het05(kractS2second_het05)];

OutputS2P_wt = OutputS2_wt/10000;%for plot
OutputS2P_het = OutputS2_het/10000;%for plot

figure(174);
boxplot([OutputS2P_wt OutputS2P_het],...
    [ones(1,length(OutputS2_wt)), 2*ones(1,length(OutputS2_het))],...
    'Labels',{'\it\bfwt','\it\bfKr^{1}/+'},'symbol','');hold on
set(gca, 'TickLabelInterpreter', 'tex');
set(findobj(gca,'type','line'),'linew',2);
s1 = swarmchart(ones(1,length(OutputS2P_wt)),OutputS2P_wt,10,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
s1.XJitterWidth = 0.3;
s2 = swarmchart(2*ones(1,length(OutputS2P_het)),OutputS2P_het,10,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
s2.XJitterWidth = 0.3;
hold off
 %% Fraction of active nuclei
stripe1ActCountSecond_wt = [size(actS1second_wt01,2)/size(M_wt01,2),size(actS1second_wt02,2)/size(M_wt02,2),size(actS1second_wt03,2)/size(M_wt03,2),...
    size(kractS1second_wt01,2)/size(krM_wt01,2),size(kractS1second_wt02,2)/size(krM_wt02,2),size(kractS1second_wt03,2)/size(krM_wt03,2),size(kractS1second_wt04,2)/size(krM_wt04,2)];
stripe1ActCountSecond_het = [size(kractS1second_het01,2)/size(krM_het01,2),size(kractS1second_het02,2)/size(krM_het02,2),size(kractS1second_het03,2)/size(krM_het03,2),...
    size(kractS1second_het04,2)/size(krM_het04,2),size(kractS1second_het05,2)/size(krM_het05,2)];

stripe2ActCountSecond_wt = [size(actS2second_wt01,2)/size(M_wt01,2),size(actS2second_wt02,2)/size(M_wt02,2),size(actS2second_wt03,2)/size(M_wt03,2),...
    size(kractS2second_wt01,2)/size(krM_wt01,2),size(kractS2second_wt02,2)/size(krM_wt02,2),size(kractS2second_wt03,2)/size(krM_wt03,2),size(kractS2second_wt04,2)/size(krM_wt04,2)];
stripe2ActCountSecond_het = [size(kractS2second_het01,2)/size(krM_het01,2),size(kractS2second_het02,2)/size(krM_het02,2),size(kractS2second_het03,2)/size(krM_het03,2),...
    size(kractS2second_het04,2)/size(krM_het04,2),size(kractS2second_het05,2)/size(krM_het05,2)];

stripe3ActCountSecond_wt = [size(actS3second_wt01,2)/size(M_wt01,2),size(actS3second_wt02,2)/size(M_wt02,2),size(actS3second_wt03,2)/size(M_wt03,2),...
    size(kractS3second_wt01,2)/size(krM_wt01,2),size(kractS3second_wt02,2)/size(krM_wt02,2),size(kractS3second_wt03,2)/size(krM_wt03,2),size(kractS3second_wt04,2)/size(krM_wt04,2)];
stripe3ActCountSecond_het = [size(kractS3second_het01,2)/size(krM_het01,2),size(kractS3second_het02,2)/size(krM_het02,2),size(kractS3second_het03,2)/size(krM_het03,2),...
    size(kractS3second_het04,2)/size(krM_het04,2),size(kractS3second_het05,2)/size(krM_het05,2)];

stripe4ActCountSecond_wt = [size(actS4second_wt01,2)/size(M_wt01,2),size(actS4second_wt02,2)/size(M_wt02,2),size(actS4second_wt03,2)/size(M_wt03,2),...
    size(kractS4second_wt01,2)/size(krM_wt01,2),size(kractS4second_wt02,2)/size(krM_wt02,2),size(kractS4second_wt03,2)/size(krM_wt03,2),size(kractS4second_wt04,2)/size(krM_wt04,2)];
stripe4ActCountSecond_het = [size(kractS4second_het01,2)/size(krM_het01,2),size(kractS4second_het02,2)/size(krM_het02,2),size(kractS4second_het03,2)/size(krM_het03,2),...
    size(kractS4second_het04,2)/size(krM_het04,2),size(kractS4second_het05,2)/size(krM_het05,2)];

stripe5ActCountSecond_wt = [size(actS5second_wt01,2)/size(M_wt01,2),size(actS5second_wt02,2)/size(M_wt02,2),size(actS5second_wt03,2)/size(M_wt03,2),...
    size(kractS5second_wt01,2)/size(krM_wt01,2),size(kractS5second_wt02,2)/size(krM_wt02,2),size(kractS5second_wt03,2)/size(krM_wt03,2),size(kractS5second_wt04,2)/size(krM_wt04,2)];
stripe5ActCountSecond_het = [size(kractS5second_het01,2)/size(krM_het01,2),size(kractS5second_het02,2)/size(krM_het02,2),size(kractS5second_het03,2)/size(krM_het03,2),...
    size(kractS5second_het04,2)/size(krM_het04,2),size(kractS5second_het05,2)/size(krM_het05,2)];

stripe6ActCountSecond_wt = [size(actS6second_wt01,2)/size(M_wt01,2),size(actS6second_wt02,2)/size(M_wt02,2),size(actS6second_wt03,2)/size(M_wt03,2),...
    size(kractS6second_wt01,2)/size(krM_wt01,2),size(kractS6second_wt02,2)/size(krM_wt02,2),size(kractS6second_wt03,2)/size(krM_wt03,2),size(kractS6second_wt04,2)/size(krM_wt04,2)];
stripe6ActCountSecond_het = [size(kractS6second_het01,2)/size(krM_het01,2),size(kractS6second_het02,2)/size(krM_het02,2),size(kractS6second_het03,2)/size(krM_het03,2),...
    size(kractS6second_het04,2)/size(krM_het04,2),size(kractS6second_het05,2)/size(krM_het05,2)];


ScountMatN = nan(7,12);
ScountMatN(1:7,1)=stripe1ActCountSecond_wt';
ScountMatN(1:5,2)=stripe1ActCountSecond_het';
ScountMatN(1:7,3)=stripe2ActCountSecond_wt';
ScountMatN(1:5,4)=stripe2ActCountSecond_het';
ScountMatN(1:7,5)=stripe3ActCountSecond_wt';
ScountMatN(1:5,6)=stripe3ActCountSecond_het';
ScountMatN(1:7,7)=stripe4ActCountSecond_wt';
ScountMatN(1:5,8)=stripe4ActCountSecond_het';
ScountMatN(1:7,9)=stripe5ActCountSecond_wt';
ScountMatN(1:5,10)=stripe5ActCountSecond_het';
ScountMatN(1:7,11)=stripe6ActCountSecond_wt';
ScountMatN(1:5,12)=stripe6ActCountSecond_het';

ScountCellN{1} = ScountMatN(1:7,[1 3 5 7 9 11]);
ScountCellN{2} = ScountMatN(1:7,[2 4 6 8 10 12]);

%% individual fraction of active nuclei

condition_names =  {'\it\bfwt','\it\bfKr^{1}/+'};
 figure(221);
h = daboxplot({ScountCellN{1}(:,1)*100;ScountCellN{2}(:,1)*100},'boxalpha',0.7,...
    'xtlabels', condition_names); 
set(findobj(gca,'type','line'),'linew',2);


condition_names =  {'\it\bfwt','\it\bfKr^{1}/+'};
 figure(222);
h = daboxplot({ScountCellN{1}(:,2)*100;ScountCellN{2}(:,2)*100},'boxalpha',0.7,...
    'xtlabels', condition_names); 
set(findobj(gca,'type','line'),'linew',2);


condition_names =  {'\it\bfwt','\it\bfKr^{1}/+'};
 figure(223);
h = daboxplot({ScountCellN{1}(:,3)*100;ScountCellN{2}(:,3)*100},'boxalpha',0.7,...
    'xtlabels', condition_names); 
set(findobj(gca,'type','line'),'linew',2);


condition_names =  {'\it\bfwt','\it\bfKr^{1}/+'};
 figure(224);
h = daboxplot({ScountCellN{1}(:,4)*100;ScountCellN{2}(:,4)*100},'boxalpha',0.7,...
    'xtlabels', condition_names);
set(findobj(gca,'type','line'),'linew',2);


condition_names =  {'\it\bfwt','\it\bfKr^{1}/+'};
 figure(225);
h = daboxplot({ScountCellN{1}(:,5)*100;ScountCellN{2}(:,5)*100},'boxalpha',0.7,...
    'xtlabels', condition_names); 
set(findobj(gca,'type','line'),'linew',2);


condition_names =  {'\it\bfwt','\it\bfKr^{1}/+'};
 figure(226);
h = daboxplot({ScountCellN{1}(:,6)*100;ScountCellN{2}(:,6)*100},'boxalpha',0.7,...
    'xtlabels', condition_names); 
set(findobj(gca,'type','line'),'linew',2);



%% active duration of all nuclei and transcription rate

actDuraI_wt01 = {};
for i = 1: size(M_wt01,2)
    if ismember(i,act2nd_wt01) == 1
        tmp = find(M_wt01(51:end,i)>thr2);
        actDura_wt01(i) = length(tmp);
        actDuraI_wt01{i} = tmp';
    else
        actDura_wt01(i) = NaN;
        actDuraI_wt01{i} = [];
    end
end

actDuraI_wt02 = {};
for i = 1: size(M_wt02,2)
    if ismember(i,act2nd_wt02) == 1
        tmp = find(M_wt02(51:end,i)>thr2);
        actDura_wt02(i) = length(tmp);
        actDuraI_wt02{i} = tmp';
    else
        actDura_wt02(i) = NaN;
        actDuraI_wt02{i} = [];
    end
end

actDuraI_wt03 = {};
for i = 1: size(M_wt03,2)
    if ismember(i,act2nd_wt03) == 1
        tmp = find(M_wt03(51:end,i)>thr2);
        actDura_wt03(i) = length(tmp);
        actDuraI_wt03{i} = tmp';
    else
        actDura_wt03(i) = NaN;
        actDuraI_wt03{i} = [];
    end
end

kractDuraI_wt01 = {};
for i = 1: size(krM_wt01,2)
    if ismember(i,kract2nd_wt01) == 1
        tmp = find(krM_wt01(51:end,i)>thr2);
        kractDura_wt01(i) = length(tmp);
        kractDuraI_wt01{i} = tmp';
    else
        kractDura_wt01(i) = NaN;
        kractDuraI_wt01{i} = [];
    end
end

kractDuraI_wt02 = {};
for i = 1: size(krM_wt02,2)
    if ismember(i,kract2nd_wt02) == 1
        tmp = find(krM_wt02(51:end,i)>thr2);
        kractDura_wt02(i) = length(tmp);
        kractDuraI_wt02{i} = tmp';
    else
        kractDura_wt02(i) = NaN;
        kractDuraI_wt02{i} = [];
    end
end

kractDuraI_wt03 = {};
for i = 1: size(krM_wt03,2)
    if ismember(i,kract2nd_wt03) == 1
        tmp = find(krM_wt03(51:end,i)>thr2);
        kractDura_wt03(i) = length(tmp);
        kractDuraI_wt03{i} = tmp';
    else
        kractDura_wt03(i) = NaN;
        kractDuraI_wt03{i} = [];
    end
end

kractDuraI_wt04 = {};
for i = 1: size(krM_wt04,2)
    if ismember(i,kract2nd_wt04) == 1
        tmp = find(krM_wt04(51:end,i)>thr2);
        kractDura_wt04(i) = length(tmp);
        kractDuraI_wt04{i} = tmp';
    else
        kractDura_wt04(i) = NaN;
        kractDuraI_wt04{i} = [];
    end
end

kractDuraI_het01 = {};
for i = 1: size(krM_het01,2)
    if ismember(i,kract2nd_het01) == 1
        tmp = find(krM_het01(51:end,i)>thr2);
        kractDura_het01(i) = length(tmp);
        kractDuraI_het01{i} = tmp';
    else
        kractDura_het01(i) = NaN;
        kractDuraI_het01{i} = [];
    end
end

kractDuraI_het02 = {};
for i = 1: size(krM_het02,2)
    if ismember(i,kract2nd_het02) == 1
        tmp = find(krM_het02(51:end,i)>thr2);
        kractDura_het02(i) = length(tmp);
        kractDuraI_het02{i} = tmp';
    else
        kractDura_het02(i) = NaN;
        kractDuraI_het02{i} = [];
    end
end

kractDuraI_het03 = {};
for i = 1: size(krM_het03,2)
    if ismember(i,kract2nd_het03) == 1
        tmp = find(krM_het03(51:end,i)>thr2);
        kractDura_het03(i) = length(tmp);
        kractDuraI_het03{i} = tmp';
    else
        kractDura_het03(i) = NaN;
        kractDuraI_het03{i} = [];
    end
end

kractDuraI_het04 = {};
for i = 1: size(krM_het04,2)
    if ismember(i,kract2nd_het04) == 1
        tmp = find(krM_het04(51:end,i)>thr2);
        kractDura_het04(i) = length(tmp);
        kractDuraI_het04{i} = tmp';
    else
        kractDura_het04(i) = NaN;
        kractDuraI_het04{i} = [];
    end
end

kractDuraI_het05 = {};
for i = 1: size(krM_het05,2)
    if ismember(i,kract2nd_het05) == 1
        tmp = find(krM_het05(51:end,i)>thr2);
        kractDura_het05(i) = length(tmp);
        kractDuraI_het05{i} = tmp';
    else
        kractDura_het05(i) = NaN;
        kractDuraI_het05{i} = [];
    end
end

%% Transcription rate of active nuclei (addapted from new_strip_analysis_smoothed_TR_thr2 M2)

for i = 1:size(M_wt01,2)
    if ismember(i,act2nd_wt01) == 1
        TRB_wt01(i) =  trapz(M_wt01(actDuraI_wt01{i}+50,i))/actDura_wt01(i);
    else
        TRB_wt01(i) = NaN;
    end
end
for i = 1:size(M_wt02,2)
    if ismember(i,act2nd_wt02) == 1
        TRB_wt02(i) =  trapz(M_wt02(actDuraI_wt02{i}+50,i))/actDura_wt02(i);
    else
        TRB_wt02(i) = NaN;
    end
end
for i = 1:size(M_wt03,2)
    if ismember(i,act2nd_wt03) == 1
        TRB_wt03(i) =  trapz(M_wt03(actDuraI_wt03{i}+50,i))/actDura_wt03(i);
    else
        TRB_wt03(i) = NaN;
    end
end

for i = 1:size(krM_wt01,2)
    if ismember(i,kract2nd_wt01) == 1
        krTRB_wt01(i) =  trapz(krM_wt01(kractDuraI_wt01{i}+50,i))/kractDura_wt01(i);
    else
        krTRB_wt01(i) = NaN;
    end
end
for i = 1:size(krM_wt02,2)
    if ismember(i,kract2nd_wt02) == 1
        krTRB_wt02(i) =  trapz(krM_wt02(kractDuraI_wt02{i}+50,i))/kractDura_wt02(i);
    else
        krTRB_wt02(i) = NaN;
    end
end
for i = 1:size(krM_wt03,2)
    if ismember(i,kract2nd_wt03) == 1
        krTRB_wt03(i) =  trapz(krM_wt03(kractDuraI_wt03{i}+50,i))/kractDura_wt03(i);
    else
        krTRB_wt03(i) = NaN;
    end
end
for i = 1:size(krM_wt04,2)
    if ismember(i,kract2nd_wt04) == 1
        krTRB_wt04(i) =  trapz(krM_wt04(kractDuraI_wt04{i}+50,i))/kractDura_wt04(i);
    else
        krTRB_wt04(i) = NaN;
    end
end

for i = 1:size(krM_het01,2)
    if ismember(i,kract2nd_het01) == 1
        krTRB_het01(i) =  trapz(krM_het01(kractDuraI_het01{i}+50,i))/kractDura_het01(i);
    else
        krTRB_het01(i) = NaN;
    end
end
for i = 1:size(krM_het02,2)
    if ismember(i,kract2nd_het02) == 1
        krTRB_het02(i) =  trapz(krM_het02(kractDuraI_het02{i}+50,i))/kractDura_het02(i);
    else
        krTRB_het02(i) = NaN;
    end
end
for i = 1:size(krM_het03,2)
    if ismember(i,kract2nd_het03) == 1
        krTRB_het03(i) =  trapz(krM_het03(kractDuraI_het03{i}+50,i))/kractDura_het03(i);
    else
        krTRB_het03(i) = NaN;
    end
end
for i = 1:size(krM_het04,2)
    if ismember(i,kract2nd_het04) == 1
        krTRB_het04(i) =  trapz(krM_het04(kractDuraI_het04{i}+50,i))/kractDura_het04(i);
    else
        krTRB_het04(i) = NaN;
    end
end
for i = 1:size(krM_het05,2)
    if ismember(i,kract2nd_het05) == 1
        krTRB_het05(i) =  trapz(krM_het05(kractDuraI_het05{i}+50,i))/kractDura_het05(i);
    else
        krTRB_het05(i) = NaN;
    end
end

%% Combine embryos for each stripe
S1TRBArray_wt = [TRB_wt01(actS1second_wt01) TRB_wt02(actS1second_wt02) TRB_wt03(actS1second_wt03)...
    krTRB_wt01(kractS1second_wt01) krTRB_wt02(kractS1second_wt02) krTRB_wt03(kractS1second_wt03) krTRB_wt04(kractS1second_wt04)];
S1TRBArray_het = [krTRB_het01(kractS1second_het01) krTRB_het02(kractS1second_het02)...
     krTRB_het03(kractS1second_het03) krTRB_het04(kractS1second_het04) krTRB_het05(kractS1second_het05)];

 S2TRBArray_wt = [TRB_wt01(actS2second_wt01) TRB_wt02(actS2second_wt02) TRB_wt03(actS2second_wt03)...
    krTRB_wt01(kractS2second_wt01) krTRB_wt02(kractS2second_wt02) krTRB_wt03(kractS2second_wt03) krTRB_wt04(kractS2second_wt04)];
S2TRBArray_het = [krTRB_het01(kractS2second_het01) krTRB_het02(kractS2second_het02)...
     krTRB_het03(kractS2second_het03) krTRB_het04(kractS2second_het04) krTRB_het05(kractS2second_het05)];
 
 S3TRBArray_wt = [TRB_wt01(actS3second_wt01) TRB_wt02(actS3second_wt02) TRB_wt03(actS3second_wt03)...
    krTRB_wt01(kractS3second_wt01) krTRB_wt02(kractS3second_wt02) krTRB_wt03(kractS3second_wt03) krTRB_wt04(kractS3second_wt04)];
S3TRBArray_het = [krTRB_het01(kractS3second_het01) krTRB_het02(kractS3second_het02)...
     krTRB_het03(kractS3second_het03) krTRB_het04(kractS3second_het04) krTRB_het05(kractS3second_het05)];
 
 S4TRBArray_wt = [TRB_wt01(actS4second_wt01) TRB_wt02(actS4second_wt02) TRB_wt03(actS4second_wt03)...
    krTRB_wt01(kractS4second_wt01) krTRB_wt02(kractS4second_wt02) krTRB_wt03(kractS4second_wt03) krTRB_wt04(kractS4second_wt04)];
S4TRBArray_het = [krTRB_het01(kractS4second_het01) krTRB_het02(kractS4second_het02)...
     krTRB_het03(kractS4second_het03) krTRB_het04(kractS4second_het04) krTRB_het05(kractS4second_het05)];
 
 S5TRBArray_wt = [TRB_wt01(actS5second_wt01) TRB_wt02(actS5second_wt02) TRB_wt03(actS5second_wt03)...
    krTRB_wt01(kractS5second_wt01) krTRB_wt02(kractS5second_wt02) krTRB_wt03(kractS5second_wt03) krTRB_wt04(kractS5second_wt04)];
S5TRBArray_het = [krTRB_het01(kractS5second_het01) krTRB_het02(kractS5second_het02)...
     krTRB_het03(kractS5second_het03) krTRB_het04(kractS5second_het04) krTRB_het05(kractS5second_het05)];
 
 S6TRBArray_wt = [TRB_wt01(actS6second_wt01) TRB_wt02(actS6second_wt02) TRB_wt03(actS6second_wt03)...
    krTRB_wt01(kractS6second_wt01) krTRB_wt02(kractS6second_wt02) krTRB_wt03(kractS6second_wt03) krTRB_wt04(kractS6second_wt04)];
S6TRBArray_het = [krTRB_het01(kractS6second_het01) krTRB_het02(kractS6second_het02)...
     krTRB_het03(kractS6second_het03) krTRB_het04(kractS6second_het04) krTRB_het05(kractS6second_het05)];
%% histograms of M2 TR
binw =30/10000;

%M2
figure(371);
histogram(S1TRBArray_wt/10000,'BinWidth',binw,'Normalization','probability');hold on
% ylabel('Nuclei count');
% ax = gca;
% ax.YAxis.FontSize = 18;
% yyaxis right
histogram(S1TRBArray_het/10000,'BinWidth',binw,'Normalization','probability');
% xlabel('mRNA production rate (a.u./frame)');
% ylabel('Nuclei count');
legend('wt','Kr^{1}/+');
 hYLabel = get(gca,'YLabel');
 set(hYLabel,'rotation',270,'VerticalAlignment','middle');
hold off

figure(372);
histogram(S2TRBArray_wt/10000,'BinWidth',binw,'Normalization','probability');
% ylabel('Nuclei count');
hold on
% yyaxis right
histogram(S2TRBArray_het/10000,'BinWidth',binw,'Normalization','probability');
% xlabel('mRNA production rate (a.u./frame)');
% ylabel('Nuclei count');
legend('wt','Kr^{1}/+');
hYLabel = get(gca,'YLabel');
 set(hYLabel,'rotation',270,'VerticalAlignment','middle');
hold off

figure(373);
histogram(S3TRBArray_wt/10000,'BinWidth',binw,'Normalization','probability');
% ylabel('Nuclei count');
hold on
% yyaxis right
histogram(S3TRBArray_het/10000,'BinWidth',binw,'Normalization','probability');
% xlabel('mRNA production rate (a.u./frame)');
% ylabel('Nuclei count');
legend('wt','Kr^{1}/+');
hYLabel = get(gca,'YLabel');
 set(hYLabel,'rotation',270,'VerticalAlignment','middle');
hold off

figure(374);
histogram(S4TRBArray_wt/10000,'BinWidth',binw,'Normalization','probability');
% ylabel('Nuclei count');
hold on
% yyaxis right
histogram(S4TRBArray_het/10000,'BinWidth',binw,'Normalization','probability');
% xlabel('mRNA production rate (a.u./frame)');
% ylabel('Nuclei count');
legend('wt','Kr^{1}/+');
hYLabel = get(gca,'YLabel');
 set(hYLabel,'rotation',270,'VerticalAlignment','middle');
hold off

figure(375);
histogram(S5TRBArray_wt/10000,'BinWidth',binw,'Normalization','probability');
% ylabel('Nuclei count');
hold on
% yyaxis right
histogram(S5TRBArray_het/10000,'BinWidth',binw,'Normalization','probability');
% xlabel('mRNA production rate (a.u./frame)');
% ylabel('Nuclei count');
legend('wt','Kr^{1}/+');
hYLabel = get(gca,'YLabel');
 set(hYLabel,'rotation',270,'VerticalAlignment','middle');
hold off

figure(376);
histogram(S6TRBArray_wt/10000,'BinWidth',binw,'Normalization','probability');
% ylabel('Nuclei count');
hold on
% yyaxis right
histogram(S6TRBArray_het/10000,'BinWidth',binw,'Normalization','probability');
% xlabel('mRNA production rate (a.u./frame)');
% ylabel('Nuclei count');
hYLabel = get(gca,'YLabel');
 set(hYLabel,'rotation',270,'VerticalAlignment','middle');
% title('Stripe 1 mRNA production rate');
legend('wt','Kr^{1}/+');
hold off

%% boxplot active duration (not excluding bottom 10% act dura nuclei)
%S1
S1SecondAD_wt = [actDura_wt01(actS1second_wt01) actDura_wt02(actS1second_wt02) actDura_wt03(actS1second_wt03) ...
    kractDura_wt01(kractS1second_wt01) kractDura_wt02(kractS1second_wt02) kractDura_wt03(kractS1second_wt03) kractDura_wt04(kractS1second_wt04)];

S1SecondAD_het = [kractDura_het01(kractS1second_het01) kractDura_het02(kractS1second_het02) kractDura_het03(kractS1second_het03)...
    kractDura_het04(kractS1second_het04) kractDura_het05(kractS1second_het05)];


colors = [ 0.7, 0.3, 0.3; 0.5, 0.5, 0.5];
    figure(391);
    boxplot([S1SecondAD_wt, S1SecondAD_het]/100,[ones(1,length(S1SecondAD_wt)), 2*ones(1,length(S1SecondAD_het))],...
        'Labels',{'\it\bfwt','\it\bfKr^{1}/+'});
    ylim([0 0.6]);
%     title('Stripe 3 active duration');
    set(gca, 'TickLabelInterpreter', 'tex');
%     ylabel('active duration (frame)');
set(findobj(gca,'type','line'),'linew',2);
    h = findobj(gca,'Tag','Box'); 
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
%[h1,p1,c1,s1] = ttest2(S1SecondAD_wt,S1SecondAD_het)

%S2
S2SecondAD_wt = [actDura_wt01(actS2second_wt01) actDura_wt02(actS2second_wt02) actDura_wt03(actS2second_wt03) ...
    kractDura_wt01(kractS2second_wt01) kractDura_wt02(kractS2second_wt02) kractDura_wt03(kractS2second_wt03) kractDura_wt04(kractS2second_wt04)];

S2SecondAD_het = [kractDura_het01(kractS2second_het01) kractDura_het02(kractS2second_het02) kractDura_het03(kractS2second_het03)...
    kractDura_het04(kractS2second_het04) kractDura_het05(kractS2second_het05)];

figure(392);
boxplot([S2SecondAD_wt, S2SecondAD_het]/100,[ones(1,length(S2SecondAD_wt)), 2*ones(1,length(S2SecondAD_het))],...
    'Labels',{'\it\bfwt','\it\bfKr^{1}/+'});
ylim([0 0.6]);
% title('S2 active duration (active nuclei)');
% ylabel('active duration (frame)');
set(gca, 'TickLabelInterpreter', 'tex');
set(findobj(gca,'type','line'),'linew',2);
h = findobj(gca,'Tag','Box'); 
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

%[h2,p2,c2,s2] = ttest2(S2SecondAD_wt,S2SecondAD_het)

%S3
S3SecondAD_wt = [actDura_wt01(actS3second_wt01) actDura_wt02(actS3second_wt02) actDura_wt03(actS3second_wt03) ...
    kractDura_wt01(kractS3second_wt01) kractDura_wt02(kractS3second_wt02) kractDura_wt03(kractS3second_wt03) kractDura_wt04(kractS3second_wt04)];

S3SecondAD_het = [kractDura_het01(kractS3second_het01) kractDura_het02(kractS3second_het02) kractDura_het03(kractS3second_het03)...
    kractDura_het04(kractS3second_het04) kractDura_het05(kractS3second_het05)];


colors = [ 0.7, 0.3, 0.3; 0.5, 0.5, 0.5];
    figure(393);
    boxplot([S3SecondAD_wt, S3SecondAD_het]/100,[ones(1,length(S3SecondAD_wt)), 2*ones(1,length(S3SecondAD_het))],...
        'Labels',{'\it\bfwt','\it\bfKr^{1}/+'});
    ylim([0 0.6]);
%     title('Stripe 3 active duration');
    set(gca, 'TickLabelInterpreter', 'tex');
    set(findobj(gca,'type','line'),'linew',2);
%     ylabel('active duration (frame)');
    h = findobj(gca,'Tag','Box'); 
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

    % [h3,p3,c3,s3] = ttest2(S3SecondAD_wt,S3SecondAD_het)

%S4
S4SecondAD_wt = [actDura_wt01(actS4second_wt01) actDura_wt02(actS4second_wt02) actDura_wt03(actS4second_wt03) ...
    kractDura_wt01(kractS4second_wt01) kractDura_wt02(kractS4second_wt02) kractDura_wt03(kractS4second_wt03) kractDura_wt04(kractS4second_wt04)];

S4SecondAD_het = [kractDura_het01(kractS4second_het01) kractDura_het02(kractS4second_het02) kractDura_het03(kractS4second_het03)...
    kractDura_het04(kractS4second_het04) kractDura_het05(kractS4second_het05)];

figure(394);
boxplot([S4SecondAD_wt, S4SecondAD_het]/100,[ones(1,length(S4SecondAD_wt)), 2*ones(1,length(S4SecondAD_het))],...
    'Labels',{'\it\bfwt','\it\bfKr^{1}/+'});
ylim([0 0.6]);
set(gca, 'TickLabelInterpreter', 'tex');
%     title('Stripe 4 active duration');
%     ylabel('active duration (frame)');
set(findobj(gca,'type','line'),'linew',2);
h = findobj(gca,'Tag','Box'); 
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
%[h4,p4,c4,s4] = ttest2(S4SecondAD_wt,S4SecondAD_het)

%S5
S5SecondAD_wt = [actDura_wt01(actS5second_wt01) actDura_wt02(actS5second_wt02) actDura_wt03(actS5second_wt03) ...
    kractDura_wt01(kractS5second_wt01) kractDura_wt02(kractS5second_wt02) kractDura_wt03(kractS5second_wt03) kractDura_wt04(kractS5second_wt04)];

S5SecondAD_het = [kractDura_het01(kractS5second_het01) kractDura_het02(kractS5second_het02) kractDura_het03(kractS5second_het03)...
    kractDura_het04(kractS5second_het04) kractDura_het05(kractS5second_het05)];

figure(395);
boxplot([S5SecondAD_wt, S5SecondAD_het]/100,[ones(1,length(S5SecondAD_wt)), 2*ones(1,length(S5SecondAD_het))],...
    'Labels',{'\it\bfwt','\it\bfKr^{1}/+'});
ylim([0 0.6]);
% title('S5 active duration (active nuclei)');
% ylabel('active duration (frame)');
set(gca, 'TickLabelInterpreter', 'tex');
set(findobj(gca,'type','line'),'linew',2);
h = findobj(gca,'Tag','Box'); 
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
%[h5,p5,c5,s5] = ttest2(S5SecondAD_wt,S5SecondAD_het)

%S6
S6SecondAD_wt = [actDura_wt01(actS6second_wt01) actDura_wt02(actS6second_wt02) actDura_wt03(actS6second_wt03) ...
    kractDura_wt01(kractS6second_wt01) kractDura_wt02(kractS6second_wt02) kractDura_wt03(kractS6second_wt03) kractDura_wt04(kractS6second_wt04)];

S6SecondAD_het = [kractDura_het01(kractS6second_het01) kractDura_het02(kractS6second_het02) kractDura_het03(kractS6second_het03)...
    kractDura_het04(kractS6second_het04) kractDura_het05(kractS6second_het05)];

figure(396);
boxplot([S6SecondAD_wt, S6SecondAD_het]/100,[ones(1,length(S6SecondAD_wt)), 2*ones(1,length(S6SecondAD_het))],...
    'Labels',{'\it\bwt','\it\bfKr^{1}/+'});
ylim([0 0.6]);
% title('S6 active duration (active nuclei)');
% ylabel('active duration (frame)');
set(gca, 'TickLabelInterpreter', 'tex');
set(findobj(gca,'type','line'),'linew',2);
h = findobj(gca,'Tag','Box'); 
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
%[h6,p6,c6,s6] = ttest2(S6SecondAD_wt,S6SecondAD_het)

%% stripe 1 kinetics (fraction of act nuclei at given frame)
for i = 1:100
    tmp = find(M_wt01(i,stripe1_wt01)>= thr);
    stripe1kinetic_wt01(i)=length(tmp)/length(stripe1_wt01);
end
for i = 1:100
    tmp = find(M_wt01(i,stripe1_wt02)>= thr);
    stripe1kinetic_wt02(i)=length(tmp)/length(stripe1_wt02);
end
for i = 1:100
    tmp = find(M_wt01(i,stripe1_wt03)>= thr);
    stripe1kinetic_wt03(i)=length(tmp)/length(stripe1_wt03);
end

for i = 1:100
    tmp = find(krM_wt01(i,krstripe1_wt01)>= thr);
    krstripe1kinetic_wt01(i)=length(tmp)/length(krstripe1_wt01);
end
for i = 1:100
    tmp = find(krM_wt02(i,krstripe1_wt02)>= thr);
    krstripe1kinetic_wt02(i)=length(tmp)/length(krstripe1_wt02);
end
for i = 1:100
    tmp = find(krM_wt03(i,krstripe1_wt03)>= thr);
    krstripe1kinetic_wt03(i)=length(tmp)/length(krstripe1_wt03);
end
for i = 1:100
    tmp = find(krM_wt04(i,krstripe1_wt04)>= thr);
    krstripe1kinetic_wt04(i)=length(tmp)/length(krstripe1_wt04);
end

for i = 1:100
    tmp = find(krM_het01(i,krstripe1_het01)>= thr);
    krstripe1kinetic_het01(i)=length(tmp)/length(krstripe1_het01);
end
for i = 1:100
    tmp = find(krM_het02(i,krstripe1_het02)>= thr);
    krstripe1kinetic_het02(i)=length(tmp)/length(krstripe1_het02);
end
for i = 1:100
    tmp = find(krM_het03(i,krstripe1_het03)>= thr);
    krstripe1kinetic_het03(i)=length(tmp)/length(krstripe1_het03);
end
for i = 1:100
    tmp = find(krM_het04(i,krstripe1_het04)>= thr);
    krstripe1kinetic_het04(i)=length(tmp)/length(krstripe1_het04);
end
for i = 1:100
    tmp = find(krM_het05(i,krstripe1_het05)>= thr);
    krstripe1kinetic_het05(i)=length(tmp)/length(krstripe1_het05);
end

for i = 1:100
   stripe1ActFrac_wt(i) = nanmean([stripe1kinetic_wt01(i) stripe1kinetic_wt02(i) stripe1kinetic_wt03(i)...
       krstripe1kinetic_wt01(i) krstripe1kinetic_wt02(i) krstripe1kinetic_wt03(i) krstripe1kinetic_wt04(i)]);
   stripe1ActFrac_het(i) = nanmean([krstripe1kinetic_het01(i) krstripe1kinetic_het02(i) krstripe1kinetic_het03(i)...
    krstripe1kinetic_het04(i) krstripe1kinetic_het05(i)]);
    stripe1ActFracSE_wt(i) = nanstd([stripe1kinetic_wt01(i) stripe1kinetic_wt02(i) stripe1kinetic_wt03(i)...
       krstripe1kinetic_wt01(i) krstripe1kinetic_wt02(i) krstripe1kinetic_wt03(i) krstripe1kinetic_wt04(i)])/sqrt(7);
   stripe1ActFracSE_het(i) = nanstd([krstripe1kinetic_het01(i) krstripe1kinetic_het02(i) krstripe1kinetic_het03(i)...
    krstripe1kinetic_het04(i) krstripe1kinetic_het05(i)])/sqrt(5);
end

figure(251);
plot(0.01:0.01:1,stripe1ActFrac_wt,'-k','LineWidth',3);hold on
plot(0.01:0.01:1,stripe1ActFrac_het,'-r','LineWidth',3);
shadedErrorBar2(0.01:0.01:1,stripe1ActFrac_wt,stripe1ActFracSE_wt,'lineprops', {'-k','LineWidth',3}); 
shadedErrorBar2(0.01:0.01:1,stripe1ActFrac_het,stripe1ActFracSE_het,'lineprops', {'-r','LineWidth',3});
% title('Average kinetics stripe 1');
legend('wt','Kr^{1}/+');
% ylabel('Fraction of active nuclei');xlabel('nc 14 (normalized)');
hold off
%% stripe 2 kinetics (fraction of act nuclei at given frame)
for i = 1:100
    tmp = find(M_wt01(i,stripe2_wt01)>= thr);
    stripe2kinetic_wt01(i)=length(tmp)/length(stripe2_wt01);
end
for i = 1:100
    tmp = find(M_wt01(i,stripe2_wt02)>= thr);
    stripe2kinetic_wt02(i)=length(tmp)/length(stripe2_wt02);
end
for i = 1:100
    tmp = find(M_wt01(i,stripe2_wt03)>= thr);
    stripe2kinetic_wt03(i)=length(tmp)/length(stripe2_wt03);
end

for i = 1:100
    tmp = find(krM_wt01(i,krstripe2_wt01)>= thr);
    krstripe2kinetic_wt01(i)=length(tmp)/length(krstripe2_wt01);
end
for i = 1:100
    tmp = find(krM_wt02(i,krstripe2_wt02)>= thr);
    krstripe2kinetic_wt02(i)=length(tmp)/length(krstripe2_wt02);
end
for i = 1:100
    tmp = find(krM_wt03(i,krstripe2_wt03)>= thr);
    krstripe2kinetic_wt03(i)=length(tmp)/length(krstripe2_wt03);
end
for i = 1:100
    tmp = find(krM_wt04(i,krstripe2_wt04)>= thr);
    krstripe2kinetic_wt04(i)=length(tmp)/length(krstripe2_wt04);
end

for i = 1:100
    tmp = find(krM_het01(i,krstripe2_het01)>= thr);
    krstripe2kinetic_het01(i)=length(tmp)/length(krstripe2_het01);
end
for i = 1:100
    tmp = find(krM_het02(i,krstripe2_het02)>= thr);
    krstripe2kinetic_het02(i)=length(tmp)/length(krstripe2_het02);
end
for i = 1:100
    tmp = find(krM_het03(i,krstripe2_het03)>= thr);
    krstripe2kinetic_het03(i)=length(tmp)/length(krstripe2_het03);
end
for i = 1:100
    tmp = find(krM_het04(i,krstripe2_het04)>= thr);
    krstripe2kinetic_het04(i)=length(tmp)/length(krstripe2_het04);
end
for i = 1:100
    tmp = find(krM_het05(i,krstripe2_het05)>= thr);
    krstripe2kinetic_het05(i)=length(tmp)/length(krstripe2_het05);
end

for i = 1:100
   stripe2ActFrac_wt(i) = nanmean([stripe2kinetic_wt01(i) stripe2kinetic_wt02(i) stripe2kinetic_wt03(i)...
       krstripe2kinetic_wt01(i) krstripe2kinetic_wt02(i) krstripe2kinetic_wt03(i) krstripe2kinetic_wt04(i)]);
   stripe2ActFrac_het(i) = nanmean([krstripe2kinetic_het01(i) krstripe2kinetic_het02(i) krstripe2kinetic_het03(i)...
    krstripe2kinetic_het04(i) krstripe2kinetic_het05(i)]);
    stripe2ActFracSE_wt(i) = nanstd([stripe2kinetic_wt01(i) stripe2kinetic_wt02(i) stripe2kinetic_wt03(i)...
       krstripe2kinetic_wt01(i) krstripe2kinetic_wt02(i) krstripe2kinetic_wt03(i) krstripe2kinetic_wt04(i)])/sqrt(7);
   stripe2ActFracSE_het(i) = nanstd([krstripe2kinetic_het01(i) krstripe2kinetic_het02(i) krstripe2kinetic_het03(i)...
    krstripe2kinetic_het04(i) krstripe2kinetic_het05(i)])/sqrt(5);
end

figure(252);
plot(0.01:0.01:1,stripe2ActFrac_wt,'-k','LineWidth',3);hold on
plot(0.01:0.01:1,stripe2ActFrac_het,'-r','LineWidth',3);
shadedErrorBar2(0.01:0.01:1,stripe2ActFrac_wt,stripe2ActFracSE_wt,'lineprops', {'-k','LineWidth',3}); 
shadedErrorBar2(0.01:0.01:1,stripe2ActFrac_het,stripe2ActFracSE_het,'lineprops', {'-r','LineWidth',3});
% title('Average kinetics stripe 2');
legend('wt','Kr^{1}/+');
% ylabel('Fraction of active nuclei');xlabel('nc 14 (normalized)');
ax = gca;
ax.YAxis.FontSize = 18;
ax.XAxis.FontSize = 16;
ax.Legend.FontSize = 18;
ax.Title.FontSize = 18;
hold off

%% stripe 3 kinetics (fraction of act nuclei at given frame)
for i = 1:100
    tmp = find(M_wt01(i,stripe3_wt01)>= thr);
    stripe3kinetic_wt01(i)=length(tmp)/length(stripe3_wt01);
end
for i = 1:100
    tmp = find(M_wt01(i,stripe3_wt02)>= thr);
    stripe3kinetic_wt02(i)=length(tmp)/length(stripe3_wt02);
end
for i = 1:100
    tmp = find(M_wt01(i,stripe3_wt03)>= thr);
    stripe3kinetic_wt03(i)=length(tmp)/length(stripe3_wt03);
end

for i = 1:100
    tmp = find(krM_wt01(i,krstripe3_wt01)>= thr);
    krstripe3kinetic_wt01(i)=length(tmp)/length(krstripe3_wt01);
end
for i = 1:100
    tmp = find(krM_wt02(i,krstripe3_wt02)>= thr);
    krstripe3kinetic_wt02(i)=length(tmp)/length(krstripe3_wt02);
end
for i = 1:100
    tmp = find(krM_wt03(i,krstripe3_wt03)>= thr);
    krstripe3kinetic_wt03(i)=length(tmp)/length(krstripe3_wt03);
end
for i = 1:100
    tmp = find(krM_wt04(i,krstripe3_wt04)>= thr);
    krstripe3kinetic_wt04(i)=length(tmp)/length(krstripe3_wt04);
end

for i = 1:100
    tmp = find(krM_het01(i,krstripe3_het01)>= thr);
    krstripe3kinetic_het01(i)=length(tmp)/length(krstripe3_het01);
end
for i = 1:100
    tmp = find(krM_het02(i,krstripe3_het02)>= thr);
    krstripe3kinetic_het02(i)=length(tmp)/length(krstripe3_het02);
end
for i = 1:100
    tmp = find(krM_het03(i,krstripe3_het03)>= thr);
    krstripe3kinetic_het03(i)=length(tmp)/length(krstripe3_het03);
end
for i = 1:100
    tmp = find(krM_het04(i,krstripe3_het04)>= thr);
    krstripe3kinetic_het04(i)=length(tmp)/length(krstripe3_het04);
end
for i = 1:100
    tmp = find(krM_het05(i,krstripe3_het05)>= thr);
    krstripe3kinetic_het05(i)=length(tmp)/length(krstripe3_het05);
end

for i = 1:100
   stripe3ActFrac_wt(i) = nanmean([stripe3kinetic_wt01(i) stripe3kinetic_wt02(i) stripe3kinetic_wt03(i)...
       krstripe3kinetic_wt01(i) krstripe3kinetic_wt02(i) krstripe3kinetic_wt03(i) krstripe3kinetic_wt04(i)]);
   stripe3ActFrac_het(i) = nanmean([krstripe3kinetic_het01(i) krstripe3kinetic_het02(i) krstripe3kinetic_het03(i)...
    krstripe3kinetic_het04(i) krstripe3kinetic_het05(i)]);
    stripe3ActFracSE_wt(i) = nanstd([stripe3kinetic_wt01(i) stripe3kinetic_wt02(i) stripe3kinetic_wt03(i)...
       krstripe3kinetic_wt01(i) krstripe3kinetic_wt02(i) krstripe3kinetic_wt03(i) krstripe3kinetic_wt04(i)])/sqrt(7);
   stripe3ActFracSE_het(i) = nanstd([krstripe3kinetic_het01(i) krstripe3kinetic_het02(i) krstripe3kinetic_het03(i)...
    krstripe3kinetic_het04(i) krstripe3kinetic_het05(i)])/sqrt(5);
end

figure(253);
plot(0.01:0.01:1,stripe3ActFrac_wt,'-k','LineWidth',3);hold on
plot(0.01:0.01:1,stripe3ActFrac_het,'-r','LineWidth',3);
shadedErrorBar2(0.01:0.01:1,stripe3ActFrac_wt,stripe3ActFracSE_wt,'lineprops', {'-k','LineWidth',3}); 
shadedErrorBar2(0.01:0.01:1,stripe3ActFrac_het,stripe3ActFracSE_het,'lineprops', {'-r','LineWidth',3});
legend('wt','Kr^{1}/+');
% title('Average kinetics stripe 3');
% legend('wt','Kr^{1}/+','Location','northwest');
% ylabel('Fraction of active nuclei');xlabel('nc 14 (normalized)');
ax = gca;
ax.YAxis.FontSize = 18;
ax.XAxis.FontSize = 16;
ax.Legend.FontSize = 18;
ax.Title.FontSize = 18;
hold off

%% stripe 4 kinetics (fraction of act nuclei at given frame)
for i = 1:100
    tmp = find(M_wt01(i,stripe4_wt01)>= thr);
    stripe4kinetic_wt01(i)=length(tmp)/length(stripe4_wt01);
end
for i = 1:100
    tmp = find(M_wt01(i,stripe4_wt02)>= thr);
    stripe4kinetic_wt02(i)=length(tmp)/length(stripe4_wt02);
end
for i = 1:100
    tmp = find(M_wt01(i,stripe4_wt03)>= thr);
    stripe4kinetic_wt03(i)=length(tmp)/length(stripe4_wt03);
end

for i = 1:100
    tmp = find(krM_wt01(i,krstripe4_wt01)>= thr);
    krstripe4kinetic_wt01(i)=length(tmp)/length(krstripe4_wt01);
end
for i = 1:100
    tmp = find(krM_wt02(i,krstripe4_wt02)>= thr);
    krstripe4kinetic_wt02(i)=length(tmp)/length(krstripe4_wt02);
end
for i = 1:100
    tmp = find(krM_wt03(i,krstripe4_wt03)>= thr);
    krstripe4kinetic_wt03(i)=length(tmp)/length(krstripe4_wt03);
end
for i = 1:100
    tmp = find(krM_wt04(i,krstripe4_wt04)>= thr);
    krstripe4kinetic_wt04(i)=length(tmp)/length(krstripe4_wt04);
end

for i = 1:100
    tmp = find(krM_het01(i,krstripe4_het01)>= thr);
    krstripe4kinetic_het01(i)=length(tmp)/length(krstripe4_het01);
end
for i = 1:100
    tmp = find(krM_het02(i,krstripe4_het02)>= thr);
    krstripe4kinetic_het02(i)=length(tmp)/length(krstripe4_het02);
end
for i = 1:100
    tmp = find(krM_het03(i,krstripe4_het03)>= thr);
    krstripe4kinetic_het03(i)=length(tmp)/length(krstripe4_het03);
end
for i = 1:100
    tmp = find(krM_het04(i,krstripe4_het04)>= thr);
    krstripe4kinetic_het04(i)=length(tmp)/length(krstripe4_het04);
end
for i = 1:100
    tmp = find(krM_het05(i,krstripe4_het05)>= thr);
    krstripe4kinetic_het05(i)=length(tmp)/length(krstripe4_het05);
end

for i = 1:100
   stripe4ActFrac_wt(i) = nanmean([stripe4kinetic_wt01(i) stripe4kinetic_wt02(i) stripe4kinetic_wt03(i)...
       krstripe4kinetic_wt01(i) krstripe4kinetic_wt02(i) krstripe4kinetic_wt03(i) krstripe4kinetic_wt04(i)]);
   stripe4ActFrac_het(i) = nanmean([krstripe4kinetic_het01(i) krstripe4kinetic_het02(i) krstripe4kinetic_het03(i)...
    krstripe4kinetic_het04(i) krstripe4kinetic_het05(i)]);
    stripe4ActFracSE_wt(i) = nanstd([stripe4kinetic_wt01(i) stripe4kinetic_wt02(i) stripe4kinetic_wt03(i)...
       krstripe4kinetic_wt01(i) krstripe4kinetic_wt02(i) krstripe4kinetic_wt03(i) krstripe4kinetic_wt04(i)])/sqrt(7);
   stripe4ActFracSE_het(i) = nanstd([krstripe4kinetic_het01(i) krstripe4kinetic_het02(i) krstripe4kinetic_het03(i)...
    krstripe4kinetic_het04(i) krstripe4kinetic_het05(i)])/sqrt(5);
end

figure(254);
plot(0.01:0.01:1,stripe4ActFrac_wt,'-k','LineWidth',3);hold on
plot(0.01:0.01:1,stripe4ActFrac_het,'-r','LineWidth',3);
shadedErrorBar2(0.01:0.01:1,stripe4ActFrac_wt,stripe4ActFracSE_wt,'lineprops', {'-k','LineWidth',3}); 
shadedErrorBar2(0.01:0.01:1,stripe4ActFrac_het,stripe4ActFracSE_het,'lineprops', {'-r','LineWidth',3});
legend('wt','Kr^{1}/+');
% title('Average kinetics stripe 4');
% legend('wt','Kr^{1}/+','Location','northwest');
% ylabel('Fraction of active nuclei');xlabel('nc 14 (normalized)');
ax = gca;
ax.YAxis.FontSize = 18;
ax.XAxis.FontSize = 16;
ax.Legend.FontSize = 18;
ax.Title.FontSize = 18;
hold off

%% stripe 5 kinetics (fraction of act nuclei at given frame)
for i = 1:100
    tmp = find(M_wt01(i,stripe5_wt01)>= thr);
    stripe5kinetic_wt01(i)=length(tmp)/length(stripe5_wt01);
end
for i = 1:100
    tmp = find(M_wt01(i,stripe5_wt02)>= thr);
    stripe5kinetic_wt02(i)=length(tmp)/length(stripe5_wt02);
end
for i = 1:100
    tmp = find(M_wt01(i,stripe5_wt03)>= thr);
    stripe5kinetic_wt03(i)=length(tmp)/length(stripe5_wt03);
end

for i = 1:100
    tmp = find(krM_wt01(i,krstripe5_wt01)>= thr);
    krstripe5kinetic_wt01(i)=length(tmp)/length(krstripe5_wt01);
end
for i = 1:100
    tmp = find(krM_wt02(i,krstripe5_wt02)>= thr);
    krstripe5kinetic_wt02(i)=length(tmp)/length(krstripe5_wt02);
end
for i = 1:100
    tmp = find(krM_wt03(i,krstripe5_wt03)>= thr);
    krstripe5kinetic_wt03(i)=length(tmp)/length(krstripe5_wt03);
end
for i = 1:100
    tmp = find(krM_wt04(i,krstripe5_wt04)>= thr);
    krstripe5kinetic_wt04(i)=length(tmp)/length(krstripe5_wt04);
end

for i = 1:100
    tmp = find(krM_het01(i,krstripe5_het01)>= thr);
    krstripe5kinetic_het01(i)=length(tmp)/length(krstripe5_het01);
end
for i = 1:100
    tmp = find(krM_het02(i,krstripe5_het02)>= thr);
    krstripe5kinetic_het02(i)=length(tmp)/length(krstripe5_het02);
end
for i = 1:100
    tmp = find(krM_het03(i,krstripe5_het03)>= thr);
    krstripe5kinetic_het03(i)=length(tmp)/length(krstripe5_het03);
end
for i = 1:100
    tmp = find(krM_het04(i,krstripe5_het04)>= thr);
    krstripe5kinetic_het04(i)=length(tmp)/length(krstripe5_het04);
end
for i = 1:100
    tmp = find(krM_het05(i,krstripe5_het05)>= thr);
    krstripe5kinetic_het05(i)=length(tmp)/length(krstripe5_het05);
end

for i = 1:100
   stripe5ActFrac_wt(i) = nanmean([stripe5kinetic_wt01(i) stripe5kinetic_wt02(i) stripe5kinetic_wt03(i)...
       krstripe5kinetic_wt01(i) krstripe5kinetic_wt02(i) krstripe5kinetic_wt03(i) krstripe5kinetic_wt04(i)]);
   stripe5ActFrac_het(i) = nanmean([krstripe5kinetic_het01(i) krstripe5kinetic_het02(i) krstripe5kinetic_het03(i)...
    krstripe5kinetic_het04(i) krstripe5kinetic_het05(i)]);
    stripe5ActFracSE_wt(i) = nanstd([stripe5kinetic_wt01(i) stripe5kinetic_wt02(i) stripe5kinetic_wt03(i)...
       krstripe5kinetic_wt01(i) krstripe5kinetic_wt02(i) krstripe5kinetic_wt03(i) krstripe5kinetic_wt04(i)])/sqrt(7);
   stripe5ActFracSE_het(i) = nanstd([krstripe5kinetic_het01(i) krstripe5kinetic_het02(i) krstripe5kinetic_het03(i)...
    krstripe5kinetic_het04(i) krstripe5kinetic_het05(i)])/sqrt(5);
end

figure(255);
plot(0.01:0.01:1,stripe5ActFrac_wt,'-k','LineWidth',3);hold on
plot(0.01:0.01:1,stripe5ActFrac_het,'-r','LineWidth',3);
shadedErrorBar2(0.01:0.01:1,stripe5ActFrac_wt,stripe5ActFracSE_wt,'lineprops', {'-k','LineWidth',3}); 
shadedErrorBar2(0.01:0.01:1,stripe5ActFrac_het,stripe5ActFracSE_het,'lineprops', {'-r','LineWidth',3});
legend('wt','Kr^{1}/+');
% title('Average kinetics stripe 5');
% legend('wt','Kr^{1}/+','Location','northwest');
% ylabel('Fraction of active nuclei');xlabel('nc 14 (normalized)');
ax = gca;
ax.YAxis.FontSize = 18;
ax.XAxis.FontSize = 16;
ax.Legend.FontSize = 18;
ax.Title.FontSize = 18;
hold off
%% stripe 6 kinetics (fraction of act nuclei at given frame)
for i = 1:100
    tmp = find(M_wt01(i,stripe6_wt01)>= thr);
    stripe6kinetic_wt01(i)=length(tmp)/length(stripe6_wt01);
end
for i = 1:100
    tmp = find(M_wt01(i,stripe6_wt02)>= thr);
    stripe6kinetic_wt02(i)=length(tmp)/length(stripe6_wt02);
end
for i = 1:100
    tmp = find(M_wt01(i,stripe6_wt03)>= thr);
    stripe6kinetic_wt03(i)=length(tmp)/length(stripe6_wt03);
end

for i = 1:100
    tmp = find(krM_wt01(i,krstripe6_wt01)>= thr);
    krstripe6kinetic_wt01(i)=length(tmp)/length(krstripe6_wt01);
end
for i = 1:100
    tmp = find(krM_wt02(i,krstripe6_wt02)>= thr);
    krstripe6kinetic_wt02(i)=length(tmp)/length(krstripe6_wt02);
end
for i = 1:100
    tmp = find(krM_wt03(i,krstripe6_wt03)>= thr);
    krstripe6kinetic_wt03(i)=length(tmp)/length(krstripe6_wt03);
end
for i = 1:100
    tmp = find(krM_wt04(i,krstripe6_wt04)>= thr);
    krstripe6kinetic_wt04(i)=length(tmp)/length(krstripe6_wt04);
end

for i = 1:100
    tmp = find(krM_het01(i,krstripe6_het01)>= thr);
    krstripe6kinetic_het01(i)=length(tmp)/length(krstripe6_het01);
end
for i = 1:100
    tmp = find(krM_het02(i,krstripe6_het02)>= thr);
    krstripe6kinetic_het02(i)=length(tmp)/length(krstripe6_het02);
end
for i = 1:100
    tmp = find(krM_het03(i,krstripe6_het03)>= thr);
    krstripe6kinetic_het03(i)=length(tmp)/length(krstripe6_het03);
end
for i = 1:100
    tmp = find(krM_het04(i,krstripe6_het04)>= thr);
    krstripe6kinetic_het04(i)=length(tmp)/length(krstripe6_het04);
end
for i = 1:100
    tmp = find(krM_het05(i,krstripe6_het05)>= thr);
    krstripe6kinetic_het05(i)=length(tmp)/length(krstripe6_het05);
end

for i = 1:100
   stripe6ActFrac_wt(i) = nanmean([stripe6kinetic_wt01(i) stripe6kinetic_wt02(i) stripe6kinetic_wt03(i)...
       krstripe6kinetic_wt01(i) krstripe6kinetic_wt02(i) krstripe6kinetic_wt03(i) krstripe6kinetic_wt04(i)]);
   stripe6ActFrac_het(i) = nanmean([krstripe6kinetic_het01(i) krstripe6kinetic_het02(i) krstripe6kinetic_het03(i)...
    krstripe6kinetic_het04(i) krstripe6kinetic_het05(i)]);
    stripe6ActFracSE_wt(i) = nanstd([stripe6kinetic_wt01(i) stripe6kinetic_wt02(i) stripe6kinetic_wt03(i)...
       krstripe6kinetic_wt01(i) krstripe6kinetic_wt02(i) krstripe6kinetic_wt03(i) krstripe6kinetic_wt04(i)])/sqrt(7);
   stripe6ActFracSE_het(i) = nanstd([krstripe6kinetic_het01(i) krstripe6kinetic_het02(i) krstripe6kinetic_het03(i)...
    krstripe6kinetic_het04(i) krstripe6kinetic_het05(i)])/sqrt(5);
end

figure(256);
plot(0.01:0.01:1,stripe6ActFrac_wt,'-k','LineWidth',3);hold on
plot(0.01:0.01:1,stripe6ActFrac_het,'-r','LineWidth',3);
shadedErrorBar2(0.01:0.01:1,stripe6ActFrac_wt,stripe6ActFracSE_wt,'lineprops', {'-k','LineWidth',3}); 
shadedErrorBar2(0.01:0.01:1,stripe6ActFrac_het,stripe6ActFracSE_het,'lineprops', {'-r','LineWidth',3});
% title('Average kinetics stripe 6');
legend('wt','Kr^{1}/+','Location','northwest');
% ylabel('Fraction of active nuclei');xlabel('nc 14 (normalized)');
ax = gca;
ax.YAxis.FontSize = 18;
ax.XAxis.FontSize = 16;
ax.Legend.FontSize = 18;
ax.Title.FontSize = 18;
hold off

for i = 1:100
    for j = 1:bincount
        binAveAmp_wt01{i,j} = mean(M_wt01(i,binEL_wt01{j}));
        binAveAmp_wt02{i,j} = mean(M_wt02(i,binEL_wt02{j}));
        binAveAmp_wt03{i,j} = mean(M_wt03(i,binEL_wt03{j}));
        
        krbinAveAmp_wt01{i,j} = mean(krM_wt01(i,krbinEL_wt01{j}));
        krbinAveAmp_wt02{i,j} = mean(krM_wt02(i,krbinEL_wt02{j}));
        krbinAveAmp_wt03{i,j} = mean(krM_wt03(i,krbinEL_wt03{j}));
        krbinAveAmp_wt04{i,j} = mean(krM_wt04(i,krbinEL_wt04{j}));
        
        krbinAveAmp_het01{i,j} = mean(krM_het01(i,krbinEL_het01{j}));
        krbinAveAmp_het02{i,j} = mean(krM_het02(i,krbinEL_het02{j}));
        krbinAveAmp_het03{i,j} = mean(krM_het03(i,krbinEL_het03{j}));
        krbinAveAmp_het04{i,j} = mean(krM_het04(i,krbinEL_het04{j}));
        krbinAveAmp_het05{i,j} = mean(krM_het05(i,krbinEL_het05{j}));
    end
end

figure(1273);
    imagesc(cell2mat(krbinAveAmp_wt04)/1000);colorbar;
    caxis([0, 1.2]);
    xlim([9 47]);
    colormap(gca, myCustomColormap);
hold off
    
figure(1274);
    imagesc(cell2mat(krbinAveAmp_het05)/1000);colorbar;
    caxis([0, 1.2]);    
    xlim([9 47]);
    ylim([0 100]);
    colormap(gca, myCustomColormap);
    xticklabels({'0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9'});
hold off

%% mRNA 2nd half NC14 - EL
period = 51:100;

for i = 1:50
    binAveEL_wt01(i) = nanmean(trapz(M_wt01(period,binEL_wt01{i})));
end
for i = 1:50
    binAveEL_wt02(i) = nanmean(trapz(M_wt02(period,binEL_wt02{i})));
end
for i = 1:50
    binAveEL_wt03(i) = nanmean(trapz(M_wt03(period,binEL_wt03{i})));
end
for i = 1:50
    krbinAveEL_wt01(i) = nanmean(trapz(krM_wt01(period,krbinEL_wt01{i})));
end
for i = 1:50
    krbinAveEL_wt02(i) = nanmean(trapz(krM_wt02(period,krbinEL_wt02{i})));
end
for i = 1:50
    krbinAveEL_wt03(i) = nanmean(trapz(krM_wt03(period,krbinEL_wt03{i})));
end
for i = 1:50
    krbinAveEL_wt04(i) = nanmean(trapz(krM_wt04(period,krbinEL_wt04{i})));
end
for i = 1:50
    krbinAveEL_het01(i) = nanmean(trapz(krM_het01(period,krbinEL_het01{i})));
end
for i = 1:50
    krbinAveEL_het02(i) = nanmean(trapz(krM_het02(period,krbinEL_het02{i})));
end
for i = 1:50
    krbinAveEL_het03(i) = nanmean(trapz(krM_het03(period,krbinEL_het03{i})));
end
for i = 1:50
    krbinAveEL_het04(i) = nanmean(trapz(krM_het04(period,krbinEL_het04{i})));
end
for i = 1:50
    krbinAveEL_het05(i) = nanmean(trapz(krM_het05(period,krbinEL_het05{i})));
end
%avereage embryo
for i = 1:bincount
    binAveEL_wt(i)= mean([binAveEL_wt01(i), binAveEL_wt02(i), binAveEL_wt03(i),...
        krbinAveEL_wt01(i), krbinAveEL_wt02(i), krbinAveEL_wt03(i), krbinAveEL_wt04(i)]);
    binAveEL_het(i)= mean([krbinAveEL_het01(i), krbinAveEL_het02(i), krbinAveEL_het03(i),...
        krbinAveEL_het04(i), krbinAveEL_het05(i)]);
    
    binAveELSTD_wt(i)= std([binAveEL_wt01(i), binAveEL_wt02(i), binAveEL_wt03(i),...
        krbinAveEL_wt01(i), krbinAveEL_wt02(i), krbinAveEL_wt03(i), krbinAveEL_wt04(i)])/sqrt(7);
    binAveELSTD_het(i)= std([krbinAveEL_het01(i), krbinAveEL_het02(i), krbinAveEL_het03(i),...
        krbinAveEL_het04(i), krbinAveEL_het05(i)])/sqrt(5);
end

%figures
figure(1284);
plot(2:2:100,binAveEL_wt/1000,'k','LineWidth',3);hold on
plot(2:2:100,binAveEL_het/1000,'r','LineWidth',3);hold on
shadedErrorBar2(2:2:100,binAveEL_wt/1000,binAveELSTD_wt/1000,'lineprops', {'-k','LineWidth',3}); hold on
shadedErrorBar2(2:2:100,binAveEL_het/1000,binAveELSTD_het/1000,'lineprops', {'-r','LineWidth',3});
% title('Total mRNA production (2nd half nc14)');
legend('wt','Kr^{1}/+');
xlim([20 90]);
% ylabel('mRNA production (a.u.)'); 
% xlabel('anterior \rightarrow posterior');
ax = gca;
ax.YAxis.FontSize = 18;
ax.XAxis.FontSize = 16;
ax.Legend.FontSize = 18;
ax.Title.FontSize = 18;
hold off

%% calculate coefficient of varaition
for i = 1:50
    binAveELCV_wt(i) = 100*binAveELSTD_wt(i)/binAveEL_wt(i);
    binAveELCV_het(i) = 100*binAveELSTD_het(i)/binAveEL_het(i);
end

%% find stripe peaks and bases of 2nd half mRNA output
peakM_wt = [];%stripe peak of 2nd half mRNA output
peakM_het = [];
baseM_wt = [];%base between two stripe of 2nd half mRNA output
baseM_het = [];
for i = 10:45
    if binAveEL_wt(i-1) < binAveEL_wt(i) && binAveEL_wt(i+1) < binAveEL_wt(i)
        peakM_wt(end+1) = i;
    end
    if binAveEL_het(i-1) < binAveEL_het(i) && binAveEL_het(i+1) < binAveEL_het(i)
        peakM_het(end+1) = i;
    end
    
    %base
    if binAveEL_wt(i-1) > binAveEL_wt(i) && binAveEL_wt(i+1) > binAveEL_wt(i)
        baseM_wt(end+1) = i;
    end
    if binAveEL_het(i-1) > binAveEL_het(i) && binAveEL_het(i+1) > binAveEL_het(i)
        baseM_het(end+1) = i;
    end
end

boundPercent = 0.5;
%define base between stripes (noise)
baseMAve_wt = mean(binAveEL_wt(baseM_wt));
baseMAve_het = mean(binAveEL_het(baseM_het));

for i = 1:7
    tmpxABS=[];
    tmpxS=[];
    tmpx=[];
    tmpy=[];

    peakBSA_wt(i,1:50) = 15000;% Boundary signal array
    [tmpx,tmpy] = intersections(1:50,binAveEL_wt(:),1:50,peakBSA_wt(i,:));
    for j = 1:length(tmpx)
        tmpxABS(j)=abs(tmpx(j)-peakM_wt(i));
    end
    [tmpxS tmpxSind]= sort(tmpxABS);
   peakMx0_wt(i,1:2)= tmpx(tmpxSind(1:2));
end

for i = 1:7
    tmpxABS=[];
    tmpxS=[];
    tmpx=[];
    tmpy=[];
    
    if i == 4
        peakBSA_het(i,1:50) = 11100;
    else
        peakBSA_het(i,1:50) = 15000;% Boundary signal array
    end
    [tmpx,tmpy] = intersections(1:50,binAveEL_het(:),1:50,peakBSA_het(i,:));
    for j = 1:length(tmpx)
        tmpxABS(j)=abs(tmpx(j)-peakM_het(i));
    end
    [tmpxS tmpxSind]= sort(tmpxABS);
    peakMx0_het(i,1:2)= tmpx(tmpxSind(1:2));
end
%% variance of mRNA production along x axis
for i = 1:bincount
    binAveELVar_wt(i)= var([binAveEL_wt01(i), binAveEL_wt02(i), binAveEL_wt03(i),...
        binAveEL_wt01(i), krbinAveEL_wt02(i), krbinAveEL_wt03(i), krbinAveEL_wt04(i)]);
    binAveELVar_het(i)= var([krbinAveEL_het01(i), krbinAveEL_het02(i), krbinAveEL_het03(i),...
        krbinAveEL_het04(i), krbinAveEL_het05(i)]);
end

yMax = 150000000;
yMax2 = 30000;
v_wt = [peakMx0_wt(1,1)*0.02 0; peakMx0_wt(1,2)*0.02 0; peakMx0_wt(1,2)*0.02 yMax2; peakMx0_wt(1,1)*0.02 yMax2];
v2_wt = [peakMx0_wt(2,1)*0.02 0; peakMx0_wt(2,2)*0.02 0; peakMx0_wt(2,2)*0.02 yMax2; peakMx0_wt(2,1)*0.02 yMax2];
v3_wt = [peakMx0_wt(3,1)*0.02 0; peakMx0_wt(3,2)*0.02 0; peakMx0_wt(3,2)*0.02 yMax2; peakMx0_wt(3,1)*0.02 yMax2];
v4_wt = [peakMx0_wt(4,1)*0.02 0; peakMx0_wt(4,2)*0.02 0; peakMx0_wt(4,2)*0.02 yMax2; peakMx0_wt(4,1)*0.02 yMax2];
v5_wt = [peakMx0_wt(5,1)*0.02 0; peakMx0_wt(5,2)*0.02 0; peakMx0_wt(5,2)*0.02 yMax2; peakMx0_wt(5,1)*0.02 yMax2];
v6_wt = [peakMx0_wt(6,1)*0.02 0; peakMx0_wt(6,2)*0.02 0; peakMx0_wt(6,2)*0.02 yMax2; peakMx0_wt(6,1)*0.02 yMax2];
v7_wt = [peakMx0_wt(7,1)*0.02 0; peakMx0_wt(7,2)*0.02 0; peakMx0_wt(7,2)*0.02 yMax2; peakMx0_wt(7,1)*0.02 yMax2];
f_wt = [1 2 3 4];


v_het = [peakMx0_het(1,1)*0.02 0; peakMx0_het(1,2)*0.02 0; peakMx0_het(1,2)*0.02 yMax2; peakMx0_het(1,1)*0.02 yMax2];
v2_het = [peakMx0_het(2,1)*0.02 0; peakMx0_het(2,2)*0.02 0; peakMx0_het(2,2)*0.02 yMax2; peakMx0_het(2,1)*0.02 yMax2];
v3_het = [peakMx0_het(3,1)*0.02 0; peakMx0_het(3,2)*0.02 0; peakMx0_het(3,2)*0.02 yMax2; peakMx0_het(3,1)*0.02 yMax2];
v4_het = [peakMx0_het(4,1)*0.02 0; peakMx0_het(4,2)*0.02 0; peakMx0_het(4,2)*0.02 yMax2; peakMx0_het(4,1)*0.02 yMax2];
v5_het = [peakMx0_het(5,1)*0.02 0; peakMx0_het(5,2)*0.02 0; peakMx0_het(5,2)*0.02 yMax2; peakMx0_het(5,1)*0.02 yMax2];
v6_het = [peakMx0_het(6,1)*0.02 0; peakMx0_het(6,2)*0.02 0; peakMx0_het(6,2)*0.02 yMax2; peakMx0_het(6,1)*0.02 yMax2];
v7_het = [peakMx0_het(7,1)*0.02 0; peakMx0_het(7,2)*0.02 0; peakMx0_het(7,2)*0.02 yMax2; peakMx0_het(7,1)*0.02 yMax2];
f_het = [1 2 3 4];

figure(292);
patch('Faces', f_wt, 'Vertices', v2_wt*100, 'Facecolor','blue','FaceAlpha',.3);hold on
patch('Faces', f_wt, 'Vertices', v3_wt*100, 'Facecolor','blue','FaceAlpha',.3);
patch('Faces', f_wt, 'Vertices', v4_wt*100, 'Facecolor','blue','FaceAlpha',.3);
patch('Faces', f_wt, 'Vertices', v5_wt*100, 'Facecolor','blue','FaceAlpha',.3);
plot(32:2:66,binAveELCV_wt(16:33),'k','LineWidth',2);hold on
ylim([5 35]);
% legend('wt CV');
xlim([34 64]);
% ylabel('Coefficient of Variation');

figure(293);
patch('Faces', f_het, 'Vertices', v2_het*100, 'Facecolor','blue','FaceAlpha',.3);hold on
patch('Faces', f_het, 'Vertices', v3_het*100, 'Facecolor','blue','FaceAlpha',.3);
patch('Faces', f_het, 'Vertices', v4_het*100, 'Facecolor','blue','FaceAlpha',.3);
patch('Faces', f_het, 'Vertices', v5_het*100, 'Facecolor','blue','FaceAlpha',.3);
plot(32:2:66,binAveELCV_het(16:33),'r','LineWidth',2); 
ylim([5 35]);
% legend('het CV');
xlim([32 62]);
% ylabel('Coefficient of Variation');

%% 2nd half NC14 output - different time frames
for j = 1:5
    for i = 1:size(M_wt01,2)
        OutputT_wt01(j,i) = trapz(M_wt01(51+(j-1)*10:60+(j-1)*10,i));
    end
    for i = 1:size(M_wt02,2)
        OutputT_wt02(j,i) = trapz(M_wt02(51+(j-1)*10:60+(j-1)*10,i));
    end
    for i = 1:size(M_wt03,2)
        OutputT_wt03(j,i) = trapz(M_wt03(51+(j-1)*10:60+(j-1)*10,i));
    end
    
    for i = 1:size(krM_wt01,2)
        krOutputT_wt01(j,i) = trapz(krM_wt01(51+(j-1)*10:60+(j-1)*10,i));
    end
    for i = 1:size(krM_wt02,2)
        krOutputT_wt02(j,i) = trapz(krM_wt02(51+(j-1)*10:60+(j-1)*10,i));
    end
    for i = 1:size(krM_wt03,2)
        krOutputT_wt03(j,i) = trapz(krM_wt03(51+(j-1)*10:60+(j-1)*10,i));
    end
    for i = 1:size(krM_wt04,2)
        krOutputT_wt04(j,i) = trapz(krM_wt04(51+(j-1)*10:60+(j-1)*10,i));
    end
    
    for i = 1:size(krM_het01,2)
        krOutputT_het01(j,i) = trapz(krM_het01(51+(j-1)*10:60+(j-1)*10,i));
    end
    for i = 1:size(krM_het02,2)
        krOutputT_het02(j,i) = trapz(krM_het02(51+(j-1)*10:60+(j-1)*10,i));
    end
    for i = 1:size(krM_het03,2)
        krOutputT_het03(j,i) = trapz(krM_het03(51+(j-1)*10:60+(j-1)*10,i));
    end
    for i = 1:size(krM_het04,2)
        krOutputT_het04(j,i) = trapz(krM_het04(51+(j-1)*10:60+(j-1)*10,i));
    end
    for i = 1:size(krM_het05,2)
        krOutputT_het05(j,i) = trapz(krM_het05(51+(j-1)*10:60+(j-1)*10,i));
    end
end

OutputTS2_wt = [OutputT_wt01(:,actS2second_wt01) OutputT_wt02(:,actS2second_wt02) OutputT_wt03(:,actS2second_wt03) ...
    krOutputT_wt01(:,kractS2second_wt01) krOutputT_wt02(:,kractS2second_wt02) krOutputT_wt03(:,kractS2second_wt03) krOutputT_wt04(:,kractS2second_wt04)];
OutputTS2_het = [krOutputT_het01(:,kractS2second_het01) krOutputT_het02(:,kractS2second_het02) krOutputT_het03(:,kractS2second_het03)...
    krOutputT_het04(:,kractS2second_het04) krOutputT_het05(:,kractS2second_het05)];

OutputTS2P_wt = OutputTS2_wt/10000;%for plot
OutputTS2P_het = OutputTS2_het/10000;%for plot

OutputTS5_wt = [OutputT_wt01(:,actS5second_wt01) OutputT_wt02(:,actS5second_wt02) OutputT_wt03(:,actS5second_wt03) ...
    krOutputT_wt01(:,kractS5second_wt01) krOutputT_wt02(:,kractS5second_wt02) krOutputT_wt03(:,kractS5second_wt03) krOutputT_wt04(:,kractS5second_wt04)];
OutputTS5_het = [krOutputT_het01(:,kractS5second_het01) krOutputT_het02(:,kractS5second_het02) krOutputT_het03(:,kractS5second_het03)...
    krOutputT_het04(:,kractS5second_het04) krOutputT_het05(:,kractS5second_het05)];

OutputTS5P_wt = OutputTS5_wt/10000;%for plot
OutputTS5P_het = OutputTS5_het/10000;%for plot
%% plot scatter of individual nuclei output along the x-axis position - S2
%calibrate S2 nuclei position

%cxELof S2 - calibrated
cxELCS2_wt01 = cxEL_wt01(actS2second_wt01) - min(cxEL_wt01(actS2second_wt01));
cxELCS2_wt02 = cxEL_wt02(actS2second_wt02) - min(cxEL_wt02(actS2second_wt02));
cxELCS2_wt03 = cxEL_wt03(actS2second_wt03) - min(cxEL_wt03(actS2second_wt03));
krcxELCS2_wt01 = krcxEL_wt01(kractS2second_wt01) - min(krcxEL_wt01(kractS2second_wt01));
krcxELCS2_wt02 = krcxEL_wt02(kractS2second_wt02) - min(krcxEL_wt02(kractS2second_wt02));
krcxELCS2_wt03 = krcxEL_wt03(kractS2second_wt03) - min(krcxEL_wt03(kractS2second_wt03));
krcxELCS2_wt04 = krcxEL_wt04(kractS2second_wt04) - min(krcxEL_wt04(kractS2second_wt04));
krcxELCS2_het01 = krcxEL_het01(kractS2second_het01) - min(krcxEL_het01(kractS2second_het01));
krcxELCS2_het02 = krcxEL_het02(kractS2second_het02) - min(krcxEL_het02(kractS2second_het02));
krcxELCS2_het03 = krcxEL_het03(kractS2second_het03) - min(krcxEL_het03(kractS2second_het03));
krcxELCS2_het04 = krcxEL_het04(kractS2second_het04) - min(krcxEL_het04(kractS2second_het04));
krcxELCS2_het05 = krcxEL_het05(kractS2second_het05) - min(krcxEL_het05(kractS2second_het05));

cxELCS2_wt = [cxELCS2_wt01 cxELCS2_wt02 cxELCS2_wt03 ...
    krcxELCS2_wt01 krcxELCS2_wt02 krcxELCS2_wt03 krcxELCS2_wt04];
cxELCS2_het = [krcxELCS2_het01 krcxELCS2_het02 ...
    krcxELCS2_het03 krcxELCS2_het04 krcxELCS2_het05];



%% plot scatter of individual nuclei output along the x-axis position - S5
%calibrate S5 nuclei position

%cxELof S5 - calibrated
cxELCS5_wt01 = cxEL_wt01(actS5second_wt01) - min(cxEL_wt01(actS5second_wt01));
cxELCS5_wt02 = cxEL_wt02(actS5second_wt02) - min(cxEL_wt02(actS5second_wt02));
cxELCS5_wt03 = cxEL_wt03(actS5second_wt03) - min(cxEL_wt03(actS5second_wt03));
krcxELCS5_wt01 = krcxEL_wt01(kractS5second_wt01) - min(krcxEL_wt01(kractS5second_wt01));
krcxELCS5_wt02 = krcxEL_wt02(kractS5second_wt02) - min(krcxEL_wt02(kractS5second_wt02));
krcxELCS5_wt03 = krcxEL_wt03(kractS5second_wt03) - min(krcxEL_wt03(kractS5second_wt03));
krcxELCS5_wt04 = krcxEL_wt04(kractS5second_wt04) - min(krcxEL_wt04(kractS5second_wt04));
krcxELCS5_het01 = krcxEL_het01(kractS5second_het01) - min(krcxEL_het01(kractS5second_het01));
krcxELCS5_het02 = krcxEL_het02(kractS5second_het02) - min(krcxEL_het02(kractS5second_het02));
krcxELCS5_het03 = krcxEL_het03(kractS5second_het03) - min(krcxEL_het03(kractS5second_het03));
krcxELCS5_het04 = krcxEL_het04(kractS5second_het04) - min(krcxEL_het04(kractS5second_het04));
krcxELCS5_het05 = krcxEL_het05(kractS5second_het05) - min(krcxEL_het05(kractS5second_het05));

cxELCS5_wt = [cxELCS5_wt01 cxELCS5_wt02 cxELCS5_wt03 ...
    krcxELCS5_wt01 krcxELCS5_wt02 krcxELCS5_wt03 krcxELCS5_wt04];
cxELCS5_het = [krcxELCS5_het01 krcxELCS5_het02 ...
    krcxELCS5_het03 krcxELCS5_het04 krcxELCS5_het05];


%% bin calibrated S2 and S5 nuclei
binsize = 0.01;

%wt
for i = 1:12
    binELCS2_wt{i} = find(cxELCS2_wt(:) > (i-1)*binsize & cxELCS2_wt(:) <= i*binsize);
    binELCS5_wt{i} = find(cxELCS5_wt(:) > (i-1)*binsize & cxELCS5_wt(:) <= i*binsize);
end
%het
for i = 1:12
    binELCS2_het{i} = find(cxELCS2_het(:) > (i-1)*binsize & cxELCS2_het(:) <= i*binsize);
    binELCS5_het{i} = find(cxELCS5_het(:) > (i-1)*binsize & cxELCS5_het(:) <= i*binsize);
end

%% calculate average output of bin (2nd half NC14)
for i = 1:12
%     binAveEL_wt01(i) = nanmean(trapz(M_wt01(cNaN_wt01+1:end,binEL_wt01{i})));
    binAveELCS2_wt(i) = nanmean(OutputS2P_wt(binELCS2_wt{i}));
    binAveELCS2_het(i) = nanmean(OutputS2P_het(binELCS2_het{i}));
    binAveELCS2SE_wt(i) = std(OutputS2P_wt(binELCS2_wt{i}))/sqrt(length(binELCS2_wt{i}));
    binAveELCS2SE_het(i) = std(OutputS2P_het(binELCS2_het{i}))/sqrt(length(binELCS2_het{i}));

    binAveELCS5_wt(i) = nanmean(OutputS5P_wt(binELCS5_wt{i}));
    binAveELCS5_het(i) = nanmean(OutputS5P_het(binELCS5_het{i}));
    binAveELCS5SE_wt(i) = std(OutputS5P_wt(binELCS5_wt{i}))/sqrt(length(binELCS5_wt{i}));
    binAveELCS5SE_het(i) = std(OutputS5P_het(binELCS5_het{i}))/sqrt(length(binELCS5_het{i}));
end

%% calculate average output of bin (2nd half NC14)
for j = 1:5
    for i = 1:12
        binAveELTCS2_wt(j,i) = nanmean(OutputTS2P_wt(j,binELCS2_wt{i}));
        binAveELTCS2_het(j,i) = nanmean(OutputTS2P_het(j,binELCS2_het{i}));
        binAveELTCS2SE_wt(j,i) = std(OutputTS2P_wt(j,binELCS2_wt{i}))/sqrt(length(binELCS2_wt{i}));
        binAveELTCS2SE_het(j,i) = std(OutputTS2P_het(j,binELCS2_het{i}))/sqrt(length(binELCS2_het{i}));
    
        binAveELTCS5_wt(j,i) = nanmean(OutputTS5P_wt(j,binELCS5_wt{i}));
        binAveELTCS5_het(j,i) = nanmean(OutputTS5P_het(j,binELCS5_het{i}));
        binAveELTCS5SE_wt(j,i) = std(OutputTS5P_wt(j,binELCS5_wt{i}))/sqrt(length(binELCS5_wt{i}));
        binAveELTCS5SE_het(j,i) = std(OutputTS5P_het(j,binELCS5_het{i}))/sqrt(length(binELCS5_het{i}));
    end
end

%plot
figure(521);
plot(0.5:1:9.5,binAveELTCS2_wt(2,1:10),'-k','LineWidth',3);hold on
plot(0.5:1:10.5,binAveELTCS2_het(2,1:11),'-r','LineWidth',3);

shadedErrorBar2(0.5:1:9.5,binAveELTCS2_wt(2,1:10),binAveELTCS2SE_wt(2,1:10),'lineprops', {'-k','LineWidth',3});
shadedErrorBar2(0.5:1:10.5,binAveELTCS2_het(2,1:11),binAveELTCS2SE_het(2,1:11),'lineprops', {'-r','LineWidth',3});
legend('wt','Kr^{1}/+');
ylim([0.1 0.7]);
% xlabel('relative position in {\iteve} stripe 2 (%EL)');
% ylabel('signal intensity (a.u.)');


figure(523);
plot(0.5:1:7.5,binAveELTCS5_wt(2,1:8),'-k','LineWidth',3);hold on
plot(0.5:1:7.5,binAveELTCS5_het(2,1:8),'-r','LineWidth',3);

shadedErrorBar2(0.5:1:7.5,binAveELTCS5_wt(2,1:8),binAveELTCS5SE_wt(2,1:8),'lineprops', {'-k','LineWidth',3});
shadedErrorBar2(0.5:1:7.5,binAveELTCS5_het(2,1:8),binAveELTCS5SE_het(2,1:8),'lineprops', {'-r','LineWidth',3});
legend('wt','Kr^{1}/+');
ylim([0.1 0.7]);


%% trajectory for snapshot

%this section use dataset with un-normalized NC14 to match the real-time
%snapshots of embryos. Here represents the code for analysis and plotting.


%get bin amp for each embryo
for i = 1:size(M_wt01,1)
    for j = 1:bincount
        binAmp_wt01{i,j} = mean(M_wt01(i,binEL_wt01{j}));
    end
end
for i = 1:size(M_wt02,1)
    for j = 1:bincount
        binAmp_wt02{i,j} = mean(M_wt02(i,binEL_wt02{j}));
    end
end
for i = 1:size(M_wt03,1)
    for j = 1:bincount
        binAmp_wt03{i,j} = mean(M_wt03(i,binEL_wt03{j}));
    end
end

for i = 1:size(krM_wt01,1)
    for j = 1:bincount
        krbinAmp_wt01{i,j} = mean(krM_wt01(i,krbinEL_wt01{j}));
    end
end
for i = 1:size(krM_wt02,1)
    for j = 1:bincount
        krbinAmp_wt02{i,j} = mean(krM_wt02(i,krbinEL_wt02{j}));
    end
end
for i = 1:size(krM_wt03,1)
    for j = 1:bincount
        krbinAmp_wt03{i,j} = mean(krM_wt03(i,krbinEL_wt03{j}));
    end
end
for i = 1:size(krM_wt04,1)
    for j = 1:bincount
        krbinAmp_wt04{i,j} = mean(krM_wt04(i,krbinEL_wt04{j}));
    end
end

for i = 1:size(krM_het01,1)
    for j = 1:bincount
        krbinAmp_het01{i,j} = mean(krM_het01(i,krbinEL_het01{j}));
    end
end
for i = 1:size(krM_het02,1)
    for j = 1:bincount
        krbinAmp_het02{i,j} = mean(krM_het02(i,krbinEL_het02{j}));
    end
end
for i = 1:size(krM_het03,1)
    for j = 1:bincount
        krbinAmp_het03{i,j} = mean(krM_het03(i,krbinEL_het03{j}));
    end
end
for i = 1:size(krM_het04,1)
    for j = 1:bincount
        krbinAmp_het04{i,j} = mean(krM_het04(i,krbinEL_het04{j}));
    end
end
for i = 1:size(krM_het05,1)
    for j = 1:bincount
        krbinAmp_het05{i,j} = mean(krM_het05(i,krbinEL_het05{j}));
    end
end

%plot
snap1 = 24;
snap2 = 31;
snap3 = 38;

figure(231);
plot(2:2:100,[krbinAmp_wt04{snap1-3,:}]/1000,'k','LineWidth',4); hold on
plot(2:2:100,[krbinAmp_het01{snap1-4,:}]/1000,'r','LineWidth',4);
legend('/itwt','/itKr^{1}/+');
ylim([0.09 1]);
ax = gca;
ax.YAxis.FontSize = 20;
ax.XAxis.FontSize = 20;
ax.Legend.FontSize = 20;
ax.Title.FontSize = 20;
hold off

figure(232);
plot(2:2:100,[krbinAmp_wt04{snap2-3,:}]/1000,'k','LineWidth',4); hold on
plot(2:2:100,[krbinAmp_het01{snap2-4,:}]/1000,'r','LineWidth',4);
ylim([0.09 1]);
ax = gca;
ax.YAxis.FontSize = 20;
ax.XAxis.FontSize = 20;
hold off

figure(233);
plot(2:2:100,[krbinAmp_wt04{snap3-3,:}]/1000,'k','LineWidth',4); hold on
plot(2:2:100,[krbinAmp_het01{snap3-4,:}]/1000,'r','LineWidth',4);
ylim([0.09 1]);
ax = gca;
ax.YAxis.FontSize = 20;
ax.XAxis.FontSize = 20;
hold off