clear
close all

%% load files and set parameters
% trajectories used in this script is smoothed
load('dataset_whole_kr_gtMS2_real_length_forSnapshotRotationOnly.mat');
videopath = 'D:\LIM LAB\kr_gtMS2\kr_gtMS2_code\kr_gtMS2_WTvsHET_Amp_NC14_100f.avi'; 

savefile = 0;

thr = 1500;
binsize = 0.02;
bincount = 1/binsize;

checkSmooth = 0;
showSEG = 0;
checkExtreme = 0;
checkFlip = 0;
checkTrajectory = 0;
runVideo = 0;
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


% %% extract max intensity
% for i=1:size(M_wt01,2) % all nuclei
%   maxM_wt01(i) = max(M_wt01(:,i)); % obtain max intensity for each nuclues
%   maxM2_wt01(i) = max(M_wt01(round(0.75*size(M_wt01,1)):end,i));%max intensity for last 1/4 frames
% end 
% 
% for i=1:size(M_wt02,2) % all nuclei
%   maxM_wt02(i) = max(M_wt02(:,i)); % obtain max intensity for each nuclues
%   maxM2_wt02(i) = max(M_wt02(round(0.75*size(M_wt02,1)):end,i));%max intensity for last 1/4 frames
% end 
% 
% for i=1:size(M_wt03,2) % all nuclei
%   maxM_wt03(i) = max(M_wt03(:,i)); % obtain max intensity for each nuclues
%   maxM2_wt03(i) = max(M_wt03(round(0.75*size(M_wt03,1)):end,i));%max intensity for last 1/4 frames
% end 
% 
% for i=1:size(M_het01,2) % all nuclei
%   maxM_het01(i) = max(M_het01(:,i)); % obtain max intensity for each nuclues
%   maxM2_het01(i) = max(M_het01(round(0.75*size(M_het01,1)):end,i));%max intensity for last 1/4 frames
% end 
% 
% for i=1:size(M_het02,2) % all nuclei
%   maxM_het02(i) = max(M_het02(:,i)); % obtain max intensity for each nuclues
%   maxM2_het02(i) = max(M_het02(round(0.75*size(M_het02,1)):end,i));%max intensity for last 1/4 frames
% end 
% 
% for i=1:size(M_het03,2) % all nuclei
%   maxM_het03(i) = max(M_het03(:,i)); % obtain max intensity for each nuclues
%   maxM2_het03(i) = max(M_het03(round(0.75*size(M_het03,1)):end,i));%max intensity for last 1/4 frames
% end 



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
lineageCx_wt02 = abs(lineageCx_wt02 - 950);
lineageCx_het03 = abs(lineageCx_het03 - 950);

%flip nuclei fixed cx
cx_wt01 = abs(cx_wt01 - 950); %flip x axis
cx_wt02 = abs(cx_wt02 - 950); %flip x axis
cx_het03 = abs(cx_het03 - 950); %flip x axis

%flip tip points
tip_wt01 = abs(tip_wt01 - 1400); %flip x axis
tip_wt02 = abs(tip_wt02 - 1400); %flip x axis
tip_het03 = abs(tip_het03 - 1400); %flip x axis

%flip nuc image
BW_wt01 = flip(BW_wt01,2);
BW_wt02 = flip(BW_wt02,2);
BW_het03 = flip(BW_het03,2);

%flip seg image
seg_wt01 = flip(seg_wt01,2);
seg_wt02 = flip(seg_wt02,2);
seg_het03 = flip(seg_het03,2);

if checkFlip == 1
    figure(31),subplot(2,2,1);imshow(BW_wt01);hold on
    xline(tip_wt01(1),'r','LineWidth',2); xline(tip_wt01(2),'r','LineWidth',2);
    plot(cx_wt01+225,cy_wt01,'ro');
    title('wt01');hold off        
    subplot(2,2,2);imshow(BW_wt02);hold on
    xline(tip_wt02(1),'r','LineWidth',2); xline(tip_wt02(2),'r','LineWidth',2);
    plot(cx_wt02+225,cy_wt02,'ro');
    title('wt02');hold off   
    subplot(2,2,3);imshow(BW_wt03);hold on
    xline(tip_wt03(1),'r','LineWidth',2); xline(tip_wt03(2),'r','LineWidth',2);
    plot(cx_wt03+225,cy_wt03,'ro');
    title('wt03');hold off

    figure(32),subplot(2,2,1);imshow(BW_het01);hold on
    xline(tip_het01(1),'r','LineWidth',2); xline(tip_het01(2),'r','LineWidth',2);
    plot(cx_het01+225,cy_het01,'ro');
    title('het01');hold off        
    subplot(2,2,2);imshow(BW_het02);hold on
    xline(tip_het02(1),'r','LineWidth',2); xline(tip_het02(2),'r','LineWidth',2);
    plot(cx_het02+225,cy_het02,'ro');
    title('het02');hold off   
    subplot(2,2,3);imshow(BW_het03);hold on
    xline(tip_het03(1),'r','LineWidth',2); xline(tip_het03(2),'r','LineWidth',2);
    plot(cx_het03+225,cy_het03,'ro');
    title('het03');hold off
end

%% calibrate x position by EL
cxEL_wt01 = (cx_wt01+225 - min(tip_wt01))/(max(tip_wt01)- min(tip_wt01));
cxEL_wt02 = (cx_wt02+225 - min(tip_wt02))/(max(tip_wt02)- min(tip_wt02));
cxEL_wt03 = (cx_wt03+225 - min(tip_wt03))/(max(tip_wt03)- min(tip_wt03));

cxEL_het01 = (cx_het01+225 - min(tip_het01))/(max(tip_het01)- min(tip_het01));
cxEL_het02 = (cx_het02+225 - min(tip_het02))/(max(tip_het02)- min(tip_het02));
cxEL_het03 = (cx_het03+225 - min(tip_het03))/(max(tip_het03)- min(tip_het03));

%% bin by 2% EL

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

%het01
for i = 1:50
    binEL_het01{i} = find(cxEL_het01(:) > (i-1)*binsize & cxEL_het01(:) <= i*binsize);
end
%het02
for i = 1:50
    binEL_het02{i} = find(cxEL_het02(:) > (i-1)*binsize & cxEL_het02(:) <= i*binsize);
end
%het03
for i = 1:50
    binEL_het03{i} = find(cxEL_het03(:) > (i-1)*binsize & cxEL_het03(:) <= i*binsize);
end

for i = 1:42
    for j = 1:bincount
        binAmp_wt01{i,j} = mean(M_wt01(i,binEL_wt01{j}));
        binAmp_wt02{i,j} = mean(M_wt02(i,binEL_wt02{j}));
        binAmp_wt03{i,j} = mean(M_wt03(i,binEL_wt03{j}));

        binAmp_het01{i,j} = mean(M_het01(i,binEL_het01{j}));
        binAmp_het02{i,j} = mean(M_het02(i,binEL_het02{j}));
        binAmp_het03{i,j} = mean(M_het03(i,binEL_het03{j}));

    end
end

%% get boundary, 50% max, at f19

fixB = 1;
i = 19;
b = 500;
blist = b*ones(1,50);

if fixB == 0
    BoundaryValue_wt01(1:50) = 0.5*max([binAmp_wt01{i,:}])*ones(1,50);% kni boundary signal array
    BoundaryValue_wt02(1:50) = 0.5*max([binAmp_wt02{i,:}])*ones(1,50);% kni boundary signal array
    BoundaryValue_wt03(1:50) = 0.5*max([binAmp_wt03{i,:}])*ones(1,50);% kni boundary signal array
    BoundaryValue_het01(1:50) = 0.5*max([binAmp_het01{i,:}])*ones(1,50);% kni boundary signal array
    BoundaryValue_het02(1:50) = 0.5*max([binAmp_het02{i,:}])*ones(1,50);% kni boundary signal array
    BoundaryValue_het03(1:50) = 0.5*max([binAmp_het03{i,:}])*ones(1,50);% kni boundary signal array
elseif fixB == 1;
    BoundaryValue_wt01(1:50) = blist;
    BoundaryValue_wt02(1:50) = blist;
    BoundaryValue_wt03(1:50) = blist;
    BoundaryValue_het01(1:50) = blist;
    BoundaryValue_het02(1:50) = blist;
    BoundaryValue_het03(1:50) = blist;
end
[Bx_wt01,By_wt01] = intersections(1:50,[binAmp_wt01{i,1:50}],1:50,BoundaryValue_wt01(1:50));%find x y of intersection
[Bx_wt02,By_wt02] = intersections(1:50,[binAmp_wt02{i,1:50}],1:50,BoundaryValue_wt02(1:50));%find x y of intersection
[Bx_wt03,By_wt03] = intersections(1:50,[binAmp_wt03{i,1:50}],1:50,BoundaryValue_wt03(1:50));%find x y of intersection

[Bx_het01,By_het01] = intersections(1:50,[binAmp_het01{i,1:50}],1:50,BoundaryValue_het01(1:50));%find x y of intersection
[Bx_het02,By_het02] = intersections(1:50,[binAmp_het02{i,1:50}],1:50,BoundaryValue_het02(1:50));%find x y of intersection
[Bx_het03,By_het03] = intersections(1:50,[binAmp_het03{i,1:50}],1:50,BoundaryValue_het03(1:50));%find x y of intersection

Bx_wt01 = 0.02*Bx_wt01;
Bx_wt02 = 0.02*Bx_wt02;
Bx_wt03 = 0.02*Bx_wt03;
Bx_het01 = 0.02*Bx_het01;
Bx_het02 = 0.02*Bx_het02;
Bx_het03 = 0.02*Bx_het03;

%% calibrate Bx to actual image pixel
Bx_wt01 = Bx_wt01*(max(tip_wt01)- min(tip_wt01)) + min(tip_wt01)-225;
Bx_wt02 = Bx_wt02*(max(tip_wt02)- min(tip_wt02)) + min(tip_wt02)-225;
Bx_wt03 = Bx_wt03*(max(tip_wt03)- min(tip_wt03)) + min(tip_wt03)-225;

Bx_het01 = Bx_het01*(max(tip_het01)- min(tip_het01)) + min(tip_het01)-225;
Bx_het02 = Bx_het02*(max(tip_het02)- min(tip_het02)) + min(tip_het02)-225;
Bx_het03 = Bx_het03*(max(tip_het03)- min(tip_het03)) + min(tip_het03)-225;

%% wt01 f19 (original wt02)
for i = 0:10
    Scale_wt01(i+1) = i*0.1*(max(tip_wt01)- min(tip_wt01)) + min(tip_wt01)-225;
end

ypoint = 225.5;
xpoint = 475.5;
deg = 6.5;

I = imread("D:\LIM LAB\figures_arrangement\gap_snapshots\gt\mean value of bin\NewHis\wt01_gt_NC14_newHis_f19.png");%wt03 of raw = wt01
% I = imread("C:\Users\Slin\OneDrive - PennO365\Desktop\figures_arrangement\gap_snapshots\gt\falsecolor_wt01_gt_notfullNC14_f16_equals_f19.png");%wt03 of raw = wt01
I1 = rotateAround(I(:,:,1), ypoint, xpoint, -1*deg, 'nearest');
I2 = rotateAround(I(:,:,2), ypoint, xpoint, -1*deg, 'nearest');
I3 = rotateAround(I(:,:,3), ypoint, xpoint, -1*deg, 'nearest');
II(:,:,1) = I1;
II(:,:,2) = I2;
II(:,:,3) = I3;
IIf = flipdim(II,2);
IIf2 = flipdim(IIf,1);
% IIf = II;

figure(12),imshow(IIf2);hold on
plot([Scale_wt01(2) Scale_wt01(10)],[ypoint ypoint],'y', 'lineWidth',6);
plot([Bx_wt01(3) Bx_wt01(3)],[min(cy_wt01) max(cy_wt01)+30],'w', 'lineWidth',6);
plot([Bx_wt01(4) Bx_wt01(4)],[min(cy_wt01)+20 max(cy_wt01)-15],'w', 'lineWidth',6);

for i = 2:10
    if ismember(i,[2 6 10])
        plot([Scale_wt01(i) Scale_wt01(i)],[ypoint ypoint-20],'y', 'lineWidth',6);
        tmp = num2str((i-1)*10);
        xt = round(Scale_wt01(i))-25;
        yt = round(ypoint+25);
        text(xt,yt,tmp,'Color','y','FontSize', 60);
    else
        plot([Scale_wt01(i) Scale_wt01(i)],[ypoint ypoint-10],'y', 'lineWidth',4);
    end
end
%% wt02 f33 (original wt03)
for i = 0:10
    Scale_wt02(i+1) = i*0.1*(max(tip_wt02)- min(tip_wt02)) + min(tip_wt02)-225;
end

ypoint = 225.5;
xpoint = 475.5;
deg = 7;

I = imread("D:\LIM LAB\figures_arrangement\gap_snapshots\gt\mean value of bin\NewHis\wt03_gt_NC14_newHis_f33.png");%wt03 of raw = wt02
I1 = rotateAround(I(:,:,1), ypoint, xpoint, -1*deg, 'nearest');
I2 = rotateAround(I(:,:,2), ypoint, xpoint, -1*deg, 'nearest');
I3 = rotateAround(I(:,:,3), ypoint, xpoint, -1*deg, 'nearest');
II(:,:,1) = I1;
II(:,:,2) = I2;
II(:,:,3) = I3;
IIf = flipdim(II,2);
% IIf = II;

figure(22),imshow(IIf);hold on
plot([Scale_wt02(2) Scale_wt02(10)],[ypoint ypoint],'y', 'lineWidth',6);
% plot([Bx_wt02(3) Bx_wt02(3)],[ypoint+30 ypoint-30],'w', 'lineWidth',6);
% plot([Bx_wt02(4) Bx_wt02(4)],[ypoint+30 ypoint-30],'w', 'lineWidth',6);

for i = 2:10
    if ismember(i,[2 6 10])
        plot([Scale_wt02(i) Scale_wt02(i)],[ypoint ypoint-20],'y', 'lineWidth',6);
        tmp = num2str((i-1)*10);
        xt = round(Scale_wt02(i))-25;
        yt = round(ypoint+25);
        text(xt,yt,tmp,'Color','y','FontSize', 60);
    else
        plot([Scale_wt02(i) Scale_wt02(i)],[ypoint ypoint-10],'y', 'lineWidth',4);
    end
end


%% het03 f19 (original het05)
for i = 0:10
    Scale_het03(i+1) = i*0.1*(max(tip_het03)- min(tip_het03)) + min(tip_het03)-225;
end

ypoint = 225.5;
xpoint = 475.5;
deg = 0;

I = imread("D:\LIM LAB\figures_arrangement\gap_snapshots\gt\mean value of bin\NewHis\het05_gt_NC14_newHis_f19.png");%het03 of raw = het05
% I = imread("C:\Users\Slin\OneDrive - PennO365\Desktop\figures_arrangement\gap_snapshots\gt\falsecolor_het05_gt_notfullNC14_f16_equals_f19.png");
I1 = rotateAround(I(:,:,1), ypoint, xpoint, -1*deg, 'nearest');
I2 = rotateAround(I(:,:,2), ypoint, xpoint, -1*deg, 'nearest');
I3 = rotateAround(I(:,:,3), ypoint, xpoint, -1*deg, 'nearest');
II(:,:,1) = I1;
II(:,:,2) = I2;
II(:,:,3) = I3;
IIf = flipdim(II,2);
% IIf = II;

figure(62),imshow(IIf);hold on
plot([Scale_het03(2) Scale_het03(10)],[ypoint ypoint],'y', 'lineWidth',6);
plot([Bx_het03(3) Bx_het03(3)],[min(cy_het03)-10 max(cy_het03)+5],'w', 'lineWidth',6);
plot([Bx_het03(4) Bx_het03(4)],[min(cy_het03) max(cy_het03)],'w', 'lineWidth',6);

for i = 2:10
    if ismember(i,[2 6 10])
        plot([Scale_het03(i) Scale_het03(i)],[ypoint ypoint-20],'y', 'lineWidth',6);
        tmp = num2str((i-1)*10);
        xt = round(Scale_het03(i))-25;
        yt = round(ypoint+25);
        text(xt,yt,tmp,'Color','y','FontSize', 60);
    else
        plot([Scale_het03(i) Scale_het03(i)],[ypoint ypoint-10],'y', 'lineWidth',4);
    end
end
