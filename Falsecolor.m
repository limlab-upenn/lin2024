clear all
close all

path = 'E:\LIM LAB\new-kr-whole embryo\het05\';

load([path 'trajectories.mat']);
load([path 'segmentation_lineage.mat']);
load('all_stripes_nuclei_index_NC14_100f.mat');

M = double(M); Mm = double(Mm);
for i=1:size(M,2)

    M1(:,i) = M(:,i)-Mm(:,i);
    M1(:,i) = M1(:,i)-min(M1(:,i));
    
    sM(:,i) = smooth(M1(:,i));
end

for i=round(0.6*size(M,1))%pick frames of interest, or loop to get all the frames
    F = zeros(size(images_seg,1),size(images_seg,2),3,'uint16');
    
    for j=1:size(M,2) % all nuclei
        c_image = convex_image{i,nuc_lineage(i,(j))};
        c_image = imresize(c_image,1);
        r1 = regionprops(c_image,'PixelList','centroid');
        
        dx = lineage_cx(i,(j))-r1.Centroid(1); dx = round(dx);
        dy = lineage_cy(i,(j))-r1.Centroid(2); dy = round(dy);
        r1.PixelList(:,1) = r1.PixelList(:,1)+dx;
        r1.PixelList(:,2) = r1.PixelList(:,2)+dy;

        ddum = find(r1.PixelList(:,1)<1);
        r1.PixelList(ddum,1) = 1;
        ddum= find(r1.PixelList(:,2)<1);
        r1.PixelList(ddum,2)=1;
        ddum = find(r1.PixelList(:,1)>size(images_ms2,2));
        r1.PixelList(ddum,1) = size(images_ms2,2);
        ddum = find(r1.PixelList(:,2)>size(images_seg,1));
        r1.PixelList(ddum,2) = size(images_ms2,1);
        
        xx = r1.PixelList(:,1);
        yy = r1.PixelList(:,2);
        
   
            tmp = [M1(i,j) M1(i-1,j)];
            for k=1:length(xx)
                %set the color of active nuclei in specific stripe
                if ismember(j,kract_wt03)==1 && ismember(j,krstripe1_wt03)==1
%                     F(yy(k),xx(k),1) = 35000;%red
%                     F(yy(k),xx(k),2) = 35000;%green
%                     F(yy(k),xx(k),1) = 35000;%grey
%                     F(yy(k),xx(k),2) = 35000;%grey
%                     F(yy(k),xx(k),3) = 35000;%grey
                    F(yy(k),xx(k),1) = 8000;%dark
                    F(yy(k),xx(k),2) = 8000;%dark
                    F(yy(k),xx(k),3) = 8000;%dark
                elseif ismember(j,kract_wt03)==1 && ismember(j,krstripe2_wt03)==1
%                     F(yy(k),xx(k),1) = 35000;
%                     F(yy(k),xx(k),2) = 35000;
%                     F(yy(k),xx(k),2) = 35000;%green
%                     F(yy(k),xx(k),1) = 35000;%red
%                     F(yy(k),xx(k),1) = 35000;%grey
%                     F(yy(k),xx(k),2) = 35000;%grey
%                     F(yy(k),xx(k),3) = 35000;%grey
%                     F(yy(k),xx(k),2) = 35000;%cyan
%                     F(yy(k),xx(k),3) = 35000;%cyan
%                     F(yy(k),xx(k),1) = 35000;%orange
%                     F(yy(k),xx(k),2) = 17300;%orange
                    F(yy(k),xx(k),1) = 8000;%dark
                    F(yy(k),xx(k),2) = 8000;%dark
                    F(yy(k),xx(k),3) = 8000;%dark
                elseif ismember(j,kract_wt03)==1 && ismember(j,krstripe3_wt03)==1
%                     F(yy(k),xx(k),2) = 35000;
%                     F(yy(k),xx(k),1) = 35000;%yellow
%                     F(yy(k),xx(k),2) = 35000;%yellow
%                     F(yy(k),xx(k),2) = 35000;%cyan
%                     F(yy(k),xx(k),3) = 35000;%cyan
%                     F(yy(k),xx(k),1) = 35000;%grey
%                     F(yy(k),xx(k),2) = 35000;%grey
%                     F(yy(k),xx(k),3) = 35000;%grey
%                     F(yy(k),xx(k),1) = 35000;%orange
%                     F(yy(k),xx(k),2) = 17300;%orange
                    F(yy(k),xx(k),1) = 8000;%dark
                    F(yy(k),xx(k),2) = 8000;%dark
                    F(yy(k),xx(k),3) = 8000;%dark
                elseif ismember(j,kract_wt03)==1 && ismember(j,krstripe4_wt03)==1
%                     F(yy(k),xx(k),2) = 35000;%green
                    F(yy(k),xx(k),1) = 35000;%yellow
                    F(yy(k),xx(k),2) = 35000;%yellow
%                     F(yy(k),xx(k),1) = 35000;%grey
%                     F(yy(k),xx(k),2) = 35000;%grey
%                     F(yy(k),xx(k),3) = 35000;%grey
%                     F(yy(k),xx(k),1) = 35000;%orange
%                     F(yy(k),xx(k),2) = 17300;%orange
%                     F(yy(k),xx(k),2) = 35000;
%                     F(yy(k),xx(k),3) = 35000;
%                     F(yy(k),xx(k),1) = 8000;%dark
%                     F(yy(k),xx(k),2) = 8000;%dark
%                     F(yy(k),xx(k),3) = 8000;%dark
                elseif ismember(j,kract_wt03)==1 && ismember(j,krstripe5_wt03)==1
%                     F(yy(k),xx(k),2) = 35000;%cyan
%                     F(yy(k),xx(k),3) = 35000;%cyan
%                     F(yy(k),xx(k),2) = 35000;%green
%                     F(yy(k),xx(k),3) = 35000;%purple
%                     F(yy(k),xx(k),1) = 35000;%purple
%                     F(yy(k),xx(k),1) = 35000;%grey
%                     F(yy(k),xx(k),2) = 35000;%grey
%                     F(yy(k),xx(k),3) = 35000;%grey
%                     F(yy(k),xx(k),2) = 8000;
%                     F(yy(k),xx(k),3) = 50000;
                    F(yy(k),xx(k),1) = 8000;%dark
                    F(yy(k),xx(k),2) = 8000;%dark
                    F(yy(k),xx(k),3) = 8000;%dark
                elseif ismember(j,kract_wt03)==1 && ismember(j,krstripe6_wt03)==1
%                     F(yy(k),xx(k),3) = 35000;%purple
%                     F(yy(k),xx(k),1) = 35000;%purple
%                     F(yy(k),xx(k),1) = 35000;%orange
%                     F(yy(k),xx(k),2) = 17300;%orange
%                     F(yy(k),xx(k),1) = 35000;%yellow
%                     F(yy(k),xx(k),2) = 35000;%yellow
%                     F(yy(k),xx(k),1) = 35000;%grey
%                     F(yy(k),xx(k),2) = 35000;%grey
%                     F(yy(k),xx(k),3) = 35000;%grey
                    F(yy(k),xx(k),1) = 8000;%dark
                    F(yy(k),xx(k),2) = 8000;%dark
                    F(yy(k),xx(k),3) = 8000;%dark
                elseif ismember(j,kract_wt03)==1 && j > 0.5*size(M,2) && ismember(j,[krstripe1_wt03,krstripe2_wt03,krstripe3_wt03...
                        krstripe4_wt03,krstripe5_wt03,krstripe6_wt03])==0
%                     F(yy(k),xx(k),1) = 35000;%grey
%                     F(yy(k),xx(k),2) = 35000;%grey
%                     F(yy(k),xx(k),3) = 35000;%grey
                    F(yy(k),xx(k),1) = 8000;%dark
                    F(yy(k),xx(k),2) = 8000;%dark
                    F(yy(k),xx(k),3) = 8000;%dark
                else %nuclei not belong to any stripe
                    F(yy(k),xx(k),1) = 8000;
                    F(yy(k),xx(k),2) = 8000;
                    F(yy(k),xx(k),3) = 8000;
                end
            end
        end        

        imwrite(F,sprintf('%s%03d.png','E:\LIM LAB\kr_project_code_cleanup\test\',i));
end

                
            