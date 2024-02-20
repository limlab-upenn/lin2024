clear
close all

path = 'D:\LIM LAB\new-kr-whole embryo\het06\';

load([path 'trajectories.mat']);
load([path 'segmentation_lineage.mat']);
load('myCustomColormap.mat');

M = double(M); Mm = double(Mm);
for i=1:size(M,2)

    M1(:,i) = M(:,i)-Mm(:,i);
    M2(:,i) = M1(:,i)-min(M1(:,i));

    M3(:,i) = interp1(1:size(M2,1), M2(:,i), linspace(1,size(M2,1),100));   
    sM(:,i) = smooth(M3(:,i));

    
end

%calculate mRNA for each nucleus
for i = 1:size(M1,2)
    Output(i) = trapz(sM(51:100,i)); % 2nd half of NC14, normalized
end

[sortOutput sortIndexes] = sort(Output); 

cmap = parula(size(Output,2));

%create blank figure
F = zeros(size(images_seg,1),size(images_seg,2),3,'uint16');

i = round(0.6*size(M,1));
    
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
        
        %find output index to get the color
        ind = find(sortIndexes == j);
        
        %color nuclei
        for k=1:length(xx)
                F(yy(k),xx(k),1) = Output(j);
        end        
    end

%check un-rotated heatmap
figure,imagesc(F(:,:,1));colormap(gca, myCustomColormap);
caxis([4000 50000]);
axis off
colorbar
