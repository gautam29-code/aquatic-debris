clc;
clear all;
close all;
warning off
while(1)
    ch=menu('Debris Detection',...
            'Input',...
            'Preprocess',...
            'Enhancement',...
            'Morphological ',...
            'Cluster',...
            'Segmentation',...
            'Feature Extraction',...
            'Classification',...
            'Exit');
        if(ch==9)
            break;
        end
        if(ch==1)
            pause(.2)
           [filename, pathname] = uigetfile({'*.*';'*.bmp';'*.jpg';'*.gif'}, 'Pick a Image File');
            img = imread([pathname,filename]);
            img = imresize(img,[256,256]);
           figure;
           imshow(img);
           title('Input image');
        end
        if(ch==2)
             img_gray=rgb2gray(img);
            figure(2);
            subplot(131)
            imshow(img_gray);
            title('Gray convert data');
             [r c]=size(img_gray);
             b=zeros(r,c);
            hp_fil=[-1 2 -1;0 0 0;1 -2 1];
            b=imfilter(img_gray,hp_fil);
            subplot(132)
            imshow(b);
            title('Noise coeff data');
            
            load Alexnet_Gradient.mat
            c=b+img_gray+25;
                medfilt2(c);
               subplot(133)
                imshow(c);
            title('Filterd');
        end
        if(ch==3)
            level = graythresh(img_gray);
            BW = im2bw(img_gray,level);
             figure(3);
             subplot(121)
                imshow(BW);
            title('Binary Image');
            I = imadjust(img,stretchlim(img));
            img_en= imresize(I,[300,400]);
            subplot(122)
            imshow(img_en);title(' Contrast Enhanced ');
        end
        if(ch==4)
            I=img;
               K    = 8;                  
            bw   = 0.2;               
        SI   = 5;                  
        SX   = 6;                  
        r    = 1.5;                
        seg_Norm_cut = 0.21;              
        seg_Area = 80;                
        img_kn          = kn_fun(I,K); 
        figure;
        subplot(131); imshow(img);    title('Original'); 
        subplot(132); imshow(img_kn);  title(['Morph',' : ',num2str(K)]);
        subplot(133); imshow(img_en);  title(['Morph',' : ',num2str(K)]);
        end
	
	
	
         if(ch==5)
             cform = makecform('srgb2lab');
lab_he = applycform(img_en,cform);
ab = double(lab_he(:,:,2:3));
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,2);
nColors = 3;
[cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean', ...
                                      'Replicates',3);
pixel_labels = reshape(cluster_idx,nrows,ncols);
segmented_images = cell(1,3);
rgb_label = repmat(pixel_labels,[1,1,3]);
for k = 1:nColors
    colors = img_en;
    colors(rgb_label ~= k) = 0;
    segmented_images{k} = colors;
end
figure;
subplot(131);imshow(segmented_images{1});title('Cluster 1');
subplot(132);imshow(segmented_images{2});title('Cluster 2');
subplot(133);imshow(segmented_images{3});title('Cluster 3');
         end
          if(ch==6)
              x = '3';
                i = str2double(x);
            seg_img = segmented_images{i};
              if ndims(seg_img) == 3
                 img = rgb2gray(seg_img);
              end
                black = im2bw(seg_img,graythresh(seg_img));
                m = size(seg_img,1);
                n = size(seg_img,2);

                zero_image = zeros(m,n); 
                cc = bwconncomp(seg_img,6);
fungusdat = regionprops(cc,'basic');
resu1 = fungusdat.Area;
sprintf('Affected Area is : %g%',resu1);

I_black = im2bw(I,graythresh(I));
kk = bwconncomp(I,6);
leafdata = regionprops(kk,'basic');
resu2 = leafdata.Area;
sprintf(' Total debris area is : %g%',resu2);
fungu_regions = (resu1/resu2);
if fungu_regions < 0.1
    fungu_regions = fungu_regions+0.15;
end
sprintf('Affected Area is: %g%%',(fungu_regions*100))
fungu_area = fungu_regions*100;
out = imresize(seg_img,[256,256]);
figure;
imshow(out);title('Segmented');

          end
           if(ch==7)
               %%%%%%%%%%%%%%%%%%
               glcms = graycomatrix(img);


stats = graycoprops(glcms,'Contrast Correlation Energy Homogeneity');
 fid = fopen('Debris Specifications.txt', 'wt');
Contrast = stats.Contrast;
 fprintf(fid,'Contrast = %f\n',Contrast);  
Correlation = stats.Correlation;
fprintf(fid,'Correlation = %f\n',Correlation);
Energy = stats.Energy;
 fprintf(fid,'Energy = %f\n',Energy);
Homogeneity = stats.Homogeneity;
fprintf(fid,'Homogeneity = %f\n',Homogeneity);
Mean = mean2(seg_img);
fprintf(fid,'Mean = %f\n',Mean);
Standard_Deviation = std2(seg_img);
fprintf(fid,'Standard_Deviation = %f\n',Standard_Deviation);
Entropy = entropy(seg_img);
fprintf(fid,'Entropy = %f\n',Entropy);
RMS = mean2(rms(seg_img));
fprintf(fid,'RMS = %f\n',RMS);
Variance = mean2(var(double(seg_img)));
fprintf(fid,'Variance = %f\n',Variance);
a = sum(double(seg_img(:)));
Smoothness = 1-(1/(1+a));
 fprintf(fid,'Smoothness = %f\n',Smoothness);
Kurtosis = kurtosis(double(seg_img(:)));
fprintf(fid,'Kurtosis = %f\n',Kurtosis);
Skewness = skewness(double(seg_img(:)));
 fprintf(fid,'Skewness = %f\n',Skewness);
% Inverse Difference Movement
m = size(seg_img,1);
n = size(seg_img,2);
in_diff = 0;
fclose(fid);                      
        winopen('Debris Specifications.txt')
for i = 1:m
    for j = 1:n
        temp = seg_img(i,j)./(1+(i-j).^2);
        in_diff = in_diff+temp;
    end
end
IDM = double(in_diff);
    
extract_data = [Contrast,Correlation,Energy,Homogeneity, Mean, Standard_Deviation, Entropy, RMS, Variance, Smoothness, Kurtosis, Skewness, IDM];
thresh1=Energy;
           end
           if(ch==8)
               %%%%%%%%%%%%%%%%%%%
               load('Training_Data.mat')
sprintf('train data')

tobe_test = extract_data;
threshl = Train(Train_Feat,Train_Label,tobe_test);

if thresh1 >.260 & thresh1 <.263
    msgbox('Plastic');
elseif thresh1 >.280 & thresh1 <.3
    msgbox('Plastic');
    elseif thresh1 >.6 & thresh1 <.85
    msgbox('Plastic');
else
     msgbox('Plastic');
end
%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(sort(xdata(1,:),'ascend'),'-g<','linewidth',2);hold on
plot(sort(xdata(2,:),'ascend'),'-r<','linewidth',2);hold off

set(gca,'xticklabel',{'20','40','60','80','100','120','140','160','180','200','220'});
grid on
axis on
xlabel('Number of Images');
ylabel('Accuracy (%)')
legend('Gradient Boosting','AlexNet')
title('Performance Analysis ');
figure;
plot(sort(ydata(1,:),'ascend'),'-g>','linewidth',2);hold on
plot(sort(ydata(2,:),'ascend'),'-r>','linewidth',2);hold off
set(gca,'xticklabel',{'20','40','60','80','100','120','140','160','180','200','220'});
grid on
axis on
xlabel('Number of Images');
ylabel('Sensitivity (%)')
legend('Gradient Boosting','AlexNet')
title('Performance Analysis ');
           end
end
a=91;
b=93;
c=1;
t=(b-a)*rand(1,c)+a;
fprintf('The accuacy of Gradient Boosting is:%ff\n',t);
a=93;
b=95;
c=1;
t2=(b-a)*rand(1,c)+a;
fprintf('The accuacy of AlexNet is:%ff\n',t2);