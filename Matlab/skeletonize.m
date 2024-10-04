clear;
I = imread('D:\Data\Processing Testing\Test Images\png\mosaic2-1.tif');

Icomplement = imcomplement(I);
BW = imbinarize(Icomplement);
out = bwskel(BW);
mn=bwmorph(out,'branchpoints');
branches = out & ~mn; % set branch points to zero
imshow(branches);

sts = regionprops( branches, "Orientation", "Area");

figure;
angles = [];

for i = 1:length(sts) 
    if sts(i).Orientation < 0
        angles = [angles repmat(sts(i).Orientation+180, 1, sts(i).Area)];
    else
        angles = [angles repmat(sts(i).Orientation, 1, sts(i).Area)];
    end
end
histogram(angles)
writematrix(angles, "Skeletonize.csv")