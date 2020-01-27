%GETCOG Find center of gravity of each AAL region and print to a text file. 
%
%   Remarks:
%       If the true center of gravity is outside the region, the voxel 
%       inside the region that is closest to this point is used as the 
%       actual center of gravity. 
%
%       MNI coordinates are reported. 

%Origin of MNI space
origin=[91 126 72];

%Voxel dimension
mm=1; 

%Read in AAL template
[hdr,data]=read('aal.nii'); 

ind_aal=setdiff(unique(data(:)),0);

coor=zeros(length(ind_aal),3);
for i=1:length(ind_aal)
    [x,y,z]=ind2sub(size(data),find(ind_aal(i)==data));
    [val,ind_min]=min(sqrt((mean(x)-x).^2+(mean(y)-y).^2+(mean(z)-z).^2)); 
    coor(i,:)=[x(ind_min(1)),y(ind_min(1)),z(ind_min(1))]; 
end
%Map voxel coordinates to MNI coordinates
coor=(coor-repmat(origin,length(ind_aal),1))*mm;

%Write to a text file
dlmwrite('aalCOG.txt',coor,'delimiter',' ','precision','%.3f');