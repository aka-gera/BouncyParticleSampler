function histmat  = hist2(x, y, xedges, yedges)
% function histmat  = hist2(x, y, xedges, yedges)
%
% Extract 2D histogram data containing the number of events
% of [x , y] pairs that fall in each bin of the grid defined by 
% xedges and yedges. The edges are vectors with monotonically 
% non-decreasing values.  
%
%EXAMPLE 
%
% events = 1000000;
% x1 = sqrt(0.05)*randn(events,1)-0.5; x2 = sqrt(0.05)*randn(events,1)+0.5;
% y1 = sqrt(0.05)*randn(events,1)+0.5; y2 = sqrt(0.05)*randn(events,1)-0.5;
% x= [x1;x2]; y = [y1;y2];
%
%For linearly spaced edges:
% xedges = linspace(-1,1,64); yedges = linspace(-1,1,64);
% histmat = hist2(x, y, xedges, yedges);
% figure; pcolor(xedges,yedges,histmat'); colorbar ; axis square tight ;
%
%For nonlinearly spaced edges:
% xedges_ = logspace(0,log10(3),64)-2; yedges_ = linspace(-1,1,64);
% histmat_ = hist2(x, y, xedges_, yedges_);
% figure; pcolor(xedges_,yedges_,histmat_'); colorbar ; axis square tight ;
% University of Debrecen, PET Center/Laszlo Balkay 2006
% email: balkay@pet.dote.hu
if nargin ~= 4
    error ('The four input arguments are required!');
    return;
end
if any(size(x) ~= size(y)) 
    error ('The size of the two first input vectors should be same!');
    return;
end
[~, xbin] = histc(x,xedges);
[~, ybin] = histc(y,yedges);
%xbin, ybin zero for out of range values 
% (see the help of histc) force this event to the 
% first bins
xbin(xbin == 0) = inf;
ybin(ybin == 0) = inf;
xnbin = length(xedges);
ynbin = length(yedges);
if xnbin >= ynbin
    xy = ybin*(xnbin) + xbin;
      indexshift =  xnbin; 
else
    xy = xbin*(ynbin) + ybin;
      indexshift =  ynbin; 
end
%[xyuni, m, n] = unique(xy);
xyuni = unique(xy);
xyuni(end) = []; 
hstres = histc(xy,xyuni);
clear xy;
histmat = zeros(ynbin,xnbin);
histmat(xyuni-indexshift) = hstres;
