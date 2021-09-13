function displayColors(rgb,scale,header);
%
% displayColors(rgb,[scale]);
%
% Displays a matrix of n color tiles.
% rgb:     3xn matrix where the columns define
%          the rgb gun values of the n tiles.
%          Tile numbers increase row-wise in the image.
%
% scale:   scale factor.  True rgb values is rgb/scale.
%          default value is 1;
%
%header:   title of image.  


%6/7/96	gmb	Wrote it.

if nargin<2	
	scale=1;
end

ncolors=size(rgb,2);

%determine the number of columns and rows
nrows=floor(sqrt(ncolors));
ncols=ceil(ncolors/nrows);

%set up image parameters
boxsize=10;
border=2;

%size of image
m=nrows*(boxsize+border)+border;
n=ncols*(boxsize+border)+border;


%position of checks
[x,y]=meshgrid((border+1):(boxsize+border):n,(border+1):(boxsize+border):m);

%make boxes
img=ncolors+1*ones(m,n);
cur_color=0;
for i=1:nrows
	for j=1:ncols
		cur_color=cur_color+1;
		if cur_color<=ncolors
			img(y(i,j):y(i,j)+boxsize-1,x(i,j):x(i,j)+boxsize-1)=cur_color*ones(boxsize);
		end	
	end
end

%set up color map

cmap=[rgb'/scale;0,0,0];

bad=find(cmap'<0 | cmap' >1);



if length(bad>0)
	disp(['Warning, image(s) ',sprintf('%3.0f',(bad+1)/3),' out of range.']);
end
cmap(cmap<0)=zeros(size(find(cmap<0)));
cmap(cmap>1)=ones(size(find(cmap>1)));

clf
hold off
image(img);
% x=x';y=y';
% for i=1:length(bad)
% 	text(x(round((bad(i)+1)/3)+boxsize/2),y(round((bad(i)+1)/3)+boxsize/2),'X','Color','k');
% end


%set(gca,'AspectRatio',[1,1])
if (nargin>2)
	title(header);
	set(gca,'Position',[0.05,0.00,0.95,0.90]);
else
	set(gca,'Position',[0.05,0.00,0.95,0.95]);
end
colormap(cmap);
axis off

