function draw_atlas_boundaries(mat,varargin)
% inputs: mat,atlascategory
% options: saveimg,reorderimg,numatlasnodes,datatype
% ** use datacategory to specify setting colorscales for 'ICC' or 'scandur'
%
% example usage
%
% % get scan duration to fair good exc
% cat=[0.4, 0.6, 0.74];
% for j=1:51
% iccsumm{j}=domatrixsummary_avg(icc{2,6}{j},'atlascategory','subnetwork','datacategory','ICC','saveimg',0,'reorderimg',1); end
% for j=1:3
% tmp=zeros(size(iccsumm{1}));
% for i=1:51
% tmp=tmp++(iccsumm{i}>cat(j)); end
% subnet_scandur{j}=tmp; end
%
% draw_atlas_boundaries(subnet_scandur{1},'subnetwork')

%% Parse input
p = inputParser;

defaultsaveimg=0;
defaultreorderimg=0;
defaultatlascategory='subnetwork'; % options include subnetwork, lobe
defaultdatacategory='none'; % options include ICC, scandur, none
defaultusingsummary=0;

addParameter(p,'saveimg',defaultsaveimg,@isnumeric);
addParameter(p,'reorderimg',defaultreorderimg,@isnumeric);
addParameter(p,'atlascategory',defaultatlascategory,@ischar);
addParameter(p,'datacategory',defaultdatacategory,@ischar);
addParameter(p,'usingsummary',defaultusingsummary,@isnumeric);

parse(p,varargin{:});

saveimg = p.Results.saveimg;
reorderimg = p.Results.reorderimg;
atlascategory = p.Results.atlascategory;
datacategory = p.Results.datacategory;
usingsummary = p.Results.usingsummary;

clearvars p varargin

%% set up matrix

% restruct and reorder if needed
if size(mat,1) ~= size(mat,2)
    mat=structure_data(mat);
    mat=reorder_matrix_by_atlas(mat,atlascategory);
    if ~reorderimg
        warning('reorderimg=0 but doing reordering anyways because input is vector.')
    end
elseif reorderimg
    mat=reorder_matrix_by_atlas(mat,atlascategory);
end
matdim=size(mat,1); matdim_new=matdim;
numatlasnodes=matdim;

% if need to split in half:
%     matdim_new=matdim/2;
%     x=[1 matdim_new matdim_new+1 matdim];
%     mat=mat(x(1):x(2),x(1):x(2))+mat(x(3):x(4),x(1):x(2))+mat(x(1):x(2),x(3):x(4))+mat(x(3):x(4),x(3):x(4));
%     mat=mat/4;

% load atlas 
if ~usingsummary 
    map=load_atlas_mapping(numatlasnodes,atlascategory);
else % load a default )default=268 nodes)
    map=load_atlas_mapping(268,atlascategory);
end

%% set up labels and divisions
% assuming correct size

if usingsummary
    lobe_mapping=1:matdim;
    xticks=1:matdim;
else
    lobe_mapping=map.category(2:end)-map.category(1:end-1);
    lobe_mapping=[0; find(lobe_mapping==1); length(lobe_mapping)+1];
    lobe_mapping=lobe_mapping+1;
    xticks=ceil(lobe_mapping(1:end-1)+(lobe_mapping(2:end)-lobe_mapping(1:end-1))/2);
end

for i=1:length(unique(map.label))
    xticklbls(i)=unique(map.label(map.category==i));
end

%% set up colors
% print min thresh as dark blue and max as orange (poor=0 for other ICC plots)
%previously: if max(max(mat)) < 1; reliabilitymap=1;
cmap='parula';
% can specify colormap(bipolar(100,0.2,'cubic'))
setcolors;

%% draw image and lines
figure
colormap(cmap)
if strcmp(datacategory,'none')
    imagesc(mat)
elseif strcmp(datacategory,'ICC') || strcmp(datacategory,'icc')
    image((mat-minthresh)*imscale)
end
hold on

if ~usingsummary
    for i=2:length(lobe_mapping)-1
        % plot horizontals
        p1 = [lobe_mapping(i)-0.5,0];
        p2 = [lobe_mapping(i)-0.5,matdim_new];
        m=plot([p1(2),p2(2)],[p1(1),p2(1)],'w-','LineWidth',2);
        m.Color(4) = 0.7; % semi-transparent
        
        % plot verticals
        p1 = [0,lobe_mapping(i)-0.5];
        p2 = [matdim_new,lobe_mapping(i)-0.5];
        m=plot([p1(2),p2(2)],[p1(1),p2(1)],'w-','LineWidth',2);
        m.Color(4) = 0.7; % semi-transparent
    end
end

%% adjust axes and colorbar

set(gca,'XTick',xticks);
set(gca,'XTickLabel',xticklbls)
set(gca,'XTickLabelRotation',45)
set(gca,'YTick',xticks);
set(gca,'YTickLabel',xticklbls)

cb=colorbar;
if ~ strcmp(datacategory,'none')
    set(get(cb,'Title'),'String',ctitle)
    set(cb,'YTick',cticks); set(cb,'YTickLabel',cticklbls)
    set(cb, 'ylim', cylim)
end

hold off
axis('square')

%% save
if saveimg
    clockinfo{1}=datestr(now, 'HH'); clockinfo{2}=datestr(now, 'MM'); clockinfo{3}=datestr(now, 'SS');
    timeID=sprintf('%sH%sM%sS',clockinfo{1},clockinfo{2},clockinfo{3});
    savefig(sprintf('matimg_%s.fig',timeID));
    saveas(gcf,sprintf('matimg_%s.png',timeID))
    close(gcf)
end


