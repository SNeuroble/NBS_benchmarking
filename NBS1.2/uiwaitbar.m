function h = uiwaitbar(varargin)
%uiwaitbar: A waitbar that can be embedded in a GUI figure.
% Syntax:
% POSITION = [20 20 200 20]; % Position of uiwaitbar in pixels.
% H = uiwaitbar(POSITION);
% for i = 1:100
% uiwaitbar(H,i/100)
% end
% written by Doug Schwarz, 11 December 2008

if ishandle(varargin{1})
    ax = varargin{1};
    value = varargin{2};
    p = get(ax,'Child');
    y = get(p,'YData');
    y(2:3) = value;
    set(p,'YData',y);
    return
end
 
pos = varargin{1};
H=varargin{2};
bg_color = 'w';
fg_color = 'r';
h = axes('Parent',H,'Units','normalized',...
    'Position',pos,...
    'XLim',[0 1],'YLim',[0 1],...
    'XTick',[],'YTick',[],...
    'Color',bg_color,...
    'XColor',bg_color,'YColor',bg_color);
patch([0 0 1 1],[0 0 0 0],fg_color,...
    'Parent',h,...
    'EdgeColor','none',...
    'EraseMode','none');