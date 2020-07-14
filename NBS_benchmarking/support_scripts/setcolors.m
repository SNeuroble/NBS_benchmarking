% set threshold parameters for colormap
% first, set 'datacategory' variable in environment ('ICC' or 'none')

if ~ strcmp(datacategory,'none')
    imscale=54;
    imscale=75; % manually changed for 0-1
    if strcmp(datacategory,'ICC') | strcmp(datacategory,'icc') % reliability map (0-1)
	% scales such that min is 0.1 and max is 0.74 (threshold of excellent) to maximize color gradations
        minthresh=0.14; % user-defined
        maxthresh=.6; % user-defined

%         minthresh=0.4; % user-defined
%         maxthresh=.74; % user-defined
        
        ctitle='ICC';
        cmin=0.1;
        cmax=0.75;
        
%         % manually changed for simple 0-1 plots
%         minthresh=0.3;
%         maxthresh=1;
%         cmin=0.3;
%         cmax=1;
        
        inc=0.1;
        cticklbls=[cmin:inc:cmax];
        
    elseif strcmp(datacategory,'scandur') % scan duration in runs (0-6)
        minthresh=1; % user-defined
        maxthresh=4.2; % user-defined
        unitinc=0.1; inc=1/unitinc; % user-defined
        maxthresh=maxthresh*inc;
        
        ctitle='min';
        cmin=0;
        cmax=6;
        
        runlength=6;
        cticklbls=[cmin:1:cmax]*runlength+runlength;
        
        cticklbls=sort(cticklbls,'descend');
        cticklbls=num2cell(cticklbls);
        cticklbls{1}='>36';
        
        cmax=cmax*inc;
        
        
    end
    
    imscale=imscale/(maxthresh-minthresh);
    cticks=([cmin:inc:cmax]-minthresh)*imscale;
    cylim=([cmin, cmax]-minthresh)*imscale;
end
