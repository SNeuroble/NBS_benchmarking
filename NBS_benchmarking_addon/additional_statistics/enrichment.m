% Enrichment analysis
%http://software.broadinstitute.org/gsea/doc/subramanian_tamayo_gsea_pnas.pdf

% 1. order edges by ?differential expression?
% rank by correlation of gene with class label
rs = -[-10:-1];
set_id=logical([ones(5,1); zeros(5,1)]); % map from edge to network (1 = in-network, 0 = out-of-network)

% 2. calculate enrichment score
%calculate weighted sum of correlations - positive in class,  negative out of class
P(1)=0;
for i=2:length(rs)+1
    P(i) = P_hit-P_miss;
    P(i) = abs(rs(i)) + P(i-1);
end

ES = max(P);

r_sign=rs>0;
isgs=r_sign && set_id;




% separate calculations for each GeneSet
out = zeros(GS.nb,9); NES_perm = zeros(GS.nb,opts.perm_nb);
LEF = cell(GS.nb,1);
if opts.show
    [~,~] = mkdir('GSEA_plots');
    fprintf('Progress:\n');
    fprintf(['\n' repmat('.',1,GS.nb) '\n\n']);
end
parfor a=1:GS.nb 
    opts_tmp = opts;
    out_tmp = zeros(1,9);
    GS_tmp = GS;
        
    %sort ranks in GS to find LEF
    RI_tmp = sort(RI{a},opts_tmp.sort_type); 
    
    %Enrichment Score calculation
    Phit = zeros(n,1); 
    Pmiss = ones(n,1) * 1/(n-GS_tmp.entrez_nb(a));
    Phit(ind_d{a}) = (abs(RI{a}).^opts_tmp.p)/sum(abs(RI{a}));
    Pmiss(ind_d{a}) = 0;
    ES = cumsum(Phit - Pmiss);    
    [~,ind] = max(abs(ES));     %find maximum ES 
    ESmax = ES(ind);

    % Find p-value by permutation test
    ES_perm = ES_perm_test(ranks_perm,data_entrez,GS_tmp.entrez{a},opts_tmp.perm_nb,opts_tmp.sort_type);
        
    out_tmp(1) = a;        % GS name  
    out_tmp(2) = GS_tmp.entrez_nb(a);     % GS Size  
    out_tmp(3) = ESmax;           % Enrichment Score 
    
    if ESmax>0    %positive ES
        [LEF_tmp,ind_LEF] = intersect_fast(GS_tmp.entrez{a},data_entrez(1:ind));   % genes from LEF 
        LEF_tmp(:,2) = RI_tmp(ind_LEF);      % ES from genes of LEF   
        ind = ES_perm > 0;
        tmp_NES = zeros(1,opts.perm_nb);
        tmp_NES(ind) = ES_perm(ind)/mean(ES_perm(ind));
        NES_perm(a,:) = tmp_NES;
        if ~isempty(ES_perm(ind))
            out_tmp(4) = max(1,sum(ES_perm >= ESmax))/length(ES_perm(ind));  % p-value
            out_tmp(5) = ESmax/mean(ES_perm(ind)); % Normalized Enrichment Score 
        else
            out_tmp(4) = NaN;  % p-value
            out_tmp(5) = NaN; % Normalized Enrichment Score 
        end
        
    else    %negative ES
        [LEF_tmp,ind_LEF] = intersect_fast(GS_tmp.entrez{a},data_entrez(ind:end));   % genes from LEF 
        LEF_tmp(:,2) = RI_tmp(ind_LEF);      % ES from genes of LEF   
        ind = ES_perm < 0;
        tmp_NES = zeros(1,opts.perm_nb);
        tmp_NES(ind) = ES_perm(ind)/abs(mean(ES_perm(ind)));
        NES_perm(a,:) = tmp_NES;
        if ~isempty(ES_perm(ind))
            out_tmp(4) = max(1,sum(ES_perm <= ESmax))/length(ES_perm(ind));  % p-value
            out_tmp(5) = ESmax/abs(mean(ES_perm(ind))); % Normalized Enrichment Score 
        else
            out_tmp(4) = NaN;  % p-value
            out_tmp(5) = NaN; % Normalized Enrichment Score 
        end
    end
    out_tmp(9) = size(LEF_tmp,1);          % observed counts in LEF
    
    out(a,:) = out_tmp;
    LEF{a} = LEF_tmp;
    if opts.show
        fprintf('\b|\n');
    end
end



% 3. estimate significance

% 4. multiple comparison correction

% may also want to try hypergeometric:
% https://www.mathworks.com/help/stats/hygecdf.html