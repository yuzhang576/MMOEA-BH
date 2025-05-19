function [Archive, LBA, LBA_SCD, PBA, PBA_SCD, idx, maxCluster] = Environmental_Selection_DBSCAN(Problem, Archive, Population, LBA, PBA, n_PBA, N, maxCluster, epsilon, minpts)
   %% construct Archive by DBSCAN
    if ~exist('epsilon','var')
        [epsilon, minpts] = deal(0.15, 5);
    end
    maxCluster_old = maxCluster;
    Archive_old = Archive;

    % remove the same individuals
    temp_p = [Population, Archive, LBA{:}];
    temp_decs = temp_p.decs;
    [~,idx] = unique(temp_decs(:,1));
    Archive = temp_p(idx);
    
    % cluster decs
    temp_Archive = Archive;
    idx_1 = dbscan(map0to1(temp_Archive.decs), epsilon, minpts);
    if ismember(-1, idx_1)        maxCluster_1 = numel(unique(idx_1))-1;    else        maxCluster_1 = numel(unique(idx_1));    end
    % cluster objs
    idx_2 = dbscan(map0to1(temp_Archive.objs), epsilon, minpts);
    if ismember(-1, idx_2)        maxCluster_2 = numel(unique(idx_2))-1;    else        maxCluster_2 = numel(unique(idx_2));    end
    % detect irregularity 1
    temp_num = tabulate(idx_1);
    if ismember(-1, idx_1)        temp_num = temp_num(2:end,2);    else        temp_num = temp_num(:,2);    end
    if max(temp_num)/min(temp_num) > 5        flag_dec = 0;        else        flag_dec = 1;                end
    % detect irregularity 2
    if maxCluster_1==1 || maxCluster_2==1
        flag_one = 1;
        if maxCluster_1==1
            flag_dec_one = 1;
        else
            flag_dec_one = 0;
        end
    else
        flag_one = 0;
    end

    if flag_one
        if flag_dec_one
            [idx, maxCluster_d] = deal(idx_2, maxCluster_2);
        else
            [idx, maxCluster_d] = deal(idx_1, maxCluster_1);
        end
    else
        if maxCluster_1 >= maxCluster_2 && flag_dec
            [idx, maxCluster_d] = deal(idx_1, maxCluster_1);
        else
            [idx, maxCluster_d] = deal(idx_2, maxCluster_2);
        end
    end

    temp_Archive_nd = [];
    for i = -1:maxCluster_d
        indv_t1 = temp_Archive(idx==i);
        temp_Archive_nd = [temp_Archive_nd,indv_t1(NDSort(indv_t1.objs,1)==1)];
    end
    Archive = temp_Archive_nd;

    % Remove outliers
    % cluster decs
    idx_1 = dbscan(map0to1(Archive.decs), epsilon, minpts);
    if ismember(-1, idx_1)        maxCluster_1 = numel(unique(idx_1))-1;    else        maxCluster_1 = numel(unique(idx_1));    end
    % cluster objs
    idx_2 = dbscan(map0to1(Archive.objs), epsilon, minpts);
    if ismember(-1, idx_2)        maxCluster_2 = numel(unique(idx_2))-1;    else        maxCluster_2 = numel(unique(idx_2));    end
    if maxCluster_1 == maxCluster_2
        idx_del = any([idx_1==-1,idx_2==-1],2);
    else
        idx_del = all([idx_1==-1,idx_2==-1],2);
    end
    if sum(idx_del==1) == numel(idx_del),        idx_del = [];    end
    Archive(idx_del) = [];

    %% Remove points when Archive exceeds N
    % Find the niche with the most points
    % Select the M+1 points with highest density in y space from this niche
    % Calculate their density in x space and remove the one with highest density
    % cluster decs
    idx_1 = dbscan(map0to1(Archive.decs), epsilon, minpts);
    % cluster objs
    idx_2 = dbscan(map0to1(Archive.objs), epsilon, minpts);
    if maxCluster_1 >= maxCluster_2
        idx = idx_1;
    else
        idx = idx_2;
    end

    while numel(Archive)>N
        num_cluster = tabulate(idx);
        [~,idx_max] = max(num_cluster(:,2));
        idx_max = num_cluster(idx_max,1);
        density_y = calc_pccs(Archive.objs);
        sel_cluster = find(idx == idx_max);
        density_yt = density_y(sel_cluster);
        [~,idx1] = sort(density_yt,'descend');
        try
            idx2 = idx1(1:Problem.M+1);
        catch
            idx2 = idx1;
        end
        density_x = calc_pccs(Archive.decs);
        density_xt = density_x(sel_cluster(idx2));
        [~,idx3] = max(density_xt);
        Archive(sel_cluster(idx2(idx3))) = [];
        idx(sel_cluster(idx2(idx3))) = [];
    end
    maxCluster = numel(unique(idx));
    if ismember(-1, idx),        idx(idx==-1) = maxCluster;    end
    
    %% update PBA
    for i = 1:numel(Population)
        [~, temp_PBA,temp_PBA_SCD] = nd_pccs_sort_2([Population(i), PBA{i}]);
        PBA{i} = temp_PBA(1:min(numel(temp_PBA),n_PBA));
        PBA_SCD{i} = temp_PBA_SCD(1:min(numel(temp_PBA),n_PBA));
    end

    %% update LBA
    for i = 1:maxCluster
        [~, LBA{i}, LBA_SCD{i}]= nd_pccs_sort_2(Archive(i==idx),ceil(N/maxCluster));
    end
    LBA(maxCluster+1:end) = [];
    LBA_SCD(maxCluster+1:end) = [];
end