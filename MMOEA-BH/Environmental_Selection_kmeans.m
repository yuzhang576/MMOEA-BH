function [Archive, LBA, LBA_SCD, PBA, PBA_SCD, idx, maxCluster] = Environmental_Selection_kmeans(Problem, Archive, Population, LBA, PBA, n_PBA, maxCluster)
   %% construct Archive
    % Remove duplicate individuals
    N = Problem.N;
    temp_p = [Population, Archive, LBA{:}];
    temp_decs = temp_p.decs;
    [~,idx_1] = unique(temp_decs(:,1));
    Archive = temp_p(idx_1);
    
    % cluster decs
    temp_Archive = Archive;
    idx = kmeans(Archive.decs, maxCluster);    

    temp_Archive_nd = [];
    temp_idx_d_nd = [];
    for i = 1:maxCluster
        indv_t1 = temp_Archive(idx==i);
        num_ND = NDSort(indv_t1.objs,1)==1;
        temp_Archive_nd = [temp_Archive_nd,indv_t1(num_ND)];
        temp_idx_d_nd = [temp_idx_d_nd;i*ones(sum(num_ND),1)];
    end
    Archive = temp_Archive_nd;
    idx = temp_idx_d_nd;

    %% Remove points when Archive exceeds N (Method 2):
    % Find the niche with most points
    % Select the M+1 points with highest density in y space from this niche
    % Calculate their density in x space and remove the one with highest density
    while numel(Archive)>N
        num_cluster = tabulate(idx);
        [~,idx_max] = max(num_cluster(:,2));
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

    %% repair idx_d and maxCluster_d
    a1 = tabulate(idx);
    idx_t = sort(find(a1(:,2)==0),'descend')';
    for i = idx_t
        idx(idx>i,:) = idx(idx>i,:)-1;
    end
    maxCluster = max(idx);

    %% update PBA
    for i = 1:N
        [~, temp_PBA,temp_PBA_SCD] = nd_pccs_sort_2([Population(i), PBA{i}]);
        PBA{i} = temp_PBA(1:min(numel(temp_PBA),n_PBA));
        PBA_SCD{i} = temp_PBA_SCD(1:min(numel(temp_PBA),n_PBA));
    end

    % update LBA
    [idx, maxCluster] = deal(idx, maxCluster);
    for i = 1:maxCluster
        [~, LBA{i}, LBA_SCD{i}]= nd_pccs_sort_2(Archive(i==idx),ceil(N/maxCluster));
    end
    LBA(maxCluster+1:end) = [];
    LBA_SCD(maxCluster+1:end) = [];
end