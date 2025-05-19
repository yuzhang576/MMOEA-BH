%%
function [pbest, lbest] = REPSelection_PSO1(Problem, PBA, PBA_SCD, LBA, LBA_SCD, idx_APC, maxCluster_APC)
% Select particles as the pbest and lbest
for i = 1:numel(PBA)
    %         idx1 = TournamentSelection(2,1,-PBA_SCD{i});
    idx1 = 1;
    pbest(i) = PBA{i}(idx1);
end

if numel(idx_APC)<Problem.N
    k = ceil(Problem.N/numel(idx_APC));
    idx_APC = repmat(idx_APC,k,1);    
end

idx_APC = idx_APC(1:Problem.N);
for i = 1:maxCluster_APC
    lbest_temp = LBA{i};
    pos = i==idx_APC;
    idx = TournamentSelection(2,sum(pos),-LBA_SCD{i});
    lbest(pos) = lbest_temp(idx);
end
end