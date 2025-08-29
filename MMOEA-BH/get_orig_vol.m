function [orig_vol, center] = get_orig_vol(Problem, maxCluster, LBA)
    for i = 1:maxCluster
        pop = LBA{i};
        pop_dec = pop.decs;
        center(i,:) = mean(pop_dec,1);
    end
    grid_num = ceil(power(maxCluster,1/Problem.D));
    len1 = (Problem.upper-Problem.lower)/grid_num;
    % Calculate the area (or volume) of the original rectangle
    orig_vol = prod(len1);
end