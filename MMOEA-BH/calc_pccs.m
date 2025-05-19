function varargout = calc_pccs(obj)
    K = size(obj,1);
    if K >= 2
        fmax = max(obj);
        fmin = min(obj);
        L = ceil(K * (obj-fmin)./(fmax-fmin));
        L(L==0) = L(L==0) + 1;
        PCD = pdist2(L,L,'cityblock');
        PCD(logical(eye(K))) = [];
        PCD = reshape(PCD,K-1,K)';
        PCD(PCD==0) = PCD(PCD==0) + 0.5;
        density = sum(1./PCD.^2,2);
        varargout = {density};
    else
        varargout = {1};
    end
end