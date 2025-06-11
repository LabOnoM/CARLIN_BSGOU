function motif = find_motif_from_bp(obj, n)
    
    assert(n >= 1 && n <= obj.width.CARLIN);
    motif = find(obj.bounds.ordered(:,1) <= n & obj.bounds.ordered(:,2) >= n);
    
end