function site = coarse_grain_motif(obj, n)

    assert(n >= 1 && n <= obj.N.motifs);
    site = max(1,ceil((n-1)/3));
    
end