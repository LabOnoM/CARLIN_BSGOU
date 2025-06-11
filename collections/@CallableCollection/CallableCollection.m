classdef (Abstract) CallableCollection < handle
    
    methods (Abstract)
        call_result = call_alleles(obj, CARLIN_def, aligned, depth);
    end
    
    methods (Static)
        [alleles, which_seq, weight_contribution] = call_alleles_exact(CARLIN_def, aligned_seqs, aligned_seq_weights, top_only);
        [alleles, which_seq, weight_contribution] = call_alleles_coarse_grain(CARLIN_def, aligned_seqs, aligned_seq_weights, top_only);
    end
    
end