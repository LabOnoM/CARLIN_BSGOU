function [sc, al] = cas9_align(obj, seq)

    assert(~isempty(seq), 'Sequence to align is empty');
    assert(~any(isgap(seq)), 'Input sequence should not have any insertions');
    
    refseq        = obj.seq.CARLIN;
    open_penalty  = obj.open_penalty;
    close_penalty = obj.close_penalty;
    match         = obj.sub_score;
    
    [sc, al] = CARLIN_amplicon.cas9_align_mex([0 nt2int(seq)], [0 nt2int(refseq)], open_penalty, [0 close_penalty], padarray(match, [1 1], 'pre'));        
    al(al==0)=16;
    al = int2nt(al);
    al_seq = al(1,:);
    al_ref = al(2,:);
   
    assert(isequal(degap(al_seq), seq), '%s\n%s', degap(al_seq), seq);
    assert(isequal(degap(al_ref), refseq), '%s\n%s', degap(al_ref), refseq);
    assert(length(al_seq) == length(al_ref));
    assert(~any(isgap(al_seq) & isgap(al_ref)));
    
    al = obj.desemble_sequence(al_seq, al_ref);
    
end