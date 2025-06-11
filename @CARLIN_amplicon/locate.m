function loc = locate(obj, bp)
    
    loc.abs = find(obj.bounds.ordered(:,1) <= bp & obj.bounds.ordered(:,2) >= bp);    
    assert(~isempty(loc.abs) && loc.abs >= 1 && loc.abs <= obj.N.motifs);
    loc.pos = bp-obj.bounds.ordered(loc.abs,1)+1;
    loc.rel = round(loc.abs/3);
    
    if (loc.abs==1)
        loc.type = 'prefix';
        loc.rel = 1;
    elseif (loc.abs==obj.N.motifs)
        loc.type = 'postfix';
        loc.rel = 1;
    elseif (mod(loc.abs,3)==0)
        loc.type = 'cutsites';
    elseif (mod(loc.abs,3)==2)
        loc.type = 'consites';
    elseif (mod(loc.abs,3)==1)
        loc.type = 'pams';
    end 
end