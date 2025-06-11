classdef (Abstract, Sealed) CARLIN_viz < handle
        
    %% Constant Definitions
    properties (Constant)        
        color = containers.Map({'N'; 'E'; 'M'; 'I'; 'D'}, ...
                               {[255, 255, 255]; [255, 220, 220]; [0, 255, 0]; [0, 0, 255]; [255, 0, 0]});        
        alpha = struct('prefix', 0.5, ...
                       'consite', 0.7, ...
                       'cutsite', 1.0, ...
                       'pam', 0.9, ...
                       'postfix', 0.5, ...
                       'overlay', 0.4);
        
    end 
    
    methods (Static)

        function alpha_val = get_alpha(CARLIN_def)

            alpha_val = zeros(1,CARLIN_def.width.CARLIN);
            alpha_val(CARLIN_def.bps.prefix)  = CARLIN_viz.alpha.prefix;
            alpha_val(CARLIN_def.bps.consite) = CARLIN_viz.alpha.consite;
            alpha_val(CARLIN_def.bps.cutsite) = CARLIN_viz.alpha.cutsite;
            alpha_val(CARLIN_def.bps.pam)     = CARLIN_viz.alpha.pam;
            alpha_val(CARLIN_def.bps.postfix) = CARLIN_viz.alpha.postfix;

        end
    
    end
           
end
