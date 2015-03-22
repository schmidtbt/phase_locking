classdef pmethods < handle
    
    methods (Static)
        
        function [plv_hilbert, plv_morlet] = average_run_plv(corrcoeff, filter, noise, freq, num_times)
            
            plv_hilbert = [];
            plv_morlet = [];
            
            for i = 1:num_times
                plvs = run_plv(corrcoeff, filter, noise, 'both', freq, 0);
                plv_hilbert = [plv_hilbert; plvs(1,:)];
                plv_morlet = [plv_morlet; plvs(2,:)];
            end
            
        end
        
        function [Hout, Mout] = varying_corr_coeff
            corrcoeffs = [0:.1:1];
            num_runs = 1;
            
            Hout = [];
            Mout = [];
            
            for r = corrcoeffs
                [plv_hilbert, plv_morlet] = pmethods.average_run_plv(r, [0,0], [0,0], [6:44], num_runs);
                Hout = [Hout; mean(plv_hilbert, 1)];
                Mout = [Mout; mean(plv_morlet, 1)];
            end
            
            figure;
            plot(Hout')
            hold on
            plot(Mout')
            
        end
        
    end
    
end