classdef pmethods < handle
    
    methods (Static)
        
        function [plv_hilbert, plv_morlet] = average_run_plv(corrcoeff, filter, noise, freq, num_times)
            
            plv_hilbert = [];
            plv_morlet = [];
            
            for i = 1:num_times
                plvs = plvm.run_plv(corrcoeff, filter, noise, 'both', freq, 0);
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
                r
                [plv_hilbert, plv_morlet] = plvm.pmethods.average_run_plv(r, [0,0], [0,0], [6:44], num_runs);
                Hout = [Hout; mean(plv_hilbert, 1)];
                Mout = [Mout; mean(plv_morlet, 1)];
            end
            
            figure;
            plot(Hout')
            hold on
            plot(Mout')
            
            leg = [];
            for i = corrcoeffs
                leg = [leg; sprintf('H, r=%.02f',i)];
            end
            for i = corrcoeffs
                leg = [leg; sprintf('M, r=%.02f',i)];
            end
            legend(leg)
            
        end
        
        function [Hout, Mout] = varying_corr_coeff_with_filter
            % Single narrow band signal over a variety of coefficients with
            % dual-narrow band signals the answer is the same as without
            % any filtering applied.
            % Also, try varying the hilbert width. For small widths, there
            % are multiple harmonics in the PLV that arise. For larger
            % widths, the effect is more spread out and lower over the
            % range of frequencies. (be sure to change frequency cutoffs as
            % well for estimation of PLV).
            % At approx width = 8, the morlet and hilbert give same graphs
            corrcoeffs = [0:.2:1];
            num_runs = 1;
            
            Hout = [];
            Mout = [];
            
            for r = corrcoeffs
                [plv_hilbert, plv_morlet] = plvm.pmethods.average_run_plv(r, [1,0], [0,0], [6:44], num_runs);
                Hout = [Hout; mean(plv_hilbert, 1)];
                Mout = [Mout; mean(plv_morlet, 1)];
            end
            
            figure;
            plot(Hout')
            hold on
            plot(Mout')
            
            leg = [];
            for i = corrcoeffs
                leg = [leg; sprintf('H, r=%.02f',i)];
            end
            for i = corrcoeffs
                leg = [leg; sprintf('M, r=%.02f',i)];
            end
            legend(leg)
            
            
        end
        
        function [Hout, Mout] = varying_corr_coeff_with_noise
            % Has the effect of reducing the overall PLV values.
            % For high enough noise (tested with 80) the effect disappears
            corrcoeffs = [0:.2:1];
            num_runs = 1;
            
            Hout = [];
            Mout = [];
            
            for r = corrcoeffs
                [plv_hilbert, plv_morlet] = plvm.pmethods.average_run_plv(r, [0,0], [3,3], [6:44], num_runs);
                Hout = [Hout; mean(plv_hilbert, 1)];
                Mout = [Mout; mean(plv_morlet, 1)];
            end
            
            figure;
            plot(Hout')
            hold on
            plot(Mout')
            
            leg = [];
            for i = corrcoeffs
                leg = [leg; sprintf('H, r=%.02f',i)];
            end
            for i = corrcoeffs
                leg = [leg; sprintf('M, r=%.02f',i)];
            end
            legend(leg)
            
            
        end
        
        
        function [Hout, Mout] = varying_corr_coeff_with_noise_filter
            % Has the effect of reducing the overall PLV values.
            % For high enough noise (tested with 80) the effect disappears
            corrcoeffs = [0:.2:1];
            num_runs = 1;
            
            Hout = [];
            Mout = [];
            
            for r = corrcoeffs
                [plv_hilbert, plv_morlet] = pmethods.average_run_plv(r, [1,0], [15,15], [6:44], num_runs);
                Hout = [Hout; mean(plv_hilbert, 1)];
                Mout = [Mout; mean(plv_morlet, 1)];
            end
            
            figure;
            plot(Hout')
            hold on
            plot(Mout')
            
            leg = []
            for i = corrcoeffs
                leg = [leg; sprintf('H, r=%.02f',i)]
            end
            for i = corrcoeffs
                leg = [leg; sprintf('M, r=%.02f',i)]
            end
            legend(leg)
            
            
        end
        
        function get_morlet_widths
            width=37; Fs=100; f=10; sf=f/width; dt = 1/Fs; st=1/(2*pi*sf);t=-3.5*st:dt:3.5*st;
            y = plvm.PLVM.morlet(f,t,width); pwelch(real(y),[],[],[],Fs)
        end
        
        
    end
    
end