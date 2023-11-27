import bisect
import os
import numpy as np
from scipy.interpolate import interp1d
from spec_tools import Spectrum
from turbospec_wrapper.turbospec_tools import SpectrumTS

class Abundometer:
    def __init__(self, turbosynth_obj, input_obs_spec, element, line,
             conv_profiles, conv_broad,
             delta_abunds=[-0.6, -0.3, 0.0, 0.3, 0.6, ],
             synth_range=None, delta_lambda=0.01):
        
        self.element = element
        self.line = line
        self.delta_abunds = delta_abunds
        self.delta_lambda = delta_lambda
        self.conv_profiles = conv_profiles
        self.conv_broad = conv_broad
        self.turbosynth_obj = turbosynth_obj

        if synth_range is None:
            synth_range = line/2000. + 1

        self.wave_range = [line - synth_range - 1, line + synth_range + 1]
        self.fit_range = [line - synth_range, line + synth_range]

        selection = (input_obs_spec.wavelength > self.wave_range[0]) & \
                    (input_obs_spec.wavelength < self.wave_range[1])
        obs_waves = input_obs_spec.wavelength[selection]
        obs_flux = input_obs_spec.flux[selection]
        if hasattr(input_obs_spec, 'variance'):
            obs_var = input_obs_spec.variance[selection]
        else:
            obs_var = 0.1**2.

        self.obs_spec = Spectrum(obs_waves, obs_flux, variance=obs_var)
        
    def fit(self, run_renorm=True, run_selrange=True, add_shift=False, **kwargs):
        renorm_kw = ['input_snr', 'npts_cont_thres']
        selrange_kw = ['obs_n_null1stder_limit', 'syn_n_null1stder_limit', 'fluxdiff_limit']
        
        renorm_dict = {kw : kwargs[kw] for kw in renorm_kw if kw in kwargs}
        selrange_dict = {kw : kwargs[kw] for kw in selrange_kw if kw in kwargs}
        
        self.synthesize()
        if run_renorm:
            self.renormalize(**renorm_dict)
        self.determine_fit_range(run_selrange=run_selrange, **selrange_dict)
        if add_shift:
            self.shift_obs_wavelength()
        else:
            self.shift = 0.
        self.measure_abund()
        
    
    def synthesize(self, delta_abunds=None, conv_profiles=None, conv_broad=None):
        
        if delta_abunds is not None:
            self.delta_abunds = delta_abunds
        if conv_profiles is not None:
            self.conv_profiles = conv_profiles
        if conv_broad is not None:
            self.conv_broad = conv_broad
    
        synth_abunds, synths = self._calculate_synths(self.turbosynth_obj, 
                                                      self.obs_spec, self.element, 
                                                      self.line, self.conv_profiles, 
                                                      self.conv_broad, self.delta_abunds,
                                                      self.wave_range, self.delta_lambda)

        synths, avg_synth = self._make_avg_synth(synths)
            
        self.synth_abunds = synth_abunds
        self.synths = synths
        self.avg_synth = avg_synth
        
    def renormalize(self, input_snr=None, npts_cont_thres=6.):
        self.snr = input_snr
        
        obs_spec, continuum, cont_coeff, snr = self._renormalization(self.obs_spec, 
                                                                     self.avg_synth, 
                                                                     self.line, input_snr, 
                                                                     npts_cont_thres)
        
        self.obs_spec = obs_spec
        self.continuum = continuum
        self.cont_coeff = cont_coeff
        self.snr = snr
        
    def determine_fit_range(self, run_selrange=True, obs_n_null1stder_limit=3, syn_n_null1stder_limit=1,
                            fluxdiff_limit=0.05):
        fit_range, shift, pseudo_line_depth_obs = self._calc_fit_range(self.obs_spec, self.avg_synth, 
                                                self.synths, self.line, self.fit_range,
                                                obs_n_null1stder_limit=obs_n_null1stder_limit,
                                                syn_n_null1stder_limit=syn_n_null1stder_limit)
        if run_selrange:
            self.fit_range = fit_range
        self.shift = shift
        self.pseudo_line_depth_obs = pseudo_line_depth_obs
        
    def shift_obs_wavelength(self):
        
        obs_spec, synths = self._shift_and_reinterpolate(self.obs_spec, self.synths, 
                                                         self.line, self.shift)
        self.obs_spec = obs_spec
        self.synths = synths
    
    def measure_abund(self):
        
        fit_abund_chi2, fit_min_chi2, synth_chi2s = self._fit_chi2(self.obs_spec, self.synths, 
                                                                   self.element, self.line, 
                                                                   self.synth_abunds, self.fit_range)
        
        self.fit_abund_chi2 = fit_abund_chi2
        self.fit_min_chi2 = fit_min_chi2
        self.synth_chi2s = synth_chi2s

        eqw, eqw_err, eqw_abund, eqw_abund_err, eqw_5sig_lim, eqw_lim_flag = self._fit_eqw(
            self.obs_spec, self.synths, self.element, self.line, 
            self.synth_abunds, self.fit_range, self.snr)

        self.eqw = eqw
        self.eqw_err = eqw_err
        self.fit_abund_eqw = eqw_abund
        self.fit_abund_eqw_err = eqw_abund_err
        self.eqw_limit = eqw_5sig_lim
        self.eqw_lim_flag = eqw_lim_flag


    def _calculate_synths(self, turbosynth_obj, obs_spec, element, line,
                         conv_profiles, conv_broad, abunds,
                         wave_range, delta_lambda):
        '''
        TODO:  Fill out the docstring

        assuming the obs_spec is a Spectrum object and that it is RV corrected

        Need obs_spec to have matching air/vac wavelengths to the used line
            lists

        TODO:  Need to define a fit range that is smaller than the synth range, so
            that it just fits the line/feature of interest

        TODO:  Need to do renormalization
        '''

        elemabund = turbosynth_obj.get_abund(element)

        synths = []
        chi2_arr = []

        synth_abund = np.array([elemabund+delta_xfe for delta_xfe in abunds])


        for delta_xfe in abunds:
            turbosynth_obj.set_abund(element, elemabund+delta_xfe)
            # TODO: Might want to make an option to run a single atmosphere
            # with the default abundances which could speed things up and keep
            # the code from breaking on random abundances
            turbosynth_obj.make_atmosphere(star='tmp')

            turbosynth_obj.synth(wave_range, synth_fname='tmp.spec',
                       delta_lambda=delta_lambda)
            os.remove(turbosynth_obj.opac_filename)
            os.remove(turbosynth_obj.opac_filename+'.mod')

            for i, (profile, kwargs) in enumerate(zip(conv_profiles,
                                                      conv_broad)):
                if i == 0:
                    fname = turbosynth_obj.synth_fname
                else:
                    fname = convol
                convol = turbosynth_obj.convol(synth_fname=fname, profile=profile,
                                     **kwargs, convol_fname='tmp.spec.convol')

            spectrum = SpectrumTS(f'{turbosynth_obj.path}/tmp.spec.convol')
            spectrum.interpolate_spectrum(obs_spec.wavelength)

            synths += [spectrum]

            turbosynth_obj.set_abund(element, elemabund)

        return synth_abund, synths

    
    def _make_avg_synth(self, synth_specs):

        # Rescale syntheses and calculate avg synthetic spectrum
        synth_flux_arr = np.zeros((len(synth_specs[0].wavelength),len(synth_specs)))
        for i in range(len(synth_specs)):
            synth_flux_arr[:,i] = synth_specs[i].flux
        stdev_synth = np.std(synth_flux_arr, axis=1)

        if stdev_synth.min() > 0.005:
            rescale_value = synth_specs[0].flux.max()
            for i, synth in enumerate(synth_specs):
                synth.flux /= rescale_value
                synth_flux_arr[:,i] = synth.flux

        avg_flux = np.mean(synth_flux_arr, axis=1)
        stdev_synth = np.std(synth_flux_arr, axis=1)
        spread_synth = np.std(synth_flux_arr - np.tile(synth_flux_arr[:,0], (len(synth_specs),1)).T, axis=1)
        avg_synth = Spectrum(synth_specs[0].wavelength, avg_flux, variance=spread_synth)

        return synth_specs, avg_synth
    
    def _renormalization(self, obs_spec, avg_synth, line, input_snr=None, npts_cont_thres=6.):

        #Estimate rough SNR
        obs_wave_temp = np.array(obs_spec.wavelength)
        obs_flux_temp = np.array(obs_spec.flux)
        obs_weight_temp = np.array(1./np.sqrt(obs_spec.variance))
        synth_wave_temp = np.array(avg_synth.wavelength)
        synth_flux_temp = np.array(avg_synth.flux)
        
        #iteratively sigma clip obs and synth
        for i in range(10):
            obs_fit = np.polyfit(obs_wave_temp, obs_flux_temp, 1, w=obs_weight_temp)
            synth_fit = np.polyfit(synth_wave_temp, synth_flux_temp, 1)
            obs_synth_diff = (obs_flux_temp/np.polyval(obs_fit,
                obs_wave_temp)) - (synth_flux_temp/np.polyval(synth_fit, synth_wave_temp))
            obs_synth_stdev = np.std(obs_synth_diff)
            selection = (obs_synth_diff < 2.5*obs_synth_stdev) & (obs_synth_diff > -2.5*obs_synth_stdev)

            if len(obs_wave_temp[selection]) < 5:
                break
            else:
                obs_wave_temp = obs_wave_temp[selection]
                obs_flux_temp = obs_flux_temp[selection]
                obs_weight_temp = obs_weight_temp[selection]
                synth_wave_temp = synth_wave_temp[selection]
                synth_flux_temp = synth_flux_temp[selection]

        os_mean_ratio = np.mean(obs_flux_temp)/np.mean(synth_flux_temp)
        os_noise_estimate = np.std(obs_flux_temp - (synth_flux_temp * os_mean_ratio))
        snr_init_estimate = os_mean_ratio/os_noise_estimate

        # Estimate SNR from continuum
        max_flux = 0.95
        for i in range(100):
            continuum_selection = avg_synth.flux > max_flux
            x_cont = np.array(avg_synth.wavelength[continuum_selection])
            y_cont = np.array(avg_synth.flux[continuum_selection])
            max_flux -= 0.01
            if len(x_cont) > 20:
                break

        for i in range(10):
            xo_clean = np.array(obs_spec.wavelength[continuum_selection])
            yo_clean = np.array(obs_spec.flux[continuum_selection])
            wo_clean = np.array(1/np.sqrt(obs_spec.variance[continuum_selection]))

            obs_cont_fit = np.polyfit(xo_clean, yo_clean, 1, w=wo_clean)
            yo_clean_temp = yo_clean/np.polyval(obs_cont_fit, xo_clean)
            sigma = np.std(yo_clean_temp)

            yo_temp = obs_spec.flux/np.polyval(obs_cont_fit, obs_spec.wavelength)
            continuum_selection &= (yo_temp > 1-1.8*sigma) & (yo_temp < 1+2.6*sigma)

        yo_ysyn_ratio = obs_spec.flux[continuum_selection]/avg_synth.flux[continuum_selection]
        mean_cont = np.mean(yo_ysyn_ratio)
        stdev_cont = np.std(yo_ysyn_ratio)
        snr_cont_estimate = mean_cont/stdev_cont
        if input_snr is not None:
            snr = input_snr
        elif (snr_cont_estimate > snr_init_estimate) & (len(yo_ysyn_ratio) > 10):
            snr = snr_cont_estimate
        else:
            snr = snr_init_estimate

        # note that instead of carrying around x and xo_clean like in BACCHUS i'm carrying around
        # the continuum pixel mask so that it's easy to go back and select points as desired
        # until now because I think we are going to heavily use xo_clean and yo_clean

        xo_clean = obs_spec.wavelength[continuum_selection]
        yo_clean = obs_spec.flux[continuum_selection]
        ysyn_clean = avg_synth.flux[continuum_selection]

        if len(xo_clean) > 0:
            lambda_index = bisect.bisect_left(xo_clean, line) # thomas calls this i_lambda

        xp = []
        yp = []
        xm = []
        ym = []
        ysynp = []
        ysynm = []


        for im, ip in zip(range(lambda_index), range(lambda_index, len(xo_clean))):
            xp = xp + [xo_clean[ip]]
            yp = yp + [yo_clean[ip]]
            ysynp = ysynp + [ysyn_clean[ip]]
            xm = [xo_clean[im]] + xm
            ym = [yo_clean[im]] + ym
            ysynm = [ysyn_clean[im]] + ysynm


        x_cont = np.array(xm+xp)
        y_cont = np.array(ym+yp)
        ysyn_cont = np.array(ysynm+ysynp)


        # Fit continuum points
        if len(x_cont) < npts_cont_thres:
            # not enough points to fit a continuum
            init_sort = np.argsort(avg_synth.flux)
            xo_temp = np.array(obs_spec.wavelength)[init_sort]
            yo_temp = np.array(obs_spec.flux)[init_sort]
            avg_flux_temp = np.array(avg_synth.flux)[init_sort]

            selec_avgflux = avg_flux_temp[-10:]
            selec_xotemp = xo_temp[-10:]
            selec_yotemp = yo_temp[-10:]

            second_sort = np.argsort(selec_yotemp)

            selec_avgflux2 = selec_avgflux[second_sort][-7:-3]
            selec_xotemp2 = selec_xotemp[second_sort][-7:-3]
            selec_yotemp2 = selec_yotemp[second_sort][-7:-3]

            a_fit_obs = 0.
            b_fit_obs = np.mean(selec_yotemp2)/np.mean(selec_avgflux2)

            x_cont = selec_xotemp2
            y_cont = selec_yotemp2
        else:
            mean_cont = np.mean(y_cont/ysyn_cont)
            continuum_fit = np.polyfit(x_cont - line, y_cont/ysyn_cont, 1)
            rms = np.std((y_cont/ysyn_cont) - np.polyval(continuum_fit, x_cont-line))

            if rms > (mean_cont/snr):
                a_fit_obs = 0
                b_fit_obs = mean_cont
            else:
                a_fit_obs = continuum_fit[0]
                b_fit_obs = continuum_fit[1]

        # Normalize observed spectrum
        normalization_factor = a_fit_obs*(obs_spec.wavelength-line)+b_fit_obs
        yo_norm = obs_spec.flux / normalization_factor
        if len(obs_spec.variance) > 1:
            yo_var_norm = obs_spec.variance/normalization_factor**2.
            var_snr = 1.0/np.mean(np.sqrt(yo_var_norm))
            if var_snr > 0.5:
                snr = var_snr
        else:
            yo_var_norm = obs_spec.variance

        normed_spectrum = Spectrum(obs_spec.wavelength, yo_norm, variance=yo_var_norm)

        continuum = Spectrum(x_cont, y_cont)

        return normed_spectrum, continuum, [a_fit_obs, b_fit_obs], snr


    def _calc_fit_range(self, obs_spec, avg_synth, synths, line, fit_range, 
                        syn_n_null1stder_limit=1, obs_n_null1stder_limit=3,
                        fluxdiff_limit=0.05):
        # the variance in avg synth gives the stdev of synths-synth1 (the lowest abund synth)

        # calculate the 1st derivative of the observations and synth
        obs_1stderiv = np.gradient(obs_spec.flux, obs_spec.wavelength)
        synth_1stderiv = np.gradient(avg_synth.flux, avg_synth.wavelength)

        lambda_index = bisect.bisect_left(obs_spec.wavelength, line)
        lambda_index_obs = lambda_index
        ip = lambda_index
        im = lambda_index - 1
        n_attempts = 0

        # NOTE:  I'm using techniques used in BACCHUS but I think it might be better to 
        # have an ip im window that moves with lambda_index_obs rather than just continues
        # to grow which I think could cause more unpredictable behaviour jumping between
        # local minima
        while ((ip < len(obs_spec.wavelength) -1) & (im > 0) & 
               ((obs_spec.flux[lambda_index_obs] > obs_spec.flux[im]) | 
                (obs_spec.flux[lambda_index_obs] >= obs_spec.flux[ip]) | n_attempts <= 3)):

            if (obs_spec.flux[lambda_index_obs] > obs_spec.flux[im]):
                n_attempts = 0
                lambda_index_obs = im
            elif (obs_spec.flux[lambda_index_obs] > obs_spec.flux[ip]):
                n_attempts = 0
                lambda_index_obs = ip
            else:
                n_attempts += 1

            ip += 1
            im -= 1

        lambda_center_obs = obs_spec.wavelength[lambda_index_obs]

        # select significant points in the obs spectrum
        ip = lambda_index_obs
        xobs_chi2 = []
        yo_signif_selec_norm = []
        n_null1stder_obsp = 0
        for ip in range(lambda_index_obs, len(obs_spec.wavelength)):
            xobs_chi2 += [obs_spec.wavelength[ip]]
            yo_signif_selec_norm += [obs_spec.flux[ip]]
            if obs_1stderiv[ip] < 0:
                n_null1stder_obsp += 1
            if n_null1stder_obsp > obs_n_null1stder_limit:
                break

        n_null1stder_obsm = 0
        for im in range(lambda_index_obs-1, -1, -1):
            xobs_chi2 = [obs_spec.wavelength[im]] + xobs_chi2
            yo_signif_selec_norm = [obs_spec.flux[im]] + yo_signif_selec_norm
            if obs_1stderiv[im] > 0:
                n_null1stder_obsm += 1
            if n_null1stder_obsm > obs_n_null1stder_limit:
                break

        #fit selected points by a polynomial
        if len(xobs_chi2) > 9:
            xo_fit_selec = xobs_chi2[3:-3] - lambda_center_obs
            yo_fit_selec = yo_signif_selec_norm[3:-3]
        elif len(xobs_chi2) > 3:
            xo_fit_selec = xobs_chi2 - lambda_center_obs
            yo_fit_selec = yo_signif_selec_norm
        else:
            print('Warning:  too few observational points selected wave range might be too low')

        fit_order = min(int((len(xo_fit_selec)-1)/2)*2+1, 9)
        if fit_order >= 3:
            fitpol = np.polyfit(xo_fit_selec, yo_fit_selec, fit_order)
            width_line_obs = xo_fit_selec.max() - xo_fit_selec.min()
            xo_fit_res = np.arange(xo_fit_selec.min(), xo_fit_selec.max(), 0.001)
            yo_fit_res = np.polyval(fitpol, xo_fit_res)
            xo_fit_res += lambda_center_obs
            yo_fit_min = yo_fit_res.min()
            xo_fit_obs = xo_fit_res[np.argmin(yo_fit_res)]
        else:
            xo_fit_res = -9.99
            yo_fit_res = -9.99
            yo_fit_min = -9.99
            xo_fit_obs = -9.99


        # Could we write a function to avoid a bit of repetition with the above

        # Calculate where the line center is 
        lambda_index = bisect.bisect_left(avg_synth.wavelength, line)
        lambda_index_syn = lambda_index
        ip = lambda_index
        im = lambda_index - 1

        while ((avg_synth.wavelength[ip] < fit_range[1]) & (avg_synth.wavelength[im] > fit_range[0]) & 
               ((avg_synth.variance[lambda_index_syn] <= avg_synth.variance[im]) | 
                (avg_synth.variance[lambda_index_syn] <= avg_synth.variance[ip]) )):

            if (avg_synth.variance[lambda_index_syn] < avg_synth.variance[im]):
                lambda_index_syn = im
            elif (avg_synth.variance[lambda_index_syn] < avg_synth.variance[ip]):
                lambda_index_syn = ip

            ip += 1
            im -= 1

        lambda_center_syn = avg_synth.wavelength[lambda_index_syn]

        # Calculate where the avg synth minimum is

        lambda_index_syn = lambda_index
        ip = lambda_index
        im = lambda_index - 1
        
        while ((avg_synth.wavelength[ip] < fit_range[1]) & (avg_synth.wavelength[im] > fit_range[0]) & 
               ((avg_synth.flux[lambda_index_syn] > avg_synth.flux[im]) | 
                (avg_synth.flux[lambda_index_syn] >= avg_synth.flux[ip]) )):


            if (avg_synth.flux[lambda_index_syn] > avg_synth.flux[ip]):
                lambda_index_syn = ip
            elif (avg_synth.flux[lambda_index_syn] > avg_synth.flux[im]):
                lambda_index_syn = im


            ip += 1
            im -= 1
            
        lambda_min_syn = avg_synth.wavelength[lambda_index_syn]
        
        # select significant points in the synth spectrum
        ip = lambda_index_syn
        xsyn_chi2 = []
        ysyn_signif_selec = [[] for synth in synths]
        fluxavg_signif_selec = []
        n_null1stder_synp = 0
        for ip in range(lambda_index_syn, len(avg_synth.wavelength)):
            if ((n_null1stder_synp >= syn_n_null1stder_limit) | 
                (avg_synth.variance[ip] <= 0.001) | 
                (avg_synth.wavelength[ip] > fit_range[1])):
                    break
            for i, synth in enumerate(synths):
               ysyn_signif_selec[i] += [synth.flux[ip]]
            fluxavg_signif_selec += [avg_synth.flux[ip]]
            xsyn_chi2 += [avg_synth.wavelength[ip]]
            if ((synth_1stderiv[ip] <= 0) & (synth_1stderiv[ip-1] > 0) & 
                (avg_synth.variance[ip] <= fluxdiff_limit)):
                    n_null1stder_synp += 1

        n_null1stder_synm = 0
        for im in range(lambda_index_syn-1, -1, -1):
            if ((n_null1stder_synm >= syn_n_null1stder_limit) | 
                (avg_synth.variance[im] <= 0.001) | 
                (avg_synth.wavelength[im] < fit_range[0])):
                    break
            for i, synth in enumerate(synths):
               ysyn_signif_selec[i] = [synth.flux[im]] + ysyn_signif_selec[i]
            fluxavg_signif_selec = [avg_synth.flux[im]] + fluxavg_signif_selec 
            xsyn_chi2 = [avg_synth.wavelength[im]] + xsyn_chi2
            if ((synth_1stderiv[im] >= 0) & (synth_1stderiv[im+1] < 0) & 
                (avg_synth.variance[im] <= fluxdiff_limit)):
                    n_null1stder_synm += 1

        xsyn_chi2 = np.array(xsyn_chi2)
        fluxavg_signif_selec = np.array(fluxavg_signif_selec)

        #fit selected points by a polynomial
        if len(xsyn_chi2) > 0:
            width_line_syn_ok = ((xsyn_chi2 >= lambda_min_syn - width_line_obs/2.) & 
                                 (xsyn_chi2 <= lambda_min_syn + width_line_obs/2.))
            x_fit_selec = (xsyn_chi2[width_line_syn_ok] - lambda_min_syn)
            y_fit_selec = fluxavg_signif_selec[width_line_syn_ok]
        else:
            width_line_syn_ok = ((avg_synth.wavelength >= lambda_min_syn - width_line_obs/2.) & 
                                 (avg_synth.wavelength <= lambda_min_syn + width_line_obs/2.))
            x_fit_selec = avg_synth.wavelength[width_line_syn_ok] - lambda_min_syn
            y_fit_selec = avg_synth.flux[width_line_syn_ok]


        fitpol = np.polyfit(x_fit_selec, y_fit_selec, 9)
        xs_fit_res = np.arange(x_fit_selec.min(), x_fit_selec.max(), 0.001)
        ys_fit_res = np.polyval(fitpol, xs_fit_res)
        xs_fit_res += lambda_min_syn
        ys_fit_min = ys_fit_res.min()
        xs_fit_syn = xs_fit_res[np.argmin(ys_fit_res)]

        self.xo_fit_res = xo_fit_res
        self.yo_fit_res = yo_fit_res
        self.xs_fit_res = xs_fit_res
        self.ys_fit_res = ys_fit_res
        

        # slightly shift obs spectrum  NOT ACTUALLY SHIFTING AT THE MOMENT
        if xo_fit_obs > 0.:
            shift = xo_fit_obs - xs_fit_syn
        else:
            shift = 0.
            
        resolution = (obs_spec.wavelength.max() - obs_spec.wavelength.min())/len(obs_spec.wavelength)
        large_shift = abs(shift)/2. > resolution

        if len(xsyn_chi2) == 0:
            minxsyn_chi2 = -9.99
            maxxsyn_chi2 = 99999999999
        else:
            minxsyn_chi2 = xsyn_chi2.min()
            maxxsyn_chi2 = xsyn_chi2.max()

        lambda_shift = (lambda_center_obs >= maxxsyn_chi2) | (lambda_center_obs <= minxsyn_chi2)

        xobs_chi2 = np.array(xobs_chi2)

        if len(xsyn_chi2) > 0:
            x_chi2 = xobs_chi2[(xobs_chi2 >= minxsyn_chi2) & (xobs_chi2 <= maxxsyn_chi2)]
            # Taken from BACCHUS, though this seems to always just give you the syn selec range
            if ((len(x_chi2) <= 5) | 
                (len(xsyn_chi2)/len(avg_synth.wavelength)*len(obs_spec.wavelength)/len(xobs_chi2)>=1.5) |
                (len(xsyn_chi2)/len(avg_synth.wavelength)*len(obs_spec.wavelength)/len(xobs_chi2)<=1.5) |
                lambda_shift
               ):
                    fit_range = [minxsyn_chi2,maxxsyn_chi2]
            else:
                fit_range = [x_chi2.min(),x_chi2.max()]
        else:
            fit_range = [xobs_chi2.min(),xobs_chi2.max()]

        # other things that we might want to return could be
        # large_shift, lambda_shift, lambda_min_syn, lambda_center_syn, lambda_center_obs, 
        # xo_fit_obs, xs_fit_syn

        return fit_range, shift, yo_fit_res.min()
    
    
    def _shift_and_reinterpolate(self, obs_spec, synths, line, shift):
        
        c = 2.99792458e5
        vel_offset = -1.*c*(shift/line)
        
        obs_spec.wavelength = obs_spec.vel_shift(vel_offset)
        for synth in synths:
            synth.interpolate_spectrum(obs_spec.wavelength)
            
        return obs_spec, synths


    def _fit_chi2(self, obs_spec, synths, element, line,
                 synth_abund, fit_range):
        '''
        TODO:  Fill out the docstring

        assuming the obs_spec is a Spectrum object and that it is RV corrected

        Need obs_spec to have matching air/vac wavelengths to the used line
            lists

        TODO:  Need to define a fit range that is smaller than the synth range, so
            that it just fits the line/feature of interest

        TODO:  Need to do renormalization
        '''
        
        abundres = 0.0001

        chi2_arr = []

        if fit_range is None:
            fit_range = np.ones(len(obs_spec.wavelength)).astype(bool)
        else:
            fit_range = (obs_spec.wavelength > fit_range[0]) & (obs_spec.wavelength < fit_range[1])


        for spectrum in synths:
            chi2_arr += [np.nanmean(((obs_spec.flux[fit_range]
                                      - spectrum.flux[fit_range])/np.sqrt(obs_spec.variance[fit_range]))**2.)]

        chi2_arr = np.array(chi2_arr)
        
        new_abunds = np.arange(synth_abund.min(),
                               synth_abund.max()+abundres, abundres)
        
        chi2_interp = interp1d(synth_abund, chi2_arr, kind='cubic', fill_value='extrapolate')

        new_chi2 = chi2_interp(new_abunds)

        min_chi2_ind = np.argmin(new_chi2)

        min_chi2 = new_chi2[min_chi2_ind]
        chi2_abund = new_abunds[min_chi2_ind]

        return chi2_abund, min_chi2, chi2_arr


    def _fit_eqw(self, obs_spec, synths, element, line,
                 synth_abund, fit_range, snr):
        '''
        TODO:  Fill out the docstring

        assuming the obs_spec is a Spectrum object and that it is RV corrected

        Need obs_spec to have matching air/vac wavelengths to the used line
            lists

        TODO:  Need to define a fit range that is smaller than the synth range, so
            that it just fits the line/feature of interest

        TODO:  Need to do renormalization
        '''
        
        abundres = 0.0001

        eqw_synth_arr = []

        if fit_range is None:
            fit_range = np.ones(len(obs_spec.wavelength)).astype(bool)
        else:
            fit_range = (obs_spec.wavelength > fit_range[0]) & (obs_spec.wavelength < fit_range[1])

        # integrated eqw
        for spectrum in synths:
             integral_arr = (1 - spectrum.flux[fit_range][:-1]) * \
                (spectrum.wavelength[fit_range][1:] - spectrum.wavelength[fit_range][:-1])
             eqw_synth_arr += [np.nansum(integral_arr)]


        integral_arr = (1 - obs_spec.flux[fit_range][:-1]) * \
                (obs_spec.wavelength[fit_range][1:] - obs_spec.wavelength[fit_range][:-1])
        eqw_obs = np.nansum(integral_arr)

        if len(obs_spec.variance) > 1:
            integral_var_arr = (obs_spec.variance[fit_range][:-1]) * \
                ((obs_spec.wavelength[fit_range][1:] - obs_spec.wavelength[fit_range][0:-1]))**2.
            eqw_obs_err = np.sqrt(np.nansum(integral_var_arr))

        else:
            integral_var_arr = (1/snr**2.) * \
                ((obs_spec.wavelength[fit_range][1:] - obs_spec.wavelength[fit_range][0:-1]))**2.
            eqw_obs_err = np.sqrt(np.nansum(integral_var_arr))


        eqw_synth_arr = np.array(eqw_synth_arr)*1000
        eqw_obs *= 1000.
        eqw_obs_err *= 1000.
        
        eqw_interp = interp1d(eqw_synth_arr, synth_abund, kind='cubic', fill_value='extrapolate')

        eqw_abund = eqw_interp(eqw_obs)
        eqw_abund_err = (eqw_interp(eqw_obs + eqw_obs_err) - eqw_interp(eqw_obs - eqw_obs_err))/2.
        eqw_abund_5sig_lim = eqw_interp(5*eqw_obs_err)

        eqw_lim_flag = int(eqw_abund > eqw_abund_5sig_lim)

        return eqw_obs, eqw_obs_err, eqw_abund, eqw_abund_err, eqw_abund_5sig_lim, eqw_lim_flag



class Eqwometer:
    def __init__(self, turbosynth_obj, input_obs_spec, element, line,
             conv_profiles, conv_broad,
             delta_abunds=[-0.6, -0.3, 0.0, 0.3, 0.6, ],
             synth_range=None, delta_lambda=0.01):
        
        self.element = element
        self.line = line
        self.delta_abunds = delta_abunds
        self.delta_lambda = delta_lambda
        self.conv_profiles = conv_profiles
        self.conv_broad = conv_broad
        self.turbosynth_obj = turbosynth_obj

        if synth_range is None:
            synth_range = line/2000. + 1

        self.wave_range = [line - synth_range - 1, line + synth_range + 1]
        self.fit_range = [line - synth_range, line + synth_range]

        selection = (input_obs_spec.wavelength > self.wave_range[0]) & \
                    (input_obs_spec.wavelength < self.wave_range[1])
        obs_waves = input_obs_spec.wavelength[selection]
        obs_flux = input_obs_spec.flux[selection]
        if hasattr(input_obs_spec, 'variance'):
            obs_var = input_obs_spec.variance[selection]
        else:
            obs_var = 0.1**2.

        self.obs_spec = Spectrum(obs_waves, obs_flux, variance=obs_var)

    def fit_ts(self):
        self.run_eqwidt()

        
    def fit(self, eqw):
        
        self.synthesize()
        self.measure_abund(eqw)

    def measure_synth_eqw(self):
        self.synthesize()
        eqws = self._measure_eqw(
            self.obs_spec, self.synths, self.element, self.line, 
            self.synth_abunds, self.fit_range)

        return eqws

    def run_eqwidt(self, delta_abunds=None, conv_profiles=None, conv_broad=None):
        
        if delta_abunds is not None:
            self.delta_abunds = delta_abunds
        if conv_profiles is not None:
            self.conv_profiles = conv_profiles
        if conv_broad is not None:
            self.conv_broad = conv_broad
    
        abund_fname = self._calculate_eqwidt_abund(self.turbosynth_obj, 
                                                      self.obs_spec, self.element, 
                                                      self.line, self.conv_profiles, 
                                                      self.conv_broad, self.delta_abunds,
                                                      self.wave_range, self.delta_lambda)

        self.abund_fname = abund_fname


    def synthesize(self, delta_abunds=None, conv_profiles=None, conv_broad=None):
        
        if delta_abunds is not None:
            self.delta_abunds = delta_abunds
        if conv_profiles is not None:
            self.conv_profiles = conv_profiles
        if conv_broad is not None:
            self.conv_broad = conv_broad
    
        synth_abunds, synths = self._calculate_synths(self.turbosynth_obj, 
                                                      self.obs_spec, self.element, 
                                                      self.line, self.conv_profiles, 
                                                      self.conv_broad, self.delta_abunds,
                                                      self.wave_range, self.delta_lambda)

        if len(self.delta_abunds) > 1:
            synths, avg_synth = self._make_avg_synth(synths)
            self.avg_synth = avg_synth
            
        self.synth_abunds = synth_abunds
        self.synths = synths

        
    def measure_abund(self, eqw):
        
        eqw_abund = self._fit_eqw(
            eqw, self.obs_spec, self.synths, self.element, self.line, 
            self.synth_abunds, self.fit_range)

        self.eqw = eqw
        self.fit_abund_eqw = eqw_abund


    def _calculate_synths(self, turbosynth_obj, obs_spec, element, line,
                         conv_profiles, conv_broad, abunds,
                         wave_range, delta_lambda):
        '''
        TODO:  Fill out the docstring

        assuming the obs_spec is a Spectrum object and that it is RV corrected

        Need obs_spec to have matching air/vac wavelengths to the used line
            lists

        TODO:  Need to define a fit range that is smaller than the synth range, so
            that it just fits the line/feature of interest

        TODO:  Need to do renormalization
        '''

        elemabund = turbosynth_obj.get_abund(element)

        synths = []
        chi2_arr = []

        synth_abund = np.array([elemabund+delta_xfe for delta_xfe in abunds])


        for delta_xfe in abunds:
            turbosynth_obj.set_abund(element, elemabund+delta_xfe)
            # TODO: Might want to make an option to run a single atmosphere
            # with the default abundances which could speed things up and keep
            # the code from breaking on random abundances
            turbosynth_obj.make_atmosphere(star='tmp')

            turbosynth_obj.synth(wave_range, synth_fname='tmp.spec',
                       delta_lambda=delta_lambda)
            os.remove(turbosynth_obj.opac_filename)
            os.remove(turbosynth_obj.opac_filename+'.mod')

            for i, (profile, kwargs) in enumerate(zip(conv_profiles,
                                                      conv_broad)):
                if i == 0:
                    fname = turbosynth_obj.synth_fname
                else:
                    fname = convol
                convol = turbosynth_obj.convol(synth_fname=fname, profile=profile,
                                     **kwargs, convol_fname='tmp.spec.convol')

            spectrum = SpectrumTS(f'{turbosynth_obj.path}/tmp.spec.convol')
            spectrum.interpolate_spectrum(obs_spec.wavelength)

            synths += [spectrum]

            turbosynth_obj.set_abund(element, elemabund)

        return synth_abund, synths


    def _calculate_eqwidt_abund(self, turbosynth_obj, obs_spec, element, line,
                         conv_profiles, conv_broad, abunds,
                         wave_range, delta_lambda):
        '''
        TODO:  Fill out the docstring

        assuming the obs_spec is a Spectrum object and that it is RV corrected

        Need obs_spec to have matching air/vac wavelengths to the used line
            lists

        TODO:  Need to define a fit range that is smaller than the synth range, so
            that it just fits the line/feature of interest

        TODO:  Need to do renormalization
        '''


        turbosynth_obj.make_atmosphere(star='tmp')

        turbosynth_obj.eqwidth(wave_range, abund_fname='tmp.findabu',
                   delta_lambda=delta_lambda)
        os.remove(turbosynth_obj.opac_filename)
        os.remove(turbosynth_obj.opac_filename+'.mod')

        fname = turbosynth_obj.abund_fname

        return fname

    def _make_avg_synth(self, synth_specs):

        # Rescale syntheses and calculate avg synthetic spectrum
        synth_flux_arr = np.zeros((len(synth_specs[0].wavelength),len(synth_specs)))
        for i in range(len(synth_specs)):
            synth_flux_arr[:,i] = synth_specs[i].flux
        stdev_synth = np.std(synth_flux_arr, axis=1)

        if stdev_synth.min() > 0.005:
            rescale_value = synth_specs[0].flux.max()
            for i, synth in enumerate(synth_specs):
                synth.flux /= rescale_value
                synth_flux_arr[:,i] = synth.flux

        avg_flux = np.mean(synth_flux_arr, axis=1)
        stdev_synth = np.std(synth_flux_arr, axis=1)
        spread_synth = np.std(synth_flux_arr - np.tile(synth_flux_arr[:,0], (len(synth_specs),1)).T, axis=1)
        avg_synth = Spectrum(synth_specs[0].wavelength, avg_flux, variance=spread_synth)

        return synth_specs, avg_synth


    
    def _fit_eqw(self, eqw_obs, obs_spec, synths, element, line,
                 synth_abund, fit_range):
        '''
        TODO:  Fill out the docstring

        assuming the obs_spec is a Spectrum object and that it is RV corrected

        Need obs_spec to have matching air/vac wavelengths to the used line
            lists

        TODO:  Need to define a fit range that is smaller than the synth range, so
            that it just fits the line/feature of interest

        TODO:  Need to do renormalization
        '''
        
        abundres = 0.0001

        eqw_synth_arr = []

        if fit_range is None:
            fit_range = np.ones(len(obs_spec.wavelength)).astype(bool)
        else:
            fit_range = (obs_spec.wavelength > fit_range[0]) & (obs_spec.wavelength < fit_range[1])

        # integrated eqw
        for spectrum in synths:
             integral_arr = (1 - spectrum.flux[fit_range][:-1]) * \
                (spectrum.wavelength[fit_range][1:] - spectrum.wavelength[fit_range][:-1])
             eqw_synth_arr += [np.nansum(integral_arr)]


        eqw_synth_arr = np.array(eqw_synth_arr)*1000
        
        eqw_interp = interp1d(eqw_synth_arr, synth_abund, kind='cubic', fill_value='extrapolate')

        eqw_abund = eqw_interp(eqw_obs)

        return eqw_abund

    def _measure_eqw(self, obs_spec, synths, element, line,
                 synth_abund, fit_range):

        eqw_synth_arr = []

        if fit_range is None:
            fit_range = np.ones(len(obs_spec.wavelength)).astype(bool)
        else:
            fit_range = (obs_spec.wavelength > fit_range[0]) & (obs_spec.wavelength < fit_range[1])

        # integrated eqw
        for spectrum in synths:
             integral_arr = (1 - spectrum.flux[fit_range][:-1]) * \
                (spectrum.wavelength[fit_range][1:] - spectrum.wavelength[fit_range][:-1])
             eqw_synth_arr += [np.nansum(integral_arr)]


        eqw_synth_arr = np.array(eqw_synth_arr)*1000

        return eqw_synth_arr
