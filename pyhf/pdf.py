import copy
import logging
log = logging.getLogger(__name__)

from . import get_backend
from . import exceptions
from . import modifiers
from . import utils
from .interpolate import _hfinterp_code1,_hfinterp_code0


class _ModelConfig(object):
    @classmethod
    def from_spec(cls,spec,poiname = 'mu', qualify_names = False):
        channels = []
        samples = []
        modifiers = []
        # hacky, need to keep track in which order we added the constraints
        # so that we can generate correctly-ordered data
        instance = cls()
        for channel in spec['channels']:
            channels.append(channel['name'])
            for sample in channel['samples']:
                samples.append(sample['name'])
                # we need to bookkeep a list of modifiers by type so that we
                # can loop over them on a type-by-type basis
                # types like histosys, normsys, etc...
                sample['modifiers_by_type'] = {}
                for modifier_def in sample['modifiers']:
                    if qualify_names:
                        fullname = '{}/{}'.format(modifier_def['type'],modifier_def['name'])
                        if modifier_def['name'] == poiname:
                            poiname = fullname
                        modifier_def['name'] = fullname
                    modifier = instance.add_or_get_modifier(channel, sample, modifier_def)
                    modifier.add_sample(channel, sample, modifier_def)
                    modifiers.append(modifier_def['name'])
                    sample['modifiers_by_type'].setdefault(modifier_def['type'],[]).append(modifier_def['name'])
        instance.channels = list(set(channels))
        instance.samples = list(set(samples))
        instance.modifiers = list(set(modifiers))
        instance.set_poi(poiname)
        return instance

    def __init__(self):
        # set up all other bookkeeping variables
        self.poi_index = None
        self.par_map = {}
        self.par_order = []
        self.auxdata = []
        self.auxdata_order = []
        self.next_index = 0

    def suggested_init(self):
        init = []
        for name in self.par_order:
            init = init + self.par_map[name]['modifier'].suggested_init
        return init

    def suggested_bounds(self):
        bounds = []
        for name in self.par_order:
            bounds = bounds + self.par_map[name]['modifier'].suggested_bounds
        return bounds

    def par_slice(self, name):
        return self.par_map[name]['slice']

    def modifier(self, name):
        return self.par_map[name]['modifier']

    def set_poi(self,name):
        if name not in self.modifiers:
            raise exceptions.InvalidModel("The paramter of interest '{0:s}' cannot be fit as it is not declared in the model specification.".format(name))
        s = self.par_slice(name)
        assert s.stop-s.start == 1
        self.poi_index = s.start

    def add_or_get_modifier(self, channel, sample, modifier_def):
        """
        Add a new modifier if it does not exist and return it
        or get the existing modifier and return it

        Args:
            channel: current channel object (e.g. from spec)
            sample: current sample object (e.g. from spec)
            modifier_def: current modifier definitions (e.g. from spec)

        Returns:
            modifier object

        """
        # get modifier class associated with modifier type
        try:
            modifier_cls = modifiers.registry[modifier_def['type']]
        except KeyError:
            log.exception('Modifier type not implemented yet (processing {0:s}). Current modifier types: {1}'.format(modifier_def['type'], modifiers.registry.keys()))
            raise exceptions.InvalidModifier()

        # if modifier is shared, check if it already exists and use it
        if modifier_cls.is_shared and modifier_def['name'] in self.par_map:
            log.info('using existing shared, {0:s}constrained modifier (name={1:s}, type={2:s})'.format('' if modifier_cls.is_constrained else 'un', modifier_def['name'], modifier_cls.__name__))
            modifier = self.par_map[modifier_def['name']]['modifier']
            if not type(modifier).__name__ == modifier_def['type']:
                raise exceptions.InvalidNameReuse('existing modifier is found, but it is of wrong type {} (instead of {}). Use unique modifier names or use qualify_names=True when constructing the pdf.'.format(type(modifier).__name__, modifier_def['type']))
            return modifier

        # did not return, so create new modifier and return it
        modifier = modifier_cls(sample['data'], modifier_def['data'])
        npars = modifier.n_parameters

        log.info('adding modifier %s (%s new nuisance parameters)', modifier_def['name'], npars)
        sl = slice(self.next_index, self.next_index + npars)
        self.next_index = self.next_index + npars
        self.par_order.append(modifier_def['name'])
        self.par_map[modifier_def['name']] = {
            'slice': sl,
            'modifier': modifier
        }
        if modifier.is_constrained:
            self.auxdata += self.modifier(modifier_def['name']).auxdata
            self.auxdata_order.append(modifier_def['name'])
        return modifier

class Model(object):
    def __init__(self, spec, **config_kwargs):
        self.spec = copy.deepcopy(spec) #may get modified by config
        self.schema = config_kwargs.get('schema', utils.get_default_schema())
        # run jsonschema validation of input specification against the (provided) schema
        log.info("Validating spec against schema: {0:s}".format(self.schema))
        utils.validate(self.spec, self.schema)
        # build up our representation of the specification
        self.config = _ModelConfig.from_spec(self.spec,**config_kwargs)
        
        _allmods = []
        _allsamples = []
        _allchannels = []
        _allmods = []
        channel_nbins = {}

        for c in self.spec['channels']:
            _allchannels.append(c['name'])
            for s in c['samples']:
                channel_nbins[c['name']] = len(s['data'])
                _allsamples.append(s['name'])
                for mod in s['modifiers']:
                    _allmods.append((mod['name'],mod['type']))
        _allmods = list(set(_allmods))
        _allsamples = list(set(_allsamples))
        _allchannels = list(set(_allchannels))
        self.do_samples  = list(sorted(_allsamples[:]))
        self.do_channels = list(sorted(_allchannels[:]))
        self.do_mods = list(sorted(_allmods[:]))
        self.channel_nbins = channel_nbins
        self._make_mega()
        self._prep_mega()
        self.prepped_constraints = self.prep_constraints()


    def _make_mega(self):
        helper = {}
        for c in self.spec['channels']:
            for s in c['samples']:
                helper.setdefault(c['name'],{})[s['name']] = (c,s)

        mega_mods = {}
        import copy
        for m,mtype in self.do_mods:
            for s in self.do_samples:
                modspec = {'type': mtype, 'name': m}
                if mtype == 'histosys':
                    modspec.setdefault('data',{})['hi_data'] = []
                    modspec.setdefault('data',{})['lo_data'] = []
                    modspec.setdefault('data',{})['mask'] = []
                elif mtype == 'normsys':
                    modspec.setdefault('data',{})['hi'] = []
                    modspec.setdefault('data',{})['lo'] = []
                    modspec.setdefault('data',{})['mask'] = []
                elif mtype == 'normfactor':
                    modspec.setdefault('data',{})['mask'] = []
                elif mtype == 'shapefactor':
                    modspec.setdefault('data',{})['mask'] = []
                elif mtype == 'shapesys':
                    modspec.setdefault('data',{})['uncrt'] = []
                    modspec.setdefault('data',{})['mask'] = []
                elif mtype == 'staterror':
                    modspec.setdefault('data',{})['uncrt'] = []
                    modspec.setdefault('data',{})['mask']  = []
                mega_mods.setdefault(s,{})[m] = copy.deepcopy(modspec)
                
        mega_samples = {}
        for s in self.do_samples:
            mega_nom = []
            for c in self.do_channels:
                defined_samp = helper.get(c,{}).get(s)
                defined_samp = None if not defined_samp else defined_samp[1]
                nom = defined_samp['data'] if defined_samp else [0.0]*self.channel_nbins[c]
                mega_nom += nom
                defined_mods = {x['name']:x for x in defined_samp['modifiers']} if defined_samp else {}
                for m,mtype in self.do_mods:
                    thismod = defined_mods.get(m)
                    if mtype == 'histosys':
                        lo_data = thismod['data']['lo_data'] if thismod else nom
                        hi_data = thismod['data']['hi_data'] if thismod else nom
                        maskval = True if thismod else False
                        mega_mods[s][m]['data']['lo_data'] += lo_data
                        mega_mods[s][m]['data']['hi_data'] += hi_data
                        mega_mods[s][m]['data']['mask']    += [maskval]*len(nom) #broadcasting
                        pass
                    elif mtype == 'normsys':
                        maskval = True if thismod else False
                        lo_factor = thismod['data']['lo'] if thismod else 1.0
                        hi_factor = thismod['data']['hi'] if thismod else 1.0
                        mega_mods[s][m]['data']['lo']   += [lo_factor]*len(nom) #broadcasting
                        mega_mods[s][m]['data']['hi']   += [hi_factor]*len(nom)
                        mega_mods[s][m]['data']['mask'] += [maskval]  *len(nom) #broadcasting
                    elif mtype == 'normfactor':
                        maskval = True if thismod else False
                        mega_mods[s][m]['data']['mask'] += [maskval]*len(nom) #broadcasting
                    elif mtype == 'shapesys':
                        uncrt = thismod['data'] if thismod else [0.0]*len(nom)
                        maskval = [True if thismod else False]*len(nom)
                        mega_mods[s][m]['data']['mask']  += maskval
                        mega_mods[s][m]['data']['uncrt'] += uncrt
                    elif mtype == 'staterror':
                        uncrt = thismod['data'] if thismod else [0.0]*len(nom)
                        maskval = [True if thismod else False]*len(nom)
                        mega_mods[s][m]['data']['mask']  += maskval
                        mega_mods[s][m]['data']['uncrt'] += uncrt
                    elif mtype == 'shapefactor':
                        maskval = True if thismod else False
                        mega_mods[s][m]['data']['mask'] += [maskval]*len(nom) #broadcasting
                    else:
                        raise RuntimeError('not sure how to combine {mtype} into the mega-channel'.format(mtype = mtype))
            sample_dict = {
                'name': 'mega_{}'.format(s),
                'nom': mega_nom,
                'modifiers': list(mega_mods[s].values())
            }
            mega_samples[s] = sample_dict
        self.mega_samples = mega_samples
        self.mega_mods    = mega_mods

    def _prep_mega(self):
        tensorlib,_ = get_backend()


        from modifiers.combined_mods import normsys_combinedmod,histosys_combinedmod,normfac_combinedmod,staterror_combined, shapesys_combined

        normsys_mods = [m for m,mtype in self.do_mods if mtype == 'normsys' ]
        self.normsys_combined = normsys_combinedmod(normsys_mods,self)

        histosys_mods = [m for m,mtype in self.do_mods if mtype == 'histosys' ]
        self.histosys_combined = histosys_combinedmod(histosys_mods,self)

        normfac_mods = [m for m,mtype in self.do_mods if mtype == 'normfactor']
        self.normfac_combined = normfac_combinedmod(normfac_mods,self)

        staterr_mods = [m for m,mtype in self.do_mods if mtype == 'staterror']
        for m in staterr_mods:
            self.config.modifier(m).finalize()
        self.staterr_combined = staterror_combined(staterr_mods,self)

        shapesys_mods = [m for m,mtype in self.do_mods if mtype == 'shapesys']
        self.shapesys_combined = shapesys_combined(shapesys_mods,self)

        thenom = tensorlib.astensor(
            [self.mega_samples[s]['nom'] for s in self.do_samples]
        )
        self.thenom = tensorlib.reshape(thenom,(
            1,
            len(self.do_samples),
            1,
            sum(list(self.channel_nbins.values()))
            )
         )


    def expected_auxdata(self, pars):
        # probably more correctly this should be the expectation value of the constraint_pdf
        # or for the constraints we are using (single par constraings with mean == mode), we can
        # just return the alphas

        tensorlib, _ = get_backend()
        # order matters! because we generated auxdata in a certain order
        auxdata = None
        for modname in self.config.auxdata_order:
            thisaux = self.config.modifier(modname).expected_data(
                pars[self.config.par_slice(modname)])
            tocat = [thisaux] if auxdata is None else [auxdata, thisaux]
            auxdata = tensorlib.concatenate(tocat)
        return auxdata

    def _modifications(self,pars):
        tensorlib, _ = get_backend()

        pars = tensorlib.astensor(pars)

        results_norm     = self.normsys_combined.apply(pars)
        results_histo    = self.histosys_combined.apply(pars)
        results_staterr  = self.staterr_combined.apply(pars)
        results_shapesys = self.shapesys_combined.apply(pars)
        results_normfac  = self.normfac_combined.apply(pars)

        deltas  = list(filter(lambda x: x is not None,[results_histo]))
        factors = list(filter(lambda x: x is not None,[
                results_norm,
                results_staterr,
                results_shapesys,
                results_normfac
        ]))
        return deltas, factors

    def expected_actualdata(self,pars):
        deltas, factors = self._modifications(pars)
        
        tensorlib, _ = get_backend()
        allsum = tensorlib.concatenate(deltas + [self.thenom])
        
        nom_plus_delta = tensorlib.sum(allsum,axis=0)
        nom_plus_delta = tensorlib.reshape(nom_plus_delta,(1,)+tensorlib.shape(nom_plus_delta))

        allfac = tensorlib.concatenate(factors + [nom_plus_delta])

        newbysample = tensorlib.product(allfac,axis=0)
        newresults = tensorlib.sum(newbysample,axis=0)
        return newresults[0] #only one alphas

    def expected_data(self, pars, include_auxdata=True):
        tensorlib, _ = get_backend()
        pars = tensorlib.astensor(pars)
        expected_actual = self.expected_actualdata(pars)

        if not include_auxdata:
            return expected_actual
        expected_constraints = self.expected_auxdata(pars)
        tocat = [expected_actual] if expected_constraints is None else [expected_actual,expected_constraints]
        return tensorlib.concatenate(tocat)
        
    def prep_constraints(self):
        tensorlib, _ = get_backend()
        # iterate over all constraints order doesn't matter....
        start_index = 0
        summands = None
        
        normal_constraint_data = []
        normal_constraint_mean_indices = []
        normal_constraint_sigmas = []
        
        poisson_constraint_rates = []
        
        poisson_constraint_data = []
        poisson_constraint_rate_indices = []

        par_indices = list(range(len(self.config.suggested_init())))
        data_indices = list(range(len(self.config.auxdata)))
        for cname in self.config.auxdata_order:
            modifier = self.config.modifier(cname)
            modslice  = self.config.par_slice(cname)

            end_index = start_index + modifier.n_parameters
            thisauxdata = data_indices[start_index:end_index]
            start_index = end_index
            if modifier.pdf_type == 'normal':
                normal_constraint_sigmas.append(modifier.sigmas)
                normal_constraint_data.append(thisauxdata)
                normal_constraint_mean_indices.append(par_indices[modslice])
            elif modifier.pdf_type == 'poisson':
                poisson_constraint_data.append(thisauxdata)
                poisson_constraint_rate_indices.append(par_indices[modslice])
            else:
                raise RuntimeError
        
        if normal_constraint_mean_indices:
            normal_mean_idc  = tensorlib.concatenate(map(lambda x: tensorlib.astensor(x,dtype = 'int'),normal_constraint_mean_indices))
            normal_sigmas    = tensorlib.concatenate(map(tensorlib.astensor,normal_constraint_sigmas))
            normal_data      = tensorlib.concatenate(map(lambda x: tensorlib.astensor(x,dtype = 'int'),normal_constraint_data))
        else:
            normal_data, normal_sigmas, normal_mean_idc = None, None, None

        if poisson_constraint_rate_indices:
            poisson_rate_idc  = tensorlib.concatenate(map(lambda x: tensorlib.astensor(x,dype = 'int'), poisson_constraint_rate_indices))
            poisson_data      = tensorlib.concatenate(map(lambda x: tensorlib.astensor(x,dtype = 'int'), poisson_constraint_data))
        else:
            poisson_rate_idc, poisson_data = None, None
        return {
            'normal': (normal_data,normal_sigmas,normal_mean_idc),
            'poisson': (poisson_data,poisson_rate_idc)
        }

    def constraint_logpdf(self, auxdata, pars):
        tensorlib, _ = get_backend()
        
        prepped = self.prepped_constraints

        toconcat = []
        if prepped['normal'][0] is not None:
            normal_data   = tensorlib.gather(auxdata,prepped['normal'][0])
            normal_means  = tensorlib.gather(pars,prepped['normal'][2])
            normal_sigmas = prepped['normal'][1]
            normal = tensorlib.normal_logpdf(normal_data,normal_means,normal_sigmas)
            toconcat.append(normal)

        if prepped['poisson'][0] is not None:
            poisson_data  = tensorlib.gather(auxdata,prepped['poisson'][0])
            poisson_rate  = tensorlib.gather(pars,prepped['poisson'][1])
            poisson = tensorlib.poisson_logpdf(poisson_data,poisson_rate)
            toconcat.append(poisson)
        if not toconcat:
            return 0
        all_constraints = tensorlib.concatenate(toconcat)
        return tensorlib.sum(all_constraints)

    def mainlogpdf(self, maindata, pars):
        tensorlib, _ = get_backend()
        lambdas_data = self.expected_actualdata(pars)
        summands   = tensorlib.poisson_logpdf(maindata, lambdas_data)
        tosum      = tensorlib.boolean_mask(summands,tensorlib.isfinite(summands))
        mainpdf    = tensorlib.sum(tosum)
        return mainpdf

    def logpdf(self, pars, data):
        try:
            tensorlib, _ = get_backend()
            pars, data = tensorlib.astensor(pars), tensorlib.astensor(data)
            cut = tensorlib.shape(data)[0] - len(self.config.auxdata)
            actual_data, aux_data = data[:cut], data[cut:]

            mainpdf    = self.mainlogpdf(actual_data,pars)
            constraint = self.constraint_logpdf(aux_data, pars)
            
            result = mainpdf + constraint
            return tensorlib.astensor(result) * tensorlib.ones((1)) #ensure (1,) array shape also for numpy
        except:
            log.error('eval failed for data {} pars: {}'.format(
                tensorlib.tolist(data),
                tensorlib.tolist(pars)
            ))
            raise

    def pdf(self, pars, data):
        tensorlib, _ = get_backend()
        return tensorlib.exp(self.logpdf(pars, data))
