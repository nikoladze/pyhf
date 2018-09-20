import copy
import logging
log = logging.getLogger(__name__)

from . import get_backend
from . import exceptions
from . import modifiers
from . import utils


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

def finalize_stats(modifier):
    tensorlib, _ = get_backend()
    inquad = tensorlib.sqrt(tensorlib.sum(tensorlib.power(tensorlib.astensor(modifier.uncertainties),2), axis=0))
    totals = tensorlib.sum(modifier.nominal_counts,axis=0)
    return tensorlib.divide(inquad,totals)

class Model(object):
    def __init__(self, spec, **config_kwargs):
        self.spec = copy.deepcopy(spec) #may get modified by config
        self.schema = config_kwargs.get('schema', utils.get_default_schema())
        # run jsonschema validation of input specification against the (provided) schema
        log.info("Validating spec against schema: {0:s}".format(self.schema))
        utils.validate(self.spec, self.schema)
        # build up our representation of the specification
        self.config = _ModelConfig.from_spec(self.spec,**config_kwargs)
        self.cube, self.hm = self._make_cube()
        from .modifiers.combined import CombinedInterpolator,TrivialCombined

        self.finalized_stats = {k:finalize_stats(self.config.modifier(k)) for k,v in self.config.par_map.items() if 'staterror' in k}
        self.allmods = self._make_mod_index()
        self.combined_mods = {k:CombinedInterpolator(self,k) for k in ['normsys','histosys']}
        self.trivial_combined = {k:TrivialCombined(self,k) for k in ['normfactor','shapesys','shapefactor','staterror']}

    def _make_cube(self):
        import  numpy as np
        nchan   = len(self.spec['channels'])
        maxsamp = max(len(c['samples']) for c in self.spec['channels'])
        maxbins = max(len(s['data']) for c in self.spec['channels'] for s in c['samples'])
        cube = np.ones((nchan,maxsamp,maxbins))*np.nan
        histoid = 0
        histomap = {}
        for i,c in enumerate(self.spec['channels']):
            for j,s in enumerate(c['samples']):
                cube[i,j,:] = s['data']
                histomap.setdefault(c['name'],{})[s['name']] = {'id': histoid, 'index': tuple([i,j])}
                histoid += 1
        return cube,histomap

    def _make_mod_index(self):
        factor_mods = ['normfactor','normsys','shapesys','shapefactor','staterror']
        delta_mods  = ['histosys']
        allmods = {}
        for channel in self.spec['channels']:
            for sample in channel['samples']:
                for m in sample['modifiers']:
                    mname, mtype = m['name'], m['type']
                    if mtype in factor_mods:
                        opcode_id = 0
                    elif mtype in delta_mods:
                        opcode_id = 1
                    allmods.setdefault(mname, [mtype,opcode_id])
        for mod_id,mname in enumerate(allmods.keys()):
            allmods[mname].append(mod_id)
        return allmods
        
    def _mtype_results(self,mtype,pars):
        if mtype in self.combined_mods.keys():
            return self.combined_mods[mtype].apply(pars)
        if mtype in ['normfactor','shapesys','shapefactor','staterror']:
            return self.trivial_combined[mtype].apply(pars)

    def _make_result_index(self,pars):
        factor_mods = ['normsys','normfactor','shapesys','shapefactor','staterror']
        delta_mods  = ['histosys']
        all_results = {}
        for mtype in factor_mods + delta_mods:
            all_results[mtype] = self._mtype_results(mtype,pars)
        result_index_list = []
        for mname,(mtype,_,_) in self.allmods.items():
            result_index_list += all_results[mtype].get(mname,[])
        return len(self.allmods),result_index_list


    def _all_modifications(self, pars):
        ntotalmods, result_index = self._make_result_index(pars)

        import numpy as np

        op_fields =  np.stack([
            np.ones((ntotalmods,) + self.cube.shape),
            np.zeros((ntotalmods,) + self.cube.shape)
        ])

        # indices, vals = list(zip(*result_index))
        # for total,r in zip(indices,vals):
        # possible vectorization with np.put?
        for res in result_index:
            total, r = res
            op_fields[total] = r
        factor_field = np.product(op_fields[0],axis=0)
        delta_field  = np.sum(op_fields[1],axis=0)
        return factor_field,delta_field

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

    def expected_actualdata(self, pars):
        import numpy as np
        tensorlib, _ = get_backend()
        pars = tensorlib.astensor(pars)
        data = []

        factor_field, delta_field = self._all_modifications(pars)
        combined = factor_field * (delta_field + self.cube)
        expected = [np.sum(combined[i][~np.isnan(combined[i])]) for i,c in enumerate(self.spec['channels'])]
        return expected

    def expected_data(self, pars, include_auxdata=True):
        tensorlib, _ = get_backend()
        pars = tensorlib.astensor(pars)
        expected_actual = self.expected_actualdata(pars)

        if not include_auxdata:
            return expected_actual
        expected_constraints = self.expected_auxdata(pars)
        tocat = [expected_actual] if expected_constraints is None else [expected_actual,expected_constraints]
        return tensorlib.concatenate(tocat)
    
    def __calculate_constraint(self,bytype):
        tensorlib, _ = get_backend()
        newsummands = None
        for k,c in bytype.items():
            c = tensorlib.astensor(c)
            #warning, call signature depends on pdf_type (2 for pois, 3 for normal)
            pdfval = getattr(tensorlib,k)(c[:,0],c[:,1],c[:,2])
            constraint_term = tensorlib.log(pdfval)
            newsummands = constraint_term if newsummands is None else tensorlib.concatenate([newsummands,constraint_term])
        return tensorlib.sum(newsummands) if newsummands is not None else 0

    def constraint_logpdf(self, auxdata, pars):
        tensorlib, _ = get_backend()
        start_index = 0
        bytype = {}
        for cname in self.config.auxdata_order:
            modifier, modslice = self.config.modifier(cname), \
                self.config.par_slice(cname)
            modalphas = modifier.alphas(pars[modslice])
            end_index = start_index + int(modalphas.shape[0])
            thisauxdata = auxdata[start_index:end_index]
            start_index = end_index
            if modifier.pdf_type=='normal':
                if modifier.__class__.__name__ in ['histosys','normsys']:
                    kwargs = {'sigma': tensorlib.astensor([1])}
                elif modifier.__class__.__name__ in ['staterror']:
                    kwargs = {'sigma': self.finalized_stats[cname]}
            else:
                kwargs = {}
            callargs = [thisauxdata,modalphas] + [kwargs['sigma'] if kwargs else []]
            bytype.setdefault(modifier.pdf_type,[]).append(callargs)
        return self.__calculate_constraint(bytype)

    def logpdf(self, pars, data):
        tensorlib, _ = get_backend()
        pars, data = tensorlib.astensor(pars), tensorlib.astensor(data)
        cut = int(data.shape[0]) - len(self.config.auxdata)
        actual_data, aux_data = data[:cut], data[cut:]
        lambdas_data = self.expected_actualdata(pars)
        summands = tensorlib.log(tensorlib.poisson(actual_data, lambdas_data))

        result = tensorlib.sum(summands) + self.constraint_logpdf(aux_data, pars)
        return tensorlib.astensor(result) * tensorlib.ones((1)) #ensure (1,) array shape also for numpy

    def pdf(self, pars, data):
        tensorlib, _ = get_backend()
        return tensorlib.exp(self.logpdf(pars, data))
