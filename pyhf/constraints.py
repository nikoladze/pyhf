from . import get_backend

class gaussian_constraint_combined(object):
    def __init__(self,pdf):
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

        par_indices = list(range(len(pdf.config.suggested_init())))
        data_indices = list(range(len(pdf.config.auxdata)))
        for cname in pdf.config.auxdata_order:
            modifier = pdf.config.modifier(cname)
            modslice  = pdf.config.par_slice(cname)

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
        self.prepped = {
            'normal': (normal_data,normal_sigmas,normal_mean_idc),
            'poisson': (poisson_data,poisson_rate_idc)
        }