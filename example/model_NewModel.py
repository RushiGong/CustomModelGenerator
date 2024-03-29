
import copy
import warnings
from symengine import exp, log, Abs, Add, And, Float, Mul, Piecewise, Pow, S, sin, StrictGreaterThan, Symbol, zoo, oo
from tinydb import where
import pycalphad.variables as v
from pycalphad.core.errors import DofError
from pycalphad.core.constants import MIN_SITE_FRACTION
from pycalphad.core.utils import unpack_components, get_pure_elements, wrap_symbol
import numpy as np
from pycalphad import Model
from collections import OrderedDict

class NewModel(Model):
    contributions = [
	("ref", "reference_energy"),
	("idmix", "ideal_mixing_energy"),
	("xsmix", "excess_mixing_energy")
	]
    def v_i(self, dbe, i, :,  , v, ., S, p, e, c, i, e, s):
    	param_query=(
			(where("phase_name") == self.phase_name) & \
			(where("parameter_type") == "NewModelV") & \
			(where("constituent_array").test(self._array_validity))
		)
		params = dbe._parameters.search(param_query)
    	#volume parameter
    	return v_i

    def z(self, dbe):
    	param_query=(
			(where("phase_name") == self.phase_name) & \
			(where("parameter_type") == "NewModelZ") & \
			(where("constituent_array").test(self._array_validity))
		)
		params = dbe._parameters.search(param_query)
    	#coordination number
    	return z

    def b_ij(self, dbe, i: v.Species, j: v.Species):
    	param_query=(
			(where("phase_name") == self.phase_name) & \
			(where("parameter_type") == "NewModelB") & \
			(where("constituent_array").test(self._array_validity))
		)
		params = dbe._parameters.search(param_query)
    	return b_ij

    
    def reference_energy(self, dbe):
        
        #Returns the weighted average of the endmember energies
        #in symbolic form.
        
        pure_param_query = (
            (where('phase_name') == self.phase_name) & \
            (where('parameter_order') == 0) & \
            (where('parameter_type') == "G") & \
            (where('constituent_array').test(self._purity_test))
        )
        phase = dbe.phases[self.phase_name]
        param_search = dbe.search
        pure_energy_term = self.redlich_kister_sum(phase, param_search,
                                                   pure_param_query)
        return pure_energy_term / self._site_ratio_normalization

    
    def ideal_mixing_energy(self, dbe):
        #pylint: disable=W0613
        
        #Returns the ideal mixing energy in symbolic form.
        
        phase = dbe.phases[self.phase_name]
        site_ratios = self.site_ratios
        ideal_mixing_term = S.Zero
        sitefrac_limit = Float(MIN_SITE_FRACTION/10.)
        for subl_index, sublattice in enumerate(phase.constituents):
            active_comps = set(sublattice).intersection(self.components)
            ratio = site_ratios[subl_index]
            for comp in active_comps:
                sitefrac = \
                    v.SiteFraction(phase.name, subl_index, comp)
                # We lose some precision here, but this makes the limit behave nicely
                # We're okay until fractions of about 1e-12 (platform-dependent)
                mixing_term = Piecewise((sitefrac*log(sitefrac),
                                         StrictGreaterThan(sitefrac, sitefrac_limit)), (0, True),
                                        )
                ideal_mixing_term += (mixing_term*ratio)
        ideal_mixing_term *= (v.R * v.T)
        return ideal_mixing_term / self._site_ratio_normalization

    #G=A+BT

    
    def excess_mixing_energy(self, dbe):

    	None
		return excess_mixing_energy

