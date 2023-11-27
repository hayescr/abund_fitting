import numpy as np
import abund_utils as au


class ChemAbund:
    def __init__(self, abundances=None, reference='Fe', solar_reference=None):
        '''
        Takes input abundances and can convert them between different reference
        elements and solar values.

        optional additions in the future
         - multiple solar references put in for different stars?
         - a way to incorporate multiple ions i.e.,
           - inputing them
           - or knowing whether to use Fe I or Fe II to calculate xfe abundances
         - A way to input uncertainties

         What are the main operations:
         - change solar references
          - for changing solar references we should just use the logeps values
          - That means that when we initially input abundances we should calculate
           logeps
          - We would also have to change any other reference values...
         - change reference element
        '''
        if abundances is None:
            raise TypeError(
                'Please enter abundances as an input dictionary in the form {element: array}')

        self.references = set()
        self.solarref = 'None'

        if reference not in abundances and reference.lower() not in ['logeps', 'h']:
            raise IndexError(
                f'{reference} was given as the reference but is not in the provided abundances.')
        else:
            self._set_data(abundances, reference.lower())

        if reference != 'logeps' and reference != 'h':
            self._ref_to_h(reference)

        if solar_reference is not None:
            self.set_solarref(solar_reference)

        if 'Fe' in self.elements and 'fe' not in self.references:
            self._h_to_ref('Fe')

    def _set_data(self, abundances, reference):
        for element, values in abundances.items():
            if element.lower() == reference:
                setattr(self, f'{element.lower()}_h', values)
            else:
                setattr(self, f'{element.lower()}_{reference.lower()}', values)
        self.elements = abundances.keys()
        self.references.add(reference.lower())

    def _ref_to_h(self, reference):
        for element in self.elements:
            if element != reference:
                x_ref = getattr(self, f'{element.lower()}_{reference.lower()}')
                ref_h = getattr(self, f'{reference.lower()}_h')
                setattr(self, f'{element.lower()}_h', x_ref + ref_h)
        self.references.add('h')

    def _h_to_ref(self, reference):
        for element in self.elements:
            if element != reference:
                x_h = getattr(self, f'{element.lower()}_h')
                ref_h = getattr(self, f'{reference.lower()}_h')
                setattr(
                    self, f'{element.lower()}_{reference.lower()}', x_h - ref_h)
        self.references.add(reference.lower())

    def _logeps_to_h(self, solar_dict):
        for element in self.elements:
            x_logeps = getattr(self, f'{element.lower()}_logeps')
            x_solar = solar_dict[element]
            setattr(self, f'{element.lower()}_h', x_logeps - x_solar)
        self.references.add('h')

    def _h_to_logeps(self, solar_dict):
        for element in self.elements:
            x_h = getattr(self, f'{element.lower()}_h')
            x_solar = solar_dict[element]
            setattr(self, f'{element.lower()}_logeps', x_h + x_solar)
        self.references.add('logeps')

    def set_solarref(self, solar_reference):
        if solar_reference is not None:
            solar_dict = au.solar_abund(reference=solar_reference)

        if 'logeps' in self.references:
            self._logeps_to_h(solar_dict)
            runlist = []
            for ref in self.references:
                if ref != 'logeps' and ref != 'h':
                    runlist += [ref]
            for ref in runlist:
                self._h_to_ref(ref.title())
            if 'Fe' in self.elements and 'fe' not in self.references:
                self._h_to_ref('Fe')
        elif 'h' in self.references:
            self._h_to_logeps(solar_dict)

        self.solarref = solar_reference.lower()

    def calculate_xref(self, reference):
        if 'h' in self.references:
            self._h_to_ref(reference)
        else:
            print(
                f'x_{reference} abundances not calculated.  Please set a solar reference first.')

    def format_abundances(self, reference):
        if reference.lower() not in self.references:
            print(f'{reference} has not been calculated for these abundances please set a solar reference or calculate relative to this reference first.')
        else:
            output = {}
            for element in self.elements:
                if element == reference:
                    output[element] = getattr(self, f'{element.lower()}_h')
                else:
                    output[element] = getattr(
                        self, f'{element.lower()}_{reference.lower()}')

            return output
