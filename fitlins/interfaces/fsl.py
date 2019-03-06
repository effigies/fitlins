from nipype.interfaces.base import BaseInterfaceInputSpec, TraitedSpec, SimpleInterface, traits


class FSLLevel1ShimInputSpec(BaseInterfaceInputSpec):
    session_info = traits.List(traits.Dict())
    contrast_info = traits.List(traits.List(traits.Dict))


class FSLLevel1ShimOutputSpec(TraitedSpec):
    interscan_interval = traits.Float
    session_info = traits.List(traits.Dict)
    bases = traits.Any
    orthogonalization = traits.Dict
    contrasts = traits.List(
        traits.Either(
            traits.Tuple(traits.Str, traits.Enum('T'), traits.List(traits.Str),
                         traits.List(traits.Float)),
            traits.Tuple(traits.Str, traits.Enum('T'), traits.List(traits.Str),
                         traits.List(traits.Float), traits.List(traits.Float)),
            traits.Tuple(traits.Str, traits.Enum('F'),
                         traits.List(
                             traits.Either(
                                 traits.Tuple(traits.Str, traits.Enum('T'),
                                              traits.List(traits.Str),
                                              traits.List(traits.Float)),
                                 traits.Tuple(traits.Str, traits.Enum('T'),
                                              traits.List(traits.Str),
                                              traits.List(traits.Float),
                                              traits.List(traits.Float)))))),
        desc="List of contrasts with each contrast being a list of the form - \
[('name', 'stat', [condition list], [weight list], [session list])]. if \
session list is None or not provided, all sessions are used. For F \
contrasts, the condition list should contain previously defined \
T-contrasts.")


class FSLLevel1Shim(SimpleInterface):
    input_spec = FSLLevel1ShimInputSpec
    output_spec = FSLLevel1ShimOutputSpec

    def _run_interface(self, runtime):
        res = [self._load_info(ses_info, con_info)
               for ses_info, con_info in zip(self.inputs.session_info,
                                             self.inputs.contrast_info)]

        isi = {subres['interscan_interval'] for subres in res}
        if len(isi) > 0:
            raise ValueError("Non-constant inter-scan interval across runs")
        self._results['interscan_interval'] = isi.pop()

        self._results['session_info'] = [subres['session_info'] for subres in res]

        contrasts = res[0]['contrasts']
        for subres in res[1:]:
            if subres['contrasts'] != contrasts:
                raise ValueError("Contrasts varying across runs not yet supported")
        self._results['contrasts'] = contrasts

        return runtime

    @staticmethod
    def _load_info(bold_file, session_info, contrast_info):
        import pandas as pd
        results = {'interscan_interval': session_info['repetition_time']}

        out_info = {}
        out_info['scans'] = bold_file
        out_info['cond'] = []
        # XXX: Need to create a model that produces sparse events to test this
        # if session_info.get('sparse'):
        #     sparse = pd.read_hdf(session_info['sparse'], 'sparse')
        #     ...
        out_info['regress'] = []
        if session_info.get('dense'):
            dense = pd.read_hdf(session_info['dense'], 'dense')
            out_info['regress'] = [{'name': name, 'val': val}
                                   for name, val in dense.to_dict('list')]

        results['info'] = out_info

        out_contrasts = []
        for contrast in contrast_info:
            if contrast['type'] == 't':
                conds = ()
                weights = ()
                for cond, weight in sorted(contrast['weights'].items()):
                    conds.append(cond)
                    weights.append(weight)
                out_contrasts.append((contrast['name'], 'T', conds, weights))
            else:
                raise ValueError(f"Unimplemented contrast type {contrast['type']!r}")

        results['contrasts'] = out_contrasts

        return results
