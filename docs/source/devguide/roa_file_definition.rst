Lens File Specification (.roa)
==============================

.. code:: ipython3

    import json
    import json_tricks
    import pprint

.. code:: ipython3

    #from rayoptics.environment import *
    from pathlib import Path
    from rayoptics.gui import roafile

Use Sasian Triplet.roa as an example
------------------------------------

.. code:: ipython3

    filename = 'Sasian Triplet.roa'

.. code:: ipython3

    file_contents = roafile.preprocess_roa(Path.cwd() / filename, roafile.module_repl_050)

The .roa file is *almost* pure json
-----------------------------------

An exception is how IEEE inf and nan are encoded.

.. code:: ipython3

    json_inputs = json.loads(file_contents)

**ray-optics** python object encoding
-------------------------------------

**ray-optics** uses the **json_tricks** package for serializing python
objects. The json_tricks package provides dump() and load() functions
that augment the default python object encoding in json by including
object metadata.

The default json encoding of a python object is to encode the python
object’s attributes as a json object. The json_tricks package encodes a
python object as a json object consisting of the object metadata info
(‘**instance_type**’) and the object ‘attributes’. The dump()
implementation recursively exports all of the object’s attributes.

Top level JSON object contains the optical model.
-------------------------------------------------

The attributes of the optical_model include the submodels (sequential,
paraxial, element, etc)

.. code:: ipython3

    pprint.pprint(json_inputs, indent=1, sort_dicts=False, depth=3)


.. parsed-literal::

    {'optical_model': {'__instance_type__': ['rayoptics.optical.opticalmodel',
                                             'OpticalModel'],
                       'attributes': {'ro_version': '0.8.5',
                                      'radius_mode': True,
                                      'specsheet': {...},
                                      'system_spec': {...},
                                      'seq_model': {...},
                                      'optical_spec': {...},
                                      'parax_model': {...},
                                      'ele_model': {...},
                                      'part_tree': {...},
                                      'profile_dict': {...},
                                      'parts_dict': {...}}}}


Optical Model attributes
------------------------

The optical model data has 3 sections:

-  the version and radius setting
-  the data for each of the submodels
-  lists of profiles and parts that are referenced in other submodels

.. code:: ipython3

    json_opm = json_inputs['optical_model']['attributes']
    pprint.pprint(json_opm, sort_dicts=False, depth=3)


.. parsed-literal::

    {'ro_version': '0.8.5',
     'radius_mode': True,
     'specsheet': {'__instance_type__': ['rayoptics.parax.specsheet', 'SpecSheet'],
                   'attributes': {'conjugate_type': 'infinite',
                                  'imager': [...],
                                  'imager_inputs': {...},
                                  'frozen_imager_inputs': [...],
                                  'etendue_inputs': {...},
                                  'etendue_values': {...}}},
     'system_spec': {'__instance_type__': ['rayoptics.optical.opticalmodel',
                                           'SystemSpec'],
                     'attributes': {'title': '',
                                    'initials': '',
                                    '_dimensions': 'mm',
                                    'temperature': 20.0,
                                    'pressure': 760.0}},
     'seq_model': {'__instance_type__': ['rayoptics.seq.sequential',
                                         'SequentialModel'],
                   'attributes': {'ifcs': [...],
                                  'gaps': [...],
                                  'z_dir': [...],
                                  'do_apertures': True,
                                  'stop_surface': 3,
                                  'cur_surface': 6}},
     'optical_spec': {'__instance_type__': ['rayoptics.raytr.opticalspec',
                                            'OpticalSpecs'],
                      'attributes': {'spectral_region': {...},
                                     'pupil': {...},
                                     'field_of_view': {...},
                                     'defocus': {...}}},
     'parax_model': {'__instance_type__': ['rayoptics.parax.paraxialdesign',
                                           'ParaxialModel'],
                     'attributes': {'seq_mapping': None,
                                    'sys': [...],
                                    'ax': [...],
                                    'pr': [...],
                                    'opt_inv': 2.274813964163764}},
     'ele_model': {'__instance_type__': ['rayoptics.elem.elements', 'ElementModel'],
                   'attributes': {}},
     'part_tree': {'__instance_type__': ['rayoptics.elem.parttree', 'PartTree'],
                   'attributes': {'root_node': {...}}},
     'profile_dict': {'5850585344': {'__instance_type__': [...],
                                     'attributes': {...}},
                      '5850592496': {'__instance_type__': [...],
                                     'attributes': {...}},
                      '5850589232': {'__instance_type__': [...],
                                     'attributes': {...}},
                      '5858434144': {'__instance_type__': [...],
                                     'attributes': {...}},
                      '5858435344': {'__instance_type__': [...],
                                     'attributes': {...}},
                      '5858436976': {'__instance_type__': [...],
                                     'attributes': {...}},
                      '5858433184': {'__instance_type__': [...],
                                     'attributes': {...}},
                      '5850587312': {'__instance_type__': [...],
                                     'attributes': {...}}},
     'parts_dict': {'5850585824': {'__instance_type__': [...], 'attributes': {...}},
                    '5850587168': {'__instance_type__': [...], 'attributes': {...}},
                    '5850593504': {'__instance_type__': [...], 'attributes': {...}},
                    '5850585680': {'__instance_type__': [...], 'attributes': {...}},
                    '5850593936': {'__instance_type__': [...], 'attributes': {...}},
                    '5850593264': {'__instance_type__': [...], 'attributes': {...}},
                    '5850590864': {'__instance_type__': [...], 'attributes': {...}},
                    '5850584144': {'__instance_type__': [...], 'attributes': {...}},
                    '5858442784': {'__instance_type__': [...],
                                   'attributes': {...}}}}


The Sequential Model
--------------------

Only the ``ifcs``, ``gaps``, and ``z_dir`` arrays are needed to fully
specify the model.

.. code:: ipython3

    json_sm = json_opm['seq_model']['attributes']
    pprint.pprint(json_sm, sort_dicts=False, depth=2)


.. parsed-literal::

    {'ifcs': [{...}, {...}, {...}, {...}, {...}, {...}, {...}, {...}],
     'gaps': [{...}, {...}, {...}, {...}, {...}, {...}, {...}],
     'z_dir': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
     'do_apertures': True,
     'stop_surface': 3,
     'cur_surface': 6}


Surface encoding example
~~~~~~~~~~~~~~~~~~~~~~~~

Note the profile is encoded in a separate json object,
json_opm[‘profile_dict’]. The ‘profile_id’ attribute is a key to the
actual profile, contained in the json_opm[‘profile_dict’].

.. code:: ipython3

    json_ifcs2 = json_sm['ifcs'][2]
    pprint.pprint(json_ifcs2, sort_dicts=False, depth=3)


.. parsed-literal::

    {'__instance_type__': ['rayoptics.elem.surface', 'Surface'],
     'attributes': {'interact_mode': 'transmit',
                    'delta_n': -0.6910020663241183,
                    'decenter': None,
                    'max_aperture': 8.948204697566771,
                    'label': '',
                    'clear_apertures': [],
                    'edge_apertures': [],
                    'profile_id': '5850589232'}}


Gap encoding example
~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    json_gaps3 = json_sm['gaps'][3]
    pprint.pprint(json_gaps3, sort_dicts=False, depth=6)


.. parsed-literal::

    {'__instance_type__': ['rayoptics.seq.gap', 'Gap'],
     'attributes': {'thi': 0.975,
                    'medium': {'__instance_type__': ['opticalglass.schott',
                                                     'SchottGlass'],
                               'slots': {},
                               'attributes': {'gname': 'N-SF5',
                                              'coefs': {'__ndarray__': [1.52481889,
                                                                        0.187085527,
                                                                        1.42729015,
                                                                        0.011254756,
                                                                        0.0588995392,
                                                                        129.141675],
                                                        'dtype': 'float64',
                                                        'shape': [6]}}}}}


Optical Spec
------------

.. code:: ipython3

    json_osp = json_opm['optical_spec']['attributes']
    pprint.pprint(json_osp, sort_dicts=False, depth=1)


.. parsed-literal::

    {'spectral_region': {...},
     'pupil': {...},
     'field_of_view': {...},
     'defocus': {...}}


spectral_region
~~~~~~~~~~~~~~~

.. code:: ipython3

    json_wvls = json_osp['spectral_region']['attributes']
    pprint.pprint(json_wvls, sort_dicts=False, depth=2)


.. parsed-literal::

    {'wavelengths': [486.1327, 587.5618, 656.2725],
     'spectral_wts': [0.5, 1.0, 0.5],
     'render_colors': ['#268bd2', '#859900', '#dc322f'],
     'reference_wvl': 1,
     'coating_wvl': 550.0}


pupil
~~~~~

.. code:: ipython3

    json_pupil = json_osp['pupil']['attributes']
    pprint.pprint(json_pupil, sort_dicts=False, depth=3)


.. parsed-literal::

    {'key': ['aperture', 'object', 'pupil'],
     'value': 12.5,
     'pupil_rays': [[0.0, 0.0], [1.0, 0.0], [-1.0, 0.0], [0.0, 1.0], [0.0, -1.0]],
     'ray_labels': ['00', '+X', '-X', '+Y', '-Y']}


field of view
~~~~~~~~~~~~~

.. code:: ipython3

    json_fov = json_osp['field_of_view']['attributes']
    pprint.pprint(json_fov, sort_dicts=False, depth=2)


.. parsed-literal::

    {'key': ['field', 'object', 'angle'],
     'value': 20.0,
     'is_relative': True,
     'fields': [{...}, {...}, {...}],
     'index_labels': ['axis', ' 0.71y', 'edge']}


individual field point
~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    json_field2 = json_fov['fields'][2]
    pprint.pprint(json_field2, sort_dicts=False, depth=2)


.. parsed-literal::

    {'__instance_type__': ['rayoptics.raytr.opticalspec', 'Field'],
     'attributes': {'x': 0.0,
                    'y': 1.0,
                    'vux': 0.0,
                    'vuy': 0.0,
                    'vlx': 0.0,
                    'vly': 0.0,
                    'wt': 1.0,
                    'aim_pt': {...}}}


defocus
~~~~~~~

.. code:: ipython3

    json_defocus = json_osp['defocus']['attributes']
    pprint.pprint(json_defocus, sort_dicts=False, depth=1)


.. parsed-literal::

    {'focus_shift': 0.0, 'defocus_range': 0.0}


Paraxial Model
--------------

.. code:: ipython3

    json_pm = json_opm['parax_model']['attributes']
    pprint.pprint(json_pm, sort_dicts=False, depth=3)


.. parsed-literal::

    {'seq_mapping': None,
     'sys': [[0.0, 10000000000.0, 1.0, 'dummy'],
             [0.029140221242530184,
              2.856885923564585,
              1.6910020663241183,
              'transmit'],
             [-9.425384275234017e-05, 5.86, 1.0, 'transmit'],
             [-0.027506830028392514,
              0.5828874868684721,
              1.6727070351743674,
              'transmit'],
             [-0.03072282769338543, 4.822, 1.0, 'transmit'],
             [0.007964615386577972,
              1.8491993961884612,
              1.6910020663241183,
              'transmit'],
             [0.03371695730129101, 41.2365, 1.0, 'transmit'],
             [0.0, 0.0, 1.0, 'dummy']],
     'ax': [[0.0, 6.249999992700494e-10],
            [6.249999992700494, -0.18212638192810443],
            [5.729685695860345, -0.1815863370335065],
            [4.665589760843997, -0.05325075249976219],
            [4.634550563545556, 0.08913574590033027],
            [5.0643631302769485, 0.0488000413897083],
            [5.154604137348769, -0.12499752621433813],
            [0.000143647611214881, -0.12499752621433813]],
     'pr': [[-3639702346.9129076, 0.36397023426620234],
            [-4.250884056091309, 0.4878419361370471],
            [-2.8571752958168855, 0.48757263638599935],
            [3.5340507098524654e-07, 0.4875726461070525],
            [0.28420034776022174, 0.49630408442169],
            [2.677378642841611, 0.47497979328721845],
            [3.5557109897900556, 0.355092037668736],
            [18.198463801116887, 0.355092037668736]],
     'opt_inv': 2.274813964163764}


Profile dictionary
------------------

.. code:: ipython3

    json_profiles = json_opm['profile_dict']
    pprint.pprint(json_profiles, sort_dicts=False, depth=3)


.. parsed-literal::

    {'5850585344': {'__instance_type__': ['rayoptics.elem.profiles', 'Spherical'],
                    'attributes': {'cv': 0.0}},
     '5850592496': {'__instance_type__': ['rayoptics.elem.profiles', 'Spherical'],
                    'attributes': {'cv': 0.042170961076202926}},
     '5850589232': {'__instance_type__': ['rayoptics.elem.profiles', 'Spherical'],
                    'attributes': {'cv': 0.00013640168003221264}},
     '5858434144': {'__instance_type__': ['rayoptics.elem.profiles', 'Spherical'],
                    'attributes': {'cv': -0.04088976120379457}},
     '5858435344': {'__instance_type__': ['rayoptics.elem.profiles', 'Spherical'],
                    'attributes': {'cv': 0.04567044208987943}},
     '5858436976': {'__instance_type__': ['rayoptics.elem.profiles', 'Spherical'],
                    'attributes': {'cv': 0.011526181721781025}},
     '5858433184': {'__instance_type__': ['rayoptics.elem.profiles', 'Spherical'],
                    'attributes': {'cv': -0.04879429301948844}},
     '5850587312': {'__instance_type__': ['rayoptics.elem.profiles', 'Spherical'],
                    'attributes': {'cv': 0.0}}}


sample profile instance
~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    json_profile2 = json_profiles['5850589232']
    pprint.pprint(json_profile2, sort_dicts=False, depth=3)


.. parsed-literal::

    {'__instance_type__': ['rayoptics.elem.profiles', 'Spherical'],
     'attributes': {'cv': 0.00013640168003221264}}


Element Model
-------------

empty by design

.. code:: ipython3

    json_em = json_opm['ele_model']['attributes']
    pprint.pprint(json_em, sort_dicts=False, depth=2)


.. parsed-literal::

    {}


Part Tree
---------

.. code:: ipython3

    json_pt = json_opm['part_tree']['attributes']
    pprint.pprint(json_pt, sort_dicts=False, depth=3)


.. parsed-literal::

    {'root_node': {'id_key': '5850583472',
                   'tag': '#group#root',
                   'name': 'root',
                   'children': [{...},
                                {...},
                                {...},
                                {...},
                                {...},
                                {...},
                                {...},
                                {...},
                                {...}]}}


Part dictionary
---------------

.. code:: ipython3

    json_parts = json_opm['parts_dict']
    pprint.pprint(json_parts, sort_dicts=False, depth=2)


.. parsed-literal::

    {'5850585824': {'__instance_type__': [...], 'attributes': {...}},
     '5850587168': {'__instance_type__': [...], 'attributes': {...}},
     '5850593504': {'__instance_type__': [...], 'attributes': {...}},
     '5850585680': {'__instance_type__': [...], 'attributes': {...}},
     '5850593936': {'__instance_type__': [...], 'attributes': {...}},
     '5850593264': {'__instance_type__': [...], 'attributes': {...}},
     '5850590864': {'__instance_type__': [...], 'attributes': {...}},
     '5850584144': {'__instance_type__': [...], 'attributes': {...}},
     '5858442784': {'__instance_type__': [...], 'attributes': {...}}}


sample Part, lens element #1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    json_part2 = json_parts['5850593504']
    pprint.pprint(json_part2, sort_dicts=False, depth=6)


.. parsed-literal::

    {'__instance_type__': ['rayoptics.elem.elements', 'Element'],
     'attributes': {'label': 'E1',
                    'tfrm': [{'__ndarray__': [[1.0, 0.0, 0.0],
                                              [0.0, 1.0, 0.0],
                                              [0.0, 0.0, 1.0]],
                              'dtype': 'float64',
                              'shape': [3, 3],
                              'Corder': True},
                             {'__ndarray__': [0.0, 0.0, 0.0],
                              'dtype': 'float64',
                              'shape': [3]}],
                    's1_indx': 1,
                    's2_indx': 2,
                    'medium_name': 'N-LAK9',
                    '_sd': 10.008722902253046,
                    'hole_sd': None,
                    'flat1': None,
                    'flat2': 8.948204697566771,
                    'do_flat1': 'if concave',
                    'do_flat2': 'if concave',
                    'edge_extent': [-10.008722902253046, 10.008722902253046],
                    'profile1_id': '5850592496',
                    'profile2_id': '5850589232'}}




