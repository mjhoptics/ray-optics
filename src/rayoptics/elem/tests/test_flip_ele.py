#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 09:40:53 2022

@author: Mike
"""

# %reset -sf
# %matplotlib widget
from rayoptics.environment import *

mainOpm = OpticalModel()

sm  = mainOpm['seq_model']
osp = mainOpm['optical_spec']
pm = mainOpm['parax_model']
em = mainOpm['ele_model']
pt = mainOpm['part_tree']
ar = mainOpm['analysis_results']

osp.pupil = PupilSpec(osp, key=['object', 'epd'], value=8.0)
osp.field_of_view = FieldSpec(osp, key=['object', 'height'],flds=[0,0.0,])

sm.gaps[0].thi=19.0
test_path = Path(__file__).parent.resolve()
mainOpm.add_from_file(test_path / "ACL3026U-Zemax(ZMX).zmx", t=50.)

# display the data in various ways
sm.list_model()     
print()
for k,el in enumerate(em.elements): print(k," -- ",el) 
print()
em.list_model()
print()
print(em.elements[3])

mainOpm.flip(em.elements[3])

mainOpm.update_model()
layout_plt = plt.figure(FigureClass=InteractiveLayout,do_draw_ray_fans=True, do_draw_rays=True, do_paraxial_layout=False, opt_model=mainOpm).plot()
