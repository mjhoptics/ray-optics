import itertools
from . import globalspec
from . import surface as srf
from . import profiles as pr
from . import gap
from . import medium as m
from . import raytrace as rt
from . import transform as trns
from glass import glassfactory as gfact
import numpy as np
import transforms3d as t3d
from math import sqrt


class SequentialModel:

    def __init__(self):
        self.radius_mode = False
        self.global_spec = globalspec.GlobalData()
        self.surfs = []
        self.gaps = []
        self.stop_surface = 1
        self.cur_surface = 0
        self.surfs.append(srf.Surface('Obj'))
        self.surfs.append(srf.Surface('Img'))
        self.gaps.append(gap.Gap())

    def reset(self):
        rdm = self.radius_mode
        self.__init__()
        self.radius_mode = rdm

    def get_num_surfaces(self):
        return len(self.surfs)

    def get_surface_and_gap(self, srf=None):
        if not srf:
            srf = self.cur_surface
        s = self.surfs[srf]
        if srf == len(self.gaps):
            g = None
        else:
            g = self.gaps[srf]
        return s, g

    def set_cur_surface(self, s):
        self.cur_surface = s

    def insert(self, surf, gap):
        """ insert surf and gap at the cur_gap edge of the sequential model
            graph """
        self.cur_surface += 1
        self.surfs.insert(self.cur_surface, surf)
        self.gaps.insert(self.cur_surface, gap)

    def add_surface(self, surf):
        """ add a surface where surf is a list that contains:
            [curvature, thickness, refractive_index, v-number] """
        s, g = self.create_surface_and_gap(surf)
        self.insert(s, g)

    def update_surface(self, surf):
        """ add a surface where surf is a list that contains:
            [curvature, thickness, refractive_index, v-number] """

        self.insert(*self.create_surface_and_gap(surf))

    def update_model(self):
        self.global_spec.update_model(self)
        self.set_clear_apertures()

    def insert_surface_and_gap(self):
        s = srf.Surface()
        g = gap.Gap()
        self.insert(s, g)
        return s, g

    def update_surface_and_gap_from_cv_input(self, dlist, idx=None):
        if isinstance(idx, int):
            s, g = self.get_surface_and_gap(idx)
        else:
            s, g = self.insert_surface_and_gap()

        if self.radius_mode:
            if dlist[0] != 0.0:
                s.profile.cv = 1.0/dlist[0]
            else:
                s.profile.cv = 0.0
        else:
            s.profile.cv = dlist[0]

        if g:
            g.thi = dlist[1]
            if len(dlist) < 3:
                g.medium = m.Air()
            elif dlist[2].upper() == 'REFL':
                s.refract_mode = 'REFL'
                g.medium = self.gaps[self.cur_surface-1].medium
            elif dlist[2].isnumeric():
                n, v = m.glass_decode(dlist[2])
                g.medium = m.Glass(n, v, '')
            else:
                name, cat = dlist[2].split('_')
                if cat.upper() == 'SCHOTT' and name[:1].upper() == 'N':
                    name = name[:1]+'-'+name[1:]
                g.medium = gfact.create_glass(cat, name)
        else:
            # at image surface, apply defocus to previous thickness
            self.gaps[idx-1].thi += dlist[1]

    def update_surface_profile_cv_input(self, profile_type, idx=None):
        if not isinstance(idx, int):
            idx = self.cur_surface
        cur_profile = self.surfs[idx].profile
        new_profile = pr.mutate_profile(cur_profile, profile_type)
        self.surfs[idx].profile = new_profile
        return self.surfs[idx].profile

    def update_surface_decenter_cv_input(self, idx=None):
        if not isinstance(idx, int):
            idx = self.cur_surface
        if not self.surfs[idx].decenter:
            self.surfs[idx].decenter = srf.DecenterData()
        return self.surfs[idx].decenter

    def create_surface_and_gap(self, surf):
        """ create a surface and gap where surf is a list that contains:
            [curvature, thickness, refractive_index, v-number] """
        s = srf.Surface()
        s.profile.cv = surf[0]

        if surf[3] == 0.0:
            if surf[2] == 1.0:
                mat = m.Air()
            else:
                mat = m.Medium(surf[2])
        else:
            mat = m.Glass(surf[2], surf[3], '')

        g = gap.Gap(surf[1], mat)
        return s, g

    def list_model(self):
        for i, sg in enumerate(itertools.zip_longest(self.surfs, self.gaps)):
            if sg[1]:
                print(i, sg[0])
                print('    ', sg[1])
            else:
                print(i, sg[0])

    def list_gaps(self):
        for i, gp in enumerate(self.gaps):
            print(i, gp)

    def list_surfaces(self):
        for i, s in enumerate(self.surfs):
            print(i, s)

    def list_surface_and_gap(self, i):
        s = self.surfs[i]
        cvr = s.profile.cv
        if self.radius_mode:
            if cvr != 0.0:
                cvr = 1.0/cvr
        sd = s.surface_od()

        if i < len(self.gaps):
            g = self.gaps[i]
            thi = g.thi
            med = g.medium.name()
        else:
            thi = ''
            med = ''
        return [cvr, thi, med, sd]

    def list_decenters(self):
        for i, sg in enumerate(itertools.zip_longest(self.surfs, self.gaps)):
            if sg[1]:
                print(i, sg[1])
                if sg[0].decenter is not None:
                    print(' ', repr(sg[0].decenter))
            else:
                if sg[0].decenter is not None:
                    print(i, repr(sg[0].decenter))

    def list_elements(self):
        for i, gp in enumerate(self.gaps):
            if gp.medium.label.lower() != 'air':
                print(self.surfs[i].profile,
                      self.surfs[i+1].profile,
                      gp)

    def trace_boundary_rays(self):
        pupil_rays = [[0., 0.], [1., 0.], [-1., 0.], [0., 1.], [0., -1.]]
        rayset = []
        fov = self.global_spec.field_of_view
        for fi, f in enumerate(fov.fields):
            rim_rays = []
            for r in pupil_rays:
                ray, op = self.global_spec.trace(self, r, fi)
                rim_rays.append([ray, op])
            rayset.append(rim_rays)
        return rayset

    def shift_start_of_rayset(self, rayset, start_offset):
        """ start_offset is positive if to left of first surface """
        s1 = self.surfs[1]
        s0 = self.surfs[0]
        g0 = gap.Gap(start_offset, self.gaps[0].medium)
        r, t = trns.reverse_transform(s1, g0, s0)
        for fi, f in enumerate(rayset):
            for ri, ray in enumerate(f):
                b4_pt = r.dot(ray[0][1][0]) + t
                b4_dir = r.dot(ray[0][0][1])
                if ri == 0:
                    # For the chief ray, use the input offset.
                    dst = -start_offset
                else:
                    pt0 = f[0][0][0][0]
                    dir0 = f[0][0][0][1]
                    # Calculate distance along ray to plane perpendicular to
                    #  the chief ray.
                    dst = -(b4_pt - pt0).dot(dir0)/b4_dir.dot(dir0)
                pt = b4_pt + dst*b4_dir
#                print("fld:", fi, "ray:", ri, dst, pt)
                ray[0][0][0] = pt
                ray[0][0][1] = b4_dir
        return r, t

    def set_clear_apertures(self):
        rayset = self.trace_boundary_rays()
        for i, s in enumerate(self.surfs):
            max_ap = -1.0e+10
            for f in rayset:
                for p in f:
                    ap = sqrt(p[0][i][0][0]**2 + p[0][i][0][1]**2)
                    if ap > max_ap:
                        max_ap = ap
            cir_ap = srf.Circular(max_ap)
            if len(s.clear_apertures):
                s.clear_apertures[0] = cir_ap
            else:
                s.clear_apertures.append(cir_ap)

    def trace(self, pt0, dir0, wl, eps=1.0e-12):
        path = itertools.zip_longest(self.surfs, self.gaps)
        return rt.trace(path, pt0, dir0, wl, eps)

    def compute_global_coords(self, glo=1):
        """ Return global surface coordinates (rot, t) wrt surface glo. """
        tfrms = []
        r, t = np.identity(3), np.array([0., 0., 0.])
        prev = r, t
        tfrms.append(prev)
#        print(glo, t, *np.rad2deg(t3d.euler.mat2euler(r)))
        if glo > 0:
            # iterate in reverse over the segments before the
            #  global reference surface
            go = glo
            path = itertools.zip_longest(self.surfs[glo::-1],
                                         self.gaps[glo-1::-1])
            after = next(path)
            # loop of remaining surfaces in path
            while True:
                try:
                    before = next(path)
                    go -= 1
                    r, t = trns.reverse_transform(before[0], after[1],
                                                  after[0])
                    t = prev[0].dot(t) + prev[1]
                    r = prev[0].dot(r)
#                    print(go, t,
#                          *np.rad2deg(trns.euler2opt(t3d.euler.mat2euler(r))))
                    prev = r, t
                    tfrms.append(prev)
                    after = before
                except StopIteration:
                    break
            tfrms.reverse()
        path = itertools.zip_longest(self.surfs[glo:], self.gaps[glo:])
        before = next(path)
        prev = np.identity(3), np.array([0., 0., 0.])
        go = glo
        # loop forward over the remaining surfaces in path
        while True:
            try:
                after = next(path)
                go += 1
                r, t = trns.forward_transform(before[0], before[1], after[0])
                t = prev[0].dot(t) + prev[1]
                r = prev[0].dot(r)
#                print(go, t,
#                      *np.rad2deg(trns.euler2opt(t3d.euler.mat2euler(r))))
                prev = r, t
                tfrms.append(prev)
                before = after
            except StopIteration:
                break
        return tfrms

    def compute_global_coords_wrt(self, rng, go=1):
        for s in rng:
            pass
