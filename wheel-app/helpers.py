### ---------------------------------------------------------------------- ###
# Bicycle Wheel App
# Web-based tool for calculating mechanical properties of bicycle wheels
#
# helpers.py : Script defining default values and utility functions
#
# Author     : Matthew Ford, Ph.D <www.dashdotrobot.com>
# Copyright  : Matthew Ford (2018)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# See LICENSE for a copy of the full license.
### ---------------------------------------------------------------------- ###

import numpy as np
from bikewheelcalc import *


RIM_SIZES = {'700C/29er': {'radius': 0.622/2},
             '20"':       {'radius': 0.400/2},
             '26"':       {'radius': 0.559/2},
             '27"':       {'radius': 0.630/2}}

RIM_MATLS = {'Alloy': {'young_mod': 69e9, 'shear_mod': 26e9, 'density': 2700.},
             'Steel': {'young_mod': 200e9, 'shear_mod': 77e9, 'density': 8000.}}

RIM_DEFAULT = 'Sun-Ringle CR18 700C, 36h'
RIM_PRESETS = {
    'Custom': {},
    'Alex ALX-295': {
        'matl': 0,
        'size': '700C/29er',
        'mass': 480,
        'EIrad': 310,
        'EIlat': 210,
        'GJ': 85
    },
    'DT Swiss R460': {
        'matl': 0,
        'size': '700C/29er',
        'mass': 460,
        'EIrad': 280,
        'EIlat': 230,
        'GJ': 100
    },
    'Sun-Ringle CR18 700C, 36h': {
        'matl': 0,
        'size': '700C/29er',
        'mass': 540,
        'EIrad': 110,
        'EIlat': 220,
        'GJ': 25
    },
    'Sun-Ringle CR18 20"': {
        'matl': 0,
        'size': '20"',
        'mass': 380,
        'EIrad': 100,
        'EIlat': 150,
        'GJ': 25
    },
    'Alex X404 27"': {
        'matl': 0,
        'size': '27"',
        'mass': 595,
        'EIrad': 130,
        'EIlat': 150,
        'GJ': 15
    },
    'Alex Y2000 26"': {
        'matl': 0,
        'size': '26"',
        'mass': 460,
        'EIrad': 110,
        'EIlat': 130,
        'GJ': 15
    },
    'Alex Y2000 700C': {
        'matl': 0,
        'size': '700C/29er',
        'mass': 550,
        'EIrad': 125,
        'EIlat': 160,
        'GJ': 20
    },
}

SPK_MATLS = {'Stainless steel': {'young_mod': 210e9, 'density': 8000.},
             'Alloy': {'young_mod': 69e9, 'density': 2700.}}


def avg_tension_side(w, side=1):
    'Calculate average tension for spokes on one side [N]'

    return np.mean([s.tension for s in w.spokes if side*s.hub_pt[2] > 0])

def print_wheel_info(w):
    'Return a formatted string of basic wheel properties'

    out = ''

    # Mass and inertia properties
    # TODO

    # Stiffness properties
    K_lat = calc_lat_stiff(w)
    K_lat_0 = calc_lat_stiff(w, tension=False, buckling=False)
    K_rad = calc_rad_stiff(w)
    K_tor = calc_tor_stiff(w)

    out += '<h4>Stiffness</h4>\n<table class="table">\n'
    out += '    <tr><td>Radial stiffness</td><td>{0:.1f} N/mm</td><td>({1:.1f} lbs/in)</td></tr>\n'\
        .format(K_rad/1000., K_rad * 0.0254/4.448)
    out += '    <tr><td>Lateral stiffness</td><td>{0:.1f} N/mm</td><td>({1:.1f} lbs/in)</td></tr>\n'\
        .format(K_lat/1000., K_lat * 0.0254/4.448)
    out += '    <tr><td>Torsional stiffness</td><td>{0:.2f} kN/deg</td><td>({1:.1f} lbs/deg)</td></tr>\n'\
        .format(K_tor*np.pi/180/1000., K_tor*np.pi/180/4.448)
    out += '</table>\n'

    # Tension properties
    Tc, nc = calc_buckling_tension(w, approx='linear')

    out += '<h4>Spoke Tension</h4><table class="table">\n'
    out += '    <tr><td>Average drive-side tension</td><td>{0:.1f} kgf</td><td>({1:.1f} lbs)</td></tr>\n'\
        .format(avg_tension_side(w, side=1)/9.81,
                avg_tension_side(w, side=1)/4.448)
    out += '    <tr><td>Average non-drive-side tension</td><td>{0:.1f} kgf</td><td>({1:.1f} lbs)</td></tr>\n'\
        .format(avg_tension_side(w, side=-1)/9.81,
                avg_tension_side(w, side=-1)/4.448)
    out += '    <tr><td>Maximum average tension</td><td>{0:.1f} kgf</td><td>({1:.1f} lbs)</td></tr>\n'\
        .format(Tc/9.81, Tc/4.448)
    out += '    <tr><td>Critical mode</td><td>{:d}</td></tr>\n'\
        .format(nc)
    out += '</table>\n'

    # Critical loads
    s0 = w.spokes[0]
    P_lat = K_lat * w.spokes[0].tension / (s0.EA/s0.length * s0.n[0])
    P_rad = K_rad * w.spokes[0].tension / (s0.EA/s0.length * s0.n[1])
    P_c_rad = K_lat_0*w.rim.radius / (1 + s0.EA/Tc*K_lat_0/K_rad)
    out += '<h4>Strength</h4><table class="table">\n'
    out += '    <tr><td>Radial load to buckle spokes</td><td>{0:.1f} kgf</td><td>({1:.1f} lbs)</td></tr>\n'\
        .format(P_rad / 9.81, P_rad / 4.448)
    out += '    <tr><td>Lateral load to buckle spokes</td><td>{0:.1f} kgf</td><td>({1:.1f} lbs)</td></tr>\n'\
        .format(P_lat / 9.81, P_lat / 4.448)
    out += '    <tr><td>Estimated peak radial load</td><td>{0:.1f} kgf</td><td>({1:.1f} lbs)</td></tr>\n'\
        .format(P_c_rad / 9.81, P_c_rad / 4.448)
    out += '</table>\n'


    return out
