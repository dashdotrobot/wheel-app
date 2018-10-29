import numpy as np
from bikewheelcalc import *


RIM_SIZES = {'700C/29er': {'radius': 0.622/2},
             '20"':       {'radius': 0.400/2},
             '26"':       {'radius': 0.559/2},
             '27"':       {'radius': 0.630/2}}

RIM_MATLS = {'Alloy': {'young_mod': 69e9, 'shear_mod': 26e9, 'density': 2700.},
             'Steel': {'young_mod': 200e9, 'shear_mod': 77e9, 'density': 8000.}}

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
	'Calculate average tension for spokes on one side'

	return np.mean([s.tension/9.81 for s in w.spokes if side*s.hub_pt[2] > 0])

def print_wheel_info(w):
	'Return a formatted string of basic wheel properties'

	out = ''

	# Mass and inertia properties

	# Stiffness properties
	out += '<h3>Stiffness</h3>\n'
	out += '<p><strong>Radial stiffness:</strong> {:.1f} [N/mm]</p>\n'\
		.format(calc_rad_stiff(w)/1000.)
	out += '<p><strong>Lateral stiffness:</strong> {:.1f} [N/mm]</p>\n'\
		.format(calc_lat_stiff(w)/1000.)
	# out += '<p><strong>Torsional stiffness [N/mm]:</strong> {:.1f}</p>\n'\
		# .format(calc_tor_stiff(w))

	# Tension properties
	out += '<h3>Spoke Tension</h3>\n'
	out += '<p><strong>Average drive-side tension</strong>: {:.1f} [kgf]</p>\n'\
		.format(avg_tension_side(w, side=1))
	out += '<p><strong>Average non-drive-side tension</strong>: {:.1f} [kgf]</p>\n'\
		.format(avg_tension_side(w, side=-1))

	Tc, nc = calc_buckling_tension(w, approx='linear')
	out += '<p><strong>Maximum average tension:</strong> {:.1f} [kgf]</p>\n'\
		.format(Tc / 9.81)
	out += '<p><strong>Critical mode:</strong> {:d}</p>\n'\
		.format(nc)


	return out
