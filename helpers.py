import numpy as np
from bikewheelcalc import *

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
