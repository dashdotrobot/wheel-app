### ---------------------------------------------------------------------- ###
# Bicycle Wheel App
# Web-based tool for calculating mechanical properties of bicycle wheels
#
# main.py   : Main routine (run by Bokeh Server)
#
# Author    : Matthew Ford, Ph.D <www.dashdotrobot.com>
# Copyright : Matthew Ford (2018)
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
from bikewheelcalc import BicycleWheel, Rim, Hub, calc_lat_stiff
from helpers import *

from bokeh.layouts import column, row, widgetbox
from bokeh.plotting import figure, curdoc
from bokeh.models import CustomJS, Range1d, ColumnDataSource
from bokeh.models.widgets import *


# -------------------------------- CALLBACKS ------------------------------- #
# Define functions to update the results based on user input.                #
# -------------------------------------------------------------------------- #
def callback_update_results():
    'Build BicycleWheel object and calculate results'

    try:
        w = build_wheel_from_UI()
    except Exception as e:
        output_div.text = 'Error building wheel: {:s}'.format(repr(e))
        return False

    # Print basic wheel information
    output_div.text = print_wheel_info(w)

    # Plot results
    plot_displacements(wheel=w)

def plot_displacements(wheel):
    'Plot displacements'

    # Calculate displacements
    mm = ModeMatrix(wheel, N=int(sim_opt_nmodes.value))

    K = (mm.K_rim(tension=(0 in sim_opts.active), r0=True) +
         mm.K_spk(tension=(0 in sim_opts.active), smeared_spokes=(1 in sim_opts.active)))

    f = [0., 0., 0., 0.]
    f[int(f1_dof.active)] = float(f1_mag.value)
    F_ext = mm.F_ext(f_theta=float(f1_loc.value) * np.pi/180., f=f)

    dm = np.linalg.solve(K, F_ext)

    Bu = mm.B_theta(theta=disp_data.data['theta'], comps=[0])
    Bv = mm.B_theta(theta=disp_data.data['theta'], comps=[1])
    Bw = mm.B_theta(theta=disp_data.data['theta'], comps=[2])

    # Update ColumnDataSource
    disp_data.data.update({'disp_u': Bu.dot(dm)*1e3,
                           'disp_v': Bv.dot(dm)*1e3,
                           'disp_w': Bw.dot(dm)*1e3})

    # Update spoke tensions
    theta_spk = np.array([s.rim_pt[1] for s in wheel.spokes])
    theta_spk = np.where(theta_spk <= np.pi, theta_spk, theta_spk - 2*np.pi)
    dT = [-s.EA/s.length/9.81 *
          np.dot(s.n, mm.B_theta(s.rim_pt[1], comps=[0, 1, 2]).dot(dm))
          for s in wheel.spokes]
    T = [s.tension/9.81 + delT for s, delT in zip(wheel.spokes, dT)]
    width = 0.5 * 2*np.pi/len(wheel.spokes) * np.ones(len(wheel.spokes))
    side = ['right' if s.hub_pt[2] > 0 else 'left' for s in wheel.spokes]
    color = ['#ff7f0e' if s == 'right' else '#1f77b4' for s in side]

    T_data.data.update({'theta': theta_spk,
                        'dT': dT, 'T': T,
                        'width': width, 'side': side, 'color': color})

def build_wheel_from_UI():
    'Create a BicycleWheel object from UI inputs'

    w = BicycleWheel()

    # Hub
    w.hub = Hub(diameter=float(hub_diam.value)/1000.,
                width=float(hub_width.value)/1000.,
                offset=float(hub_offset.value)/1000.)

    # Rim
    r_matl = list(RIM_MATLS)[rim_matl.active]
    radius = RIM_SIZES[rim_size.value]['radius']
    density = RIM_MATLS[r_matl]['density']
    mass = float(rim_mass.value) / 1000.
    area = mass / (2*np.pi*radius*density)
    rim_young_mod = RIM_MATLS[r_matl]['young_mod']
    rim_shear_mod = RIM_MATLS[r_matl]['shear_mod']

    w.rim = Rim(radius=radius,
                area=area,
                I11=float(rim_GJ.value) / rim_shear_mod,
                I22=float(rim_EI2.value) / rim_young_mod,
                I33=float(rim_EI1.value) / rim_young_mod, Iw=0.0,
                young_mod=rim_young_mod, shear_mod=rim_shear_mod)

    # Spokes
    if spk_pattern.value == 'Radial':
        n_cross = 0
    elif spk_pattern.value.endswith('-cross'):
        n_cross = int(spk_pattern.value[0])
    else:
        raise ValueError('Undefined spoke pattern: {:s}'.format(spk_pattern.value))

    s_matl = list(SPK_MATLS)[spk_matl.active]
    w.lace_cross(n_spokes=int(spk_num.value),
                 n_cross=n_cross,
                 diameter=float(spk_diam.value)/1000.,
                 young_mod=SPK_MATLS[s_matl]['young_mod'],
                 offset=0.)  # Implement this later

    if spk_pattern.value == 'Radial' and float(spk_tension.value) == 0.:
        # Apply a vanishingly small tension to make stiffness matrix invertable
        w.apply_tension(w.spokes[0].EA*1e-6)
    else:
        w.apply_tension(9.81*float(spk_tension.value))

    return w


# -------------------------------- CONTROLLER ------------------------------ #
# Define widgets, output canvases, and some native JS callbacks.             #
# -------------------------------------------------------------------------- #

# Create rim controls
rim_preset = Select(title='Preset', options=list(RIM_PRESETS),
                    value=RIM_DEFAULT)

callback_reset_preset = CustomJS(args=dict(rim_preset=rim_preset), code="""
    rim_preset.value = 'Custom'
""")

rim_matl = RadioButtonGroup(labels=list(RIM_MATLS),
                            active=RIM_PRESETS[RIM_DEFAULT]['matl'],
                            callback=callback_reset_preset)
rim_size = Select(title='Wheel size', options=list(RIM_SIZES),
                  value=list(RIM_SIZES)[0],
                  callback=callback_reset_preset)
rim_mass = Slider(title='Mass [grams',
                  start=5, end=2000, step=5,
                  value=RIM_PRESETS[RIM_DEFAULT]['mass'],
                  callback=callback_reset_preset)
rim_EI1 = Slider(title='Radial stiffness [N m^2]',
                 start=10, end=1000, step=10,
                 value=RIM_PRESETS[RIM_DEFAULT]['EIrad'],
                 callback=callback_reset_preset)
rim_EI2 = Slider(title='Lateral stiffness [N m^2]',
                 start=10, end=500, step=10,
                 value=RIM_PRESETS[RIM_DEFAULT]['EIlat'],
                 callback=callback_reset_preset)
rim_GJ = Slider(title='Torsional stiffness [N m^2]',
                start=10, end=200, step=5,
                value=RIM_PRESETS[RIM_DEFAULT]['GJ'],
                callback=callback_reset_preset)

rim_preset.callback = CustomJS(args=dict(rim_matl=rim_matl,
                                         rim_mass=rim_mass,
                                         rim_size=rim_size,
                                         rim_EI1=rim_EI1,
                                         rim_EI2=rim_EI2,
                                         rim_GJ=rim_GJ,
                                         RIM_PRESETS=RIM_PRESETS),
                               code="""

    rim_preset = cb_obj.value

    if (rim_preset != 'Custom') {
        rim_matl.active = RIM_PRESETS[rim_preset]['matl']
        rim_size.value = RIM_PRESETS[rim_preset]['size']
        rim_mass.value = RIM_PRESETS[rim_preset]['mass']
        rim_GJ.value = RIM_PRESETS[rim_preset]['GJ']
        rim_EI1.value = RIM_PRESETS[rim_preset]['EIrad']
        rim_EI2.value = RIM_PRESETS[rim_preset]['EIlat']
    }
""")

# Create hub controls
hub_width = Slider(title='Flange separation [mm]',
                   start=10, end=80, step=1, value=50)
hub_diam = Slider(title='Flange diameter [mm]',
                  start=10, end=80, step=1, value=50)
hub_offset = Slider(title='Dish offset [mm]',
                    start=-25, end=25, step=1, value=0)

hub_width.callback = CustomJS(args=dict(hub_offset=hub_offset), code="""
    hub_offset.value = 0
    hub_offset.start = -Math.floor(cb_obj.value/2)
    hub_offset.end = Math.floor(cb_obj.value/2)
""")

# Create spoke controls
spk_matl = RadioButtonGroup(labels=list(SPK_MATLS), active=0)
spk_num = Slider(title='Number of spokes',
                 start=8, end=64, step=4, value=36)
spk_diam = Slider(title='Diameter [mm]',
                  start=1., end=3., step=0.1, value=2.)
spk_tension = Slider(title='Average spoke tension [kgf]',
                     start=0., end=200, step=2, value=100)
spk_pattern = Select(title='Spoke pattern', value='3-cross',
                     options=['Radial', '1-cross', '2-cross', '3-cross', '4-cross'])

# Forces pane
f1_dof = RadioButtonGroup(labels=['Lateral', 'Radial', 'Tangential'], active=1)
f1_loc = TextInput(title='Location [degrees]:', value='0')
f1_mag = TextInput(title='Magnitude [N]:', value='1000')


# Computation option controls
sim_opts = CheckboxButtonGroup(labels=['Tension effects', 'Smeared spokes'],
                               active=[0])
sim_opt_nmodes = Slider(title='Accuracy', start=2, end=40, step=1, value=20)

# Update Results control
button_update = Button(label='Update Results', button_type='success')
button_update.on_click(callback_update_results)


# ---------------------------- PLOTS AND RESULTS --------------------------- #
# Define data sources and result canvases.                                   #
# -------------------------------------------------------------------------- #

# Panel to display text results
output_div = Div(text='', width=500)

# Displacement plot
disp_data = ColumnDataSource(data={'theta': np.linspace(-np.pi, np.pi, 501),
                                   'disp_u': np.zeros(501),
                                   'disp_v': np.zeros(501),
                                   'disp_w': np.zeros(501),
                                   'disp_Rphi': np.zeros(501)})

plot_disp = figure(plot_height=240,
                   tools='ypan,box_zoom,reset,save',
                   tooltips=[('value', '@$name')])
plot_disp.x_range = Range1d(-np.pi, np.pi, bounds=(-np.pi, np.pi))
plot_disp.yaxis.axis_label = 'Displacement [mm]'

plot_disp.line('theta', 'disp_u', legend='lateral', name='disp_u',
               line_width=2, color='#1f77b4', source=disp_data)
plot_disp.line('theta', 'disp_v', legend='radial', name='disp_v',
               line_width=2, color='#ff7f0e', source=disp_data)
plot_disp.line('theta', 'disp_w', legend='tangential', name='disp_w',
               line_width=2, color='#2ca02c', source=disp_data)

plot_disp.legend.location = 'top_left'
plot_disp.legend.click_policy = 'hide'

# Spoke tension plot
T_data = ColumnDataSource(data={'theta': [], 'T': [], 'dT': [],
                                'width': [], 'side': [], 'color': []})

plot_tension = figure(plot_height=240, x_range=plot_disp.x_range,
                      tools='ypan,box_zoom,reset,save',
                      tooltips=[('T', '@T{0.0} [kgf]'),
                                ('deltaT', '@dT{+0.0} [kgf]')])
plot_tension.yaxis.axis_label = 'Spoke tension [kgf]'
plot_tension.vbar(x='theta', top='dT', color='color',
                  width='width', legend='side', source=T_data)

plot_tension.legend.location = 'bottom_left'


# ----------------------------- USER INTERFACE ----------------------------- #
# Define widgets, output canvases, and some native JS callbacks.             #
# -------------------------------------------------------------------------- #

rim_pane = widgetbox(rim_preset,
                     Div(text='<hr/>'),
                     rim_matl, rim_size, rim_mass,
                     rim_EI1, rim_EI2, rim_GJ)

hub_pane = widgetbox(hub_width, hub_diam, hub_offset)

spk_pane = widgetbox(spk_matl, spk_num, spk_diam, spk_tension, spk_pattern)

force_pane = widgetbox(f1_dof, f1_loc, f1_mag)

plot_pane = column(row(sim_opts, sim_opt_nmodes),
                   plot_disp, plot_tension)

# Combine Wheelbuilding and Forces pane
tool_panel = Tabs(tabs=[Panel(child=rim_pane, title='Rim'),
                        Panel(child=hub_pane, title='Hub'),
                        Panel(child=spk_pane, title='Spokes'),
                        Panel(child=force_pane, title='Forces')])

result_panel = Tabs(tabs=[Panel(child=output_div, title='Results'),
                          Panel(child=plot_pane, title='Plots')])

footer = Div(text="""
<small>Powered by <a href="https://bokehplots.com/">Bokeh Server</a> and <a href="https://www.python.org/">Python</a>.<br/>
Copyright Matthew Ford 2018</small>""")

layout = row(column(tool_panel, button_update, footer), result_panel)


# Initialize the plots with results
callback_update_results()

# Render the document
curdoc().add_root(layout)
curdoc().title = 'Wheel App'
