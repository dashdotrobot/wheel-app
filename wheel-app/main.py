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
from bokeh.models import CustomJS, Range1d, ColumnDataSource, FixedTicker
from bokeh.models.widgets import *


# -------------------------------- CALLBACKS ------------------------------- #
# Define functions to update the results based on user input.                #
# -------------------------------------------------------------------------- #
def callback_update_results():
    'Build BicycleWheel object and calculate results'

    button_update.label = 'Please wait...'
    button_update.button_type = 'warning'
    status_div.text = ''

    try:
        w = build_wheel_from_UI()

        # Determine active tab
        active_tab = result_panel.active

        # Update active tab first
        result_panel_callbacks[active_tab](w)

        # Update inactive tabs
        for a, u in enumerate(result_panel_callbacks):
            if a != active_tab:
                u(w)

    except Exception as e:
        status_div.text = '<small style="color: red">Error: {:s}</small>'.format(repr(e))

    button_update.label = 'Update Results'
    button_update.button_type = 'success'

def callback_add_force():
    'Add a new force to the forces table'

    try:
        dof = force_table_src.data['dof'] + [FORCE_DOFS[f_dof.active]]
        loc = force_table_src.data['loc'] + [float(f_loc.value)]
        mag = force_table_src.data['mag'] + [float(f_mag.value)]

        force_table_src.data.update({'dof': dof, 'loc': loc, 'mag': mag})

    except Exception as e:
        status_div.text = '<small style="color: red">Error: {:s}</small>'.format(repr(e))

def callback_clear_forces():
    'Clear all forces from the table'

    force_table_src.data = {'dof': [], 'loc': [], 'mag': []}

def update_plots(wheel):
    'Update Plots tab'

    # Calculate displacements
    mm = ModeMatrix(wheel, N=SIM_OPTS_NMODES)

    # Pre-compute with and without tension effects and smeared spokes

    # Add up external forces
    F_ext = mm.F_ext(f_theta=0., f=[0., 0., 0., 0.])

    df = force_table_src.data
    for dof, loc, mag in zip(df['dof'], df['loc'], df['mag']):
        f = [0., 0., 0., 0.]
        f[FORCE_DOFS.index(dof)] = 9.81 * float(mag)
        F_ext = F_ext +\
            mm.F_ext(f_theta=float(loc) * np.pi/180, f=f)

    Bu = mm.B_theta(theta=disp_data.data['theta'], comps=[0])
    Bv = mm.B_theta(theta=disp_data.data['theta'], comps=[1])
    Bw = mm.B_theta(theta=disp_data.data['theta'], comps=[2])

    new_disp_data = \
        dict.fromkeys([u+'_'+o for o in SIM_OPTS.keys() for u in ['u', 'v', 'w']])
    new_T_data = dict.fromkeys(['dT_'+o for o in SIM_OPTS.keys()])
    for o, opts in SIM_OPTS.items():
        try:
            K = (mm.K_rim(tension=opts[0], r0=True) +
                 mm.K_spk(tension=opts[0], smeared_spokes=opts[1]))

            dm = np.linalg.solve(K, F_ext)

        except Exception as e:
            dm = np.zeros(F_ext.shape)
            status_div.text = '<small style="color: red">Error: {:s}</small>'.format(repr(e))

        new_disp_data['u_'+o] = Bu.dot(dm)*1e3
        new_disp_data['v_'+o] = Bv.dot(dm)*1e3
        new_disp_data['w_'+o] = Bw.dot(dm)*1e3

        # Calculate spoke tensions
        dT = [-s.EA/s.length/9.81 *
              np.dot(s.n, mm.B_theta(s.rim_pt[1], comps=[0, 1, 2]).dot(dm))
             for s in wheel.spokes]

        new_T_data['dT_'+o] = dT


    # Update displacements ColumnDataSource
    disp_data.data.update(new_disp_data)

    # Update grid spacing to match new number of spokes
    plot_disp.xgrid.ticker = FixedTicker(ticks=np.linspace(-np.pi, np.pi,
                                                           int(spk_num.value)+1))

    # Update spoke tensions
    theta_spk = np.array([s.rim_pt[1] for s in wheel.spokes])
    theta_spk = np.where(theta_spk <= np.pi, theta_spk, theta_spk - 2*np.pi)
    T0 = np.array([s.tension/9.81 for s in wheel.spokes])
    width = 0.5 * 2*np.pi/len(wheel.spokes) * np.ones(len(wheel.spokes))
    side = ['left' if s.hub_pt[2] > 0 else 'right' for s in wheel.spokes]
    color = ['#ff7f0e' if s == 'right' else '#1f77b4' for s in side]

    # Calculate bar heights based on sim options
    opt_str = (('T' if 0 in sim_opts.active else '_') +
               ('S' if 1 in sim_opts.active else '_'))
    if 2 in sim_opts.active:
        y = new_T_data['dT_'+opt_str]
    else:
        y = new_T_data['dT_'+opt_str] + T0

    new_T_data.update({'theta': theta_spk, 'width': width, 'side': side, 'color': color,
                       'y': y, 'T0': T0})

    T_data.data.update(new_T_data)

    # Update grid spacing to match new number of spokes
    plot_tension.xgrid.ticker = FixedTicker(ticks=np.linspace(-np.pi, np.pi,
                                                              int(spk_num.value)+1))

def update_results(wheel):
    'Update Results tab'
    result_div.text = print_wheel_info(wheel)

def build_wheel_from_UI():
    'Create a BicycleWheel object from UI inputs'

    w = BicycleWheel()

    # Hub
    w.hub = Hub(diameter=float(hub_diam.value)/1000.,
                width_nds=np.abs(float(hub_width.value[1]))/1000.,
                width_ds=np.abs(float(hub_width.value[0]))/1000.)

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

    # Drive-side
    if spk_pat_ds.value == 'Radial':
        n_cross_ds = 0
    elif spk_pat_ds.value.endswith('-cross'):
        n_cross_ds = int(spk_pat_ds.value[0])
    else:
        raise ValueError('Undefined drive-side spoke pattern: {:s}'.format(spk_pat_ds.value))

    s_matl_ds = list(SPK_MATLS)[spk_matl_ds.active]
    w.lace_cross_ds(n_spokes=int(spk_num.value)//2,
                    n_cross=n_cross_ds,
                    diameter=float(spk_diam_ds.value)/1000.,
                    young_mod=SPK_MATLS[s_matl_ds]['young_mod'],
                    offset=0.)  # Implement this later

    # Non-drive-side
    if spk_pat_nds.value == 'Radial':
        n_cross_nds = 0
    elif spk_pat_nds.value.endswith('-cross'):
        n_cross_nds = int(spk_pat_nds.value[0])
    else:
        raise ValueError('Undefined drive-side spoke pattern: {:s}'.format(spk_pat_nds.value))

    s_matl_nds = list(SPK_MATLS)[spk_matl_nds.active]
    w.lace_cross_nds(n_spokes=int(spk_num.value)//2,
                     n_cross=n_cross_nds,
                     diameter=float(spk_diam_nds.value)/1000.,
                     young_mod=SPK_MATLS[s_matl_nds]['young_mod'],
                     offset=0.)  # Implement this later


    if (spk_pat_ds.value == 'Radial' and spk_pat_nds.value == 'Radial'
        and float(spk_T_ds.value) == 0.):
        # Apply a vanishingly small tension to make stiffness matrix invertable
        w.apply_tension(T_avg=w.spokes[0].EA*1e-6)
    else:
        w.apply_tension(T_right=9.81*float(spk_T_ds.value))

    return w

callback_tension_ratio = """
    // Calculate spoke tension ratio

    R = rim_size_data.data[rim_size.value]
    r = hub_diam.value/2000

    // Get DS and NDS spoke patterns
    var n_cross_ds = 0
    if (spk_pat_ds.value == 'Radial') {
        n_cross_ds = 0
    } else {
        n_cross_ds = parseInt(spk_pat_ds.value.slice(0, 1))
    }

    var n_cross_nds = 0
    if (spk_pat_nds.value == 'Radial') {
        n_cross_nds = 0
    } else {
        n_cross_nds = parseInt(spk_pat_nds.value.slice(0, 1))
    }

    theta_h_ds = 4*3.1415/spk_num.value * n_cross_ds
    theta_h_nds = 4*3.1415/spk_num.value * n_cross_nds

    // Drive-side spoke vector
    n_ds_1 = hub_width.value[1]/1000
    n_ds_2 = R - r*Math.cos(theta_h_ds)
    n_ds_3 = r*Math.sin(theta_h_nds)
    l_ds = Math.sqrt(n_ds_1**2 + n_ds_2**2 + n_ds_3**2)

    // Non-drive-side spoke vector
    n_nds_1 = Math.abs(hub_width.value[0]/1000)
    n_nds_2 = R - r*Math.cos(theta_h_nds)
    n_nds_3 = r*Math.sin(theta_h_nds)
    l_nds = Math.sqrt(n_ds_1**2 + n_ds_2**2 + n_ds_3**2)

    c1_ds = n_ds_1 / l_ds
    c1_nds = n_nds_1 / l_nds

    spk_T_nds.value = 2*Math.round(c1_ds / c1_nds * cb_obj.value / 2)
"""

callback_sim_opts = """
    var opt_str = ((cb_obj.active.indexOf(0) >= 0 ? 'T' : '_') +
                   (cb_obj.active.indexOf(1) >= 0 ? 'S' : '_'))

    u.glyph.y.field = 'u_' + opt_str
    v.glyph.y.field = 'v_' + opt_str
    w.glyph.y.field = 'w_' + opt_str

    disp_data.change.emit()

    var data = T_data.data
    var T0 = data['T0']
    var dT = data['dT_' + opt_str]
    var y = data['y']

    for(var i=0; i < y.length; i++) {
        if (cb_obj.active.indexOf(2) >= 0) {
            y[i] = dT[i]
        } else {
            y[i] = dT[i] + T0[i]
        }
    }

    T_data.change.emit()
"""

# ---------------------------- PLOTS AND RESULTS --------------------------- #
# Define data sources and result canvases.                                   #
# -------------------------------------------------------------------------- #

# Data source to store wheel parameters for client-side calculations
RIM_SIZE_DATA = ColumnDataSource(data=dict([(s, [RIM_SIZES[s]['radius']]) for s in RIM_SIZES.keys()]))

# Panel to display text results
result_div = Div(text='Click Update Results to calculate wheel properties.', width=500)
status_div = Div(text='')

# Displacement plot
disp_names = [u+'_'+o for o in SIM_OPTS.keys() for u in ['u', 'v', 'w']]
disp_data_dict = {'theta': np.linspace(-np.pi, np.pi, 501)}
disp_data_dict.update(dict.fromkeys(disp_names, np.zeros(501)))
disp_data = ColumnDataSource(data=disp_data_dict)

plot_disp = figure(plot_height=230,
                   tools='ypan,box_zoom,reset,save',
                   tooltips=[('value', '@$name')])
plot_disp.x_range = Range1d(-np.pi, np.pi, bounds=(-np.pi, np.pi))
plot_disp.xaxis.major_tick_line_color = None
plot_disp.xaxis.minor_tick_line_color = None
plot_disp.xaxis.major_label_text_font_size = '0pt'
plot_disp.yaxis.axis_label = 'Displacement [mm]'

line_u = plot_disp.line('theta', 'u_T_', legend='lateral', name='disp_u',
                        line_width=2, color='#1f77b4', source=disp_data)
line_v = plot_disp.line('theta', 'v_T_', legend='radial', name='disp_v',
                        line_width=2, color='#ff7f0e', source=disp_data)
line_w = plot_disp.line('theta', 'w_T_', legend='tangential', name='disp_w',
                        line_width=2, color='#2ca02c', source=disp_data)

plot_disp.legend.location = 'top_left'
plot_disp.legend.click_policy = 'hide'

# Spoke tension plot
T_data_dict = {'theta': [], 'width': [], 'side': [], 'color': [], 'y': [], 'T0': []}
T_data_dict.update(dict.fromkeys(['dT_'+o for o in SIM_OPTS.keys()], []))
T_data = ColumnDataSource(data=T_data_dict)

plot_tension = figure(plot_height=230, x_range=plot_disp.x_range,
                      tools='ypan,box_zoom,reset,save',
                      tooltips=[('T', '@y{0.0} [kgf]')])
plot_tension.xaxis.major_tick_line_color = None
plot_tension.xaxis.minor_tick_line_color = None
plot_tension.xaxis.major_label_text_font_size = '0pt'
plot_tension.yaxis.axis_label = 'Spoke tension [kgf]'
bar_T = plot_tension.vbar(x='theta', top='y', color='color',
                          width='width', legend='side', source=T_data)

plot_tension.legend.location = 'bottom_left'


# ------------------------------- CONTROLLER ------------------------------- #
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
hub_symm = RadioButtonGroup(labels=['Symmetric', 'Asymmetric (e.g. rear wheel)'], active=0)
hub_width = RangeSlider(title='Rim-to-flange distance [mm]',
                        start=-50, end=50, step=1, value=(-25,25))
hub_diam = Slider(title='Flange diameter [mm]',
                  start=10, end=80, step=1, value=50)

hub_width.callback = CustomJS(args=dict(hub_symm=hub_symm), code="""
    if (hub_symm.active == 0) {
        cb_obj.value = [-cb_obj.value[1], cb_obj.value[1]]
    }
""")
hub_symm.callback = CustomJS(args=dict(hub_width=hub_width), code="""
    if (cb_obj.active == 0) {
        hub_width.value = [-hub_width.value[1], hub_width.value[1]]
    }
""")

# Create spoke controls
spk_num = Slider(title='Number of spokes',
                 start=8, end=64, step=4, value=36)

spk_matl_ds = RadioButtonGroup(labels=list(SPK_MATLS), active=0)
spk_pat_ds = Select(title='Drive-side pattern', value='3-cross',
                     options=['Radial', '1-cross', '2-cross', '3-cross', '4-cross'])
spk_diam_ds = Slider(title='Drive-side diameter [mm]',
                  start=1., end=3., step=0.1, value=2.)
spk_T_ds = Slider(title='Drive-side tension [kgf]',
                  start=0., end=200, step=2, value=100)

spk_matl_nds = RadioButtonGroup(labels=list(SPK_MATLS), active=0)
spk_pat_nds = Select(title='Non-drive-side pattern', value='3-cross',
                     options=['Radial', '1-cross', '2-cross', '3-cross', '4-cross'])
spk_diam_nds = Slider(title='Non-drive-side diameter [mm]',
                  start=1., end=3., step=0.1, value=2.)
spk_T_nds = Slider(title='Non-drive-side tension [kgf]',
                   start=0., end=200, step=2, value=100, disabled=True)

spk_T_ds.callback = CustomJS(args=dict(spk_T_nds=spk_T_nds, spk_num=spk_num,
                                       spk_pat_nds=spk_pat_nds, spk_pat_ds=spk_pat_ds,
                                       hub_width=hub_width, hub_diam=hub_diam,
                                       rim_size=rim_size, rim_size_data=RIM_SIZE_DATA),
                             code=callback_tension_ratio)

# Forces pane
force_table_src = ColumnDataSource(data=dict({'dof': ['Radial'], 'loc': [0], 'mag': [50.]}))

force_table = DataTable(source=force_table_src,
                        columns=[TableColumn(field='dof', title='DOF'),
                                 TableColumn(field='loc', title='Location [deg]'),
                                 TableColumn(field='mag', title='Magnitude [kgf]')],
                        width=270, height=120,
                        sortable=False, editable=True, reorderable=False)

f_dof = RadioButtonGroup(labels=FORCE_DOFS, active=FORCE_DOFS.index('Lateral'))
f_loc = TextInput(title='Location [deg]:', value='0')
f_mag = TextInput(title='Magnitude [kgf]:', value='0')

f_add = Button(label='Add new force')
f_add.on_click(callback_add_force)

f_clear = Button(label='Remove all forces', button_type='danger')
f_clear.on_click(callback_clear_forces)

# Computation option controls
sim_opts_list = [SIM_OPTS_TENSION, SIM_OPTS_SMEARED, SIM_OPTS_DT]
sim_opts = CheckboxButtonGroup(labels=['Tension effects', 'Smeared spokes', 'Show tension change'],
                               active=[i for i in range(len(sim_opts_list))
                                       if sim_opts_list[i]])

sim_opts.callback = CustomJS(args=dict(u=line_u, v=line_v, w=line_w, T=bar_T,
                                       disp_data=disp_data, T_data=T_data),
                             code=callback_sim_opts)

# Update Results control
button_update = Button(label='Update Results', button_type='success')
button_update.on_click(callback_update_results)


# ----------------------------- USER INTERFACE ----------------------------- #
# Define widgets, output canvases, and some native JS callbacks.             #
# -------------------------------------------------------------------------- #

rim_pane = widgetbox(rim_preset,
                     Div(text='<hr/>'),
                     rim_matl, rim_size, rim_mass,
                     rim_EI1, rim_EI2, rim_GJ)

hub_pane = widgetbox(hub_symm, hub_width, hub_diam)

spk_pane = widgetbox(spk_num,
                     Div(text='<strong>Drive-side spokes</strong>'),
                     spk_matl_ds, spk_pat_ds, spk_diam_ds, spk_T_ds,
                     Div(text='<strong>Non-drive-side spokes</strong>'),
                     spk_matl_nds, spk_pat_nds, spk_diam_nds, spk_T_nds)

force_pane = widgetbox(Div(text=('Apply multiple forces to the wheel. '
                                 'Double click on a cell to edit. Hit Enter when done.')),
                       f_clear, force_table, f_dof, f_loc, f_mag, f_add)

plot_pane = column(widgetbox(sim_opts, width=500),
                   plot_disp, plot_tension)

# Combine Wheelbuilding and Forces pane
tool_panel = Tabs(tabs=[Panel(child=rim_pane, title='Rim'),
                        Panel(child=hub_pane, title='Hub'),
                        Panel(child=spk_pane, title='Spokes'),
                        Panel(child=force_pane, title='Forces')])

result_panel = Tabs(tabs=[Panel(child=result_div, title='Results'),
                          Panel(child=plot_pane, title='Plots')])
result_panel_callbacks = [update_results, update_plots]

layout = row(column(tool_panel, button_update, status_div), result_panel, name='app')

# Render the document
curdoc().add_root(layout)
curdoc().title = 'Bicycle Wheel App'
