import numpy as np
from bikewheelcalc import BicycleWheel, Rim, Hub, calc_lat_stiff
from helpers import *

from bokeh.io import output_file, show
from bokeh.layouts import column, row, widgetbox
from bokeh.plotting import figure, curdoc
from bokeh.models.widgets import Button, RadioButtonGroup, Select, Slider, Paragraph, Div, TextInput, Panel, Tabs


RIM_SIZES = {'700C/29er': {'radius': 0.622/2},
             '20"':       {'radius': 0.400/2},
             '26"':       {'radius': 0.559/2},
             '27"':       {'radius': 0.630/2}}

RIM_MATLS = {'Alloy': {'young_mod': 69e9, 'shear_mod': 26e9, 'density': 2700.},
             'Steel': {'young_mod': 200e9, 'shear_mod': 77e9, 'density': 8000.}}

SPK_MATLS = {'Stainless steel': {'young_mod': 210e9, 'density': 8000.},
             'Alloy': {'young_mod': 69e9, 'density': 2700.}}


# Callbacks
def callback_Update_Results():
    'Build BicycleWheel object and calculate results'

    try:
        w = build_wheel_from_UI()
    except Exception as e:
        output_div.text = 'Error building wheel: {:s}'.format(repr(e))
        return False

    # Print basic wheel information
    output_div.text = print_wheel_info(w)


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


# Define output file

# Create presets controls
# TODO

# Create rim controls
rim_matl = RadioButtonGroup(labels=list(RIM_MATLS), active=0)
rim_size = Select(title='Wheel size', value=list(RIM_SIZES)[0],
                  options=list(RIM_SIZES))
rim_mass = Slider(title='Mass [grams',
                  start=5, end=2000, value=500, step=5)
rim_EI1 = Slider(title='Radial stiffness [N m^2]',
                 start=10, end=1000, value=100, step=10)
rim_EI2 = Slider(title='Lateral stiffness [N m^2]',
                 start=10, end=500, value=100, step=10)
rim_GJ = Slider(title='Torsional stiffness [N m^2]',
                start=10, end=500, value=100, step=10)

rim_pane = widgetbox(Div(text='<strong>Rim</strong>'),
                     rim_matl, rim_size, rim_mass,
                     rim_EI1, rim_EI2, rim_GJ)

# Create hub controls
hub_width = Slider(title='Flange separation [mm]', start=10, end=80, step=1, value=50)
hub_diam = Slider(title='Flange diameter [mm]', start=10, end=80, step=1, value=50)
hub_offset = Slider(title='Dish offset [mm]', start=-30, end=30, step=1, value=0)

hub_pane = widgetbox(Div(text='<hr/><strong>Hub</strong>'),
                     hub_width, hub_diam, hub_offset)

# Create spoke controls
spk_matl = RadioButtonGroup(labels=list(SPK_MATLS), active=0)
spk_num = Slider(title='Number of spokes', start=8, end=64, step=4, value=36)
spk_diam = Slider(title='Diameter [mm]', start=1., end=3., step=0.1, value=2.)
spk_tension = Slider(title='Average spoke tension [kgf]', start=0., end=200, step=2, value=100)
spk_pattern = Select(title='Spoke pattern', value='3-cross',
                     options=['Radial', '1-cross', '2-cross', '3-cross', '4-cross'])

spk_pane = widgetbox(Div(text='<hr/><strong>Spokes</strong>'),
                     spk_matl, spk_num, spk_diam, spk_tension, spk_pattern)

# Forces pane
f1_dof = RadioButtonGroup(labels=['Radial', 'Lateral', 'Tangential'], active=0)
f1_loc = TextInput(title='Location [degrees]:', value='0')
f1_mag = TextInput(title='Magnitude [N]:', value='0')
force_pane = widgetbox(Div(text='<strong>Forces</strong>'),
                       f1_dof, f1_loc, f1_mag)

# Combine Wheelbuilding and Forces pane
tool_panel = Tabs(tabs=[Panel(child=rim_pane, title='Rim'),
                        Panel(child=hub_pane, title='Hub'),
                        Panel(child=spk_pane, title='Spokes'),
                        Panel(child=force_pane, title='Forces')])


# Update Results control
button_update = Button(label='Update Results', button_type='success')
button_update.on_click(callback_Update_Results)

# Text results
output_div = Div(text='>>>')
text_pane = column(Div(text='<strong>Console</strong>'),
                   output_div)

# Plot results
p1 = figure(plot_height=250)
p2 = figure(plot_height=250)
plot_pane = column(p1, p2)

result_panel = Tabs(tabs=[Panel(child=text_pane, title='Results'),
                    Panel(child=plot_pane, title='Plots')])


# Render the document
layout = row(column(tool_panel,button_update), result_panel)

curdoc().add_root(layout)
