from bokeh.io import output_file, show
from bokeh.layouts import column, row, widgetbox
from bokeh.plotting import figure
from bokeh.models.widgets import Button, RadioButtonGroup, Select, Slider, Paragraph, Div, TextInput, Panel, Tabs
from bikewheelcalc import BicycleWheel, Rim, Hub


RIM_SIZES = {'700C/29er': {'radius': 0.622/2},
             '20"':       {'radius': 0.400/2},
             '26"':       {'radius': 0.559/2},
             '27"':       {'radius': 0.630/2}}

def build_wheel_from_UI():
    'Create a BicycleWheel object from UI inputs'

    w = BicycleWheel()

    # Hub
    w.hub = Hub(diameter=hub_diam.value/1000.,
                width=hub_width/1000.,
                offset=hub_offset/1000.)

    # Rim
    w.rim = Rim(radius=RIM_SIZES[rim_size.value]['radius'],
                area=100e-6,
                I11=25./26e9, I22=200./69e9, I33=100./69e9, Iw=0.0,
                young_mod=69e9, shear_mod=26e9)

    w.lace_cross(n_spokes=36, n_cross=n_cross, diameter=1.8e-3, young_mod=210e9, offset=0.)


# Define output file
output_file('layout_test.html')


# Create rim controls
rim_matl = RadioButtonGroup(labels=['Alloy', 'Steel', 'Carbon'], active=0)
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

# Create hub controls
hub_width = Slider(title='Flange separation [mm]', start=10, end=80, step=1, value=50)
hub_diam = Slider(title='Flange diameter [mm]', start=10, end=80, step=1, value=50)
hub_offset = Slider(title='Dish offset [mm]', start=-30, end=30, step=1, value=0)

# Create spoke controls
spk_matl = RadioButtonGroup(labels=['Steel', 'Alloy', 'Carbon'], active=0)
spk_num = Slider(title='Number of spokes', start=8, end=64, step=4, value=36)
spk_diam = Slider(title='Diameter [mm]', start=1., end=3., step=0.1, value=2.)
spk_tension = Slider(title='Average spoke tension [kgf]', start=0., end=200, step=2, value=100)
spk_pattern = Select(title='Spoke pattern', value='3-cross',
                     options=['Radial', '1-cross', '2-cross', '3-cross', '4-cross'])

wheel_pane = widgetbox(Div(text='<strong>Rim</strong>'),
                       rim_matl, rim_size, rim_mass,
                       rim_EI1, rim_EI2, rim_GJ,
                       Div(text='<hr/><strong>Hub</strong>'),
                       hub_width, hub_diam, hub_offset,
                       Div(text='<hr/><strong>Spokes</strong>'),
                       spk_matl, spk_num, spk_diam, spk_tension, spk_pattern)

# Forces pane
f1_dof = RadioButtonGroup(labels=['Radial', 'Lateral', 'Tangential'], active=0)
f1_loc = TextInput(title='Location [degrees]:', value='0')
f1_mag = TextInput(title='Magnitude [N]:', value='0')
force_pane = widgetbox(Div(text='<strong>Forces</strong>'),
                       f1_dof, f1_loc, f1_mag)

# Combine Wheelbuilding and Forces pane
tool_panel = Tabs(tabs=[Panel(child=wheel_pane, title='Wheel'),
                        Panel(child=force_pane, title='Forces')])

# Result plots
p1 = figure(plot_height=250)
p2 = figure(plot_height=250)


# Render the document
show(row(tool_panel, column(Button(label='Update Results', button_type='success'),
                             p1, p2)))
