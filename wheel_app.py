from bokeh.io import output_file, show
from bokeh.layouts import column, row, widgetbox
from bokeh.plotting import figure
from bokeh.models.widgets import Button, RadioButtonGroup, Select, Slider, Paragraph, Div, TextInput

# Define output file
output_file('layout_test.html')


# Create rim controls
rim_matl = RadioButtonGroup(labels=['Alloy', 'Steel', 'Carbon'], active=0)
rim_size = Select(title='Wheel size', value='700C/29er',
                  options=['20"', '26"', '700C/29er', '27"'])
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
spk_matl = RadioButtonGroup(labels=['Alloy', 'Steel', 'Carbon'], active=0)
spk_num = Slider(title='Number of spokes', start=8, end=64, step=4, value=36)
spk_diam = Slider(title='Diameter [mm]', start=1., end=3., step=0.1, value=2.)
spk_tension = Slider(title='Average spoke tension [kgf]', start=0., end=200, step=2, value=100)
spk_pattern = Select(title='Spoke pattern', value='3-cross',
                     options=['Radial', '1-cross', '2-cross', '3-cross', '4-cross'])

# Widget layout
panel = widgetbox(Div(text='<strong>Rim</strong>'),
                  rim_matl, rim_size, rim_mass,
                  rim_EI1, rim_EI2, rim_GJ,
                  Div(text='<hr/><strong>Hub</strong>'),
                  hub_width, hub_diam, hub_offset,
                  Div(text='<hr/><strong>Spokes</strong>'),
                  spk_matl, spk_num, spk_diam, spk_tension, spk_pattern)


# Result plots
p1 = figure(plot_height=250)
p2 = figure(plot_height=250)


# Render the document
show(column(row(Div(text='One'), Button(label='Update Results', button_type='success')),
            row(panel, column(p1, p2))))

