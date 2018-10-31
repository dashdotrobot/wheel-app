# Bicycle Wheel App

Interactive tool for designing virtual bicycle wheels and calculating and visualizing their mechanical behavior.

## Getting Started

Try the live version of the [Bike Wheel App hosted on heroku](https://bike-wheel-app.herokuapp.com). __BE VERY PATIENT__. It takes a long time to load the first time, and each time you press "Update Results."

## Running on your Local Machine

Alternatively, you can download and run the server from your own computer. This will run considerably faster than the live version.

### Installation

The server requires Python 3 with the following packages installed: `numpy, scipy, matplotlib, bokeh, pscript`. The easiest way is to install [Miniconda](https://conda.io/miniconda.html) and create a new environment by running

```
conda create --name wheel-app python=3.7.0 numpy scipy matplotlib bokeh pscript
conda activate wheel-app
```

Download or clone this repository to your local machine.

### Launching the app.

Navigate to the root of this repository and run
 
```
bokeh serve --show wheel-app
```

This will start up the server and automatically open a browser window to view the app.
