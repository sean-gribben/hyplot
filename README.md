# Hyplot

A simple plotting tool for hysplit trajectories. Uses cartopy for plotting functionality. This tool accepts hysplit trajecotry files and parses the coordinate data to determine the trajectories altitude and lat/lon.

## Structure overview

### abstract_plot

This is an abstract class inherited by **traj_plot** and **freq_plot**. It essentially reads trajecotry files and stores the data from such in a `self.trajectories` property, where each time step of the trajectory is stored as a dictionary of `lat`, `lon`, `altitude`, `pressure` and `index`. 

When adding a trajectory, the `append` function is used and the name of the trajectory file specified.

This class also includes useful utility functions used for plotting such as `determine_extent`

### traj_plot

This class plots individual trajectories on a world map. If no axes is specified, hyplot will generate a Plate Caree projection plot and plots the trajectories on it, otherwise it is up to the user to create an axes with the appropriate projection.

Additionally, an altitude plot can be created as a time series plot of the altitude aver the span over the trajectory.

### freq_plot

This class is used to create KDE (Kernel Density Estimation) AKA contour plots of the trajectory data. There are a few options here:

 - **contour_plot**: This will plot the trajectories on a world map as a representation of the density of the trajectories. This is performed by counting the collisions in cells at the specified resolution (default 0.5 deg) and using this data as weights to a KDE calculated by `scipy`. It is recomended you set the `manual_bandwidth` if you are doing multipler contnour plots of related, but different data, as the bandwidth used for the `scipy.guassian_kde` function will change between datasets, which effects how the density map is calculated. [For an example of why this might be a problem, see here](https://stackoverflow.com/questions/63068193/scipy-gaussian-kde-produces-different-results-depending-on-method-used-weights)
  - **altitude_contour_plot**: This will create a contour plot of the altitudes for all trajectories. This will be plotted as a standard 2d plot (not on any map)
  - **pressure_contour_plot**: Same as the **altitude_contour_plot** but with pressure instead.
  - **kml_contour_plot**: This will create a contour plot using the same method as **contour_plot** (NOTE: Not updated to support weights and manual_bandwidth), but instead output the results to a .kml file (a file that can be read by Google Earth)