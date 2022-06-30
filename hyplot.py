import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import datetime as dt
import os
from math import inf, log
import numpy as np
import simplekml
from scipy.stats import gaussian_kde
from matplotlib.colors import ListedColormap
from random import random
from math import exp, pi, log1p, log
import warnings


# WARNING: THIS SCRIPT EXPECTS BACKWARDS TRAJECTORIES ONLY (for now). FORWARDS TRAJECTORIES WILL PRODUCE INCORRECT ALTITUDE PLOTS

class abstract_plot():

	def __init__(self):
		self.press_to_alt = lambda press: ((press/1013.25)**(1/5.2561)-1)/(-2.25577*10**-5)
		self.alt_to_press = lambda alt: 1013.25*((100000-2.2577*alt)/(100000))**(431/82)
		self.trajectories = []
		self.markers = 'v,.1x_*HP'
		self.met_data = None

	def better_round(self, x, base=1.0):
		"""
		A simple rounding function that lets you set the base for more accurate/broader rounding
		x: Int
		base: int
		returns: float or int
		"""
		return base * np.round(x/base)

	def _input_hysplit_traj(self, file):
		"""
		Returns a dictionary of points and altitudes with time index column
		file: String of file name
		returns: dictionary of lists with keys 'lat', 'lon', 'altitude' and 'index'
		"""
		f = open(file, 'r')
		lines = f.readlines()
		f.close()
		dat_begin = 0
		for n, line in enumerate(lines):
			if 'PRESSURE' in line:
				dat_begin = 1
				data = {'lat': [], 'lon': [], 'altitude': [], 'pressure':[], 'index': []}
				continue
			elif n == 1 and self.met_data is None:
				self.met_data = list(filter(('').__ne__, line.split(' ')))[0]
			if dat_begin:
				line = line.split(' ')
				line = list(filter(('').__ne__, line))[2:]
				yr = int('20' + line[0] if int(line[0]) < 50 else '19' + line[0])
				index = dt.datetime(yr, int(line[1]), int(line[2]), int(line[3]))
				data['index'].append(index)

				data['lat'].append(float(line[7]))
				data['lon'].append(float(line[8]))
				data['altitude'].append(float(line[9]))
				data['pressure'].append(float(line[10]))

		return data

	def append(self, file):
		dat = self._input_hysplit_traj(file)
		self.trajectories.append(dat)

	def determine_extent(self, traj):
		mins = (min(traj['lat']), min(traj['lon']))
		maxs = (max(traj['lat']), max(traj['lon']))
		# Determine the max range
		lat_range = maxs[0] - mins[0]
		lon_range = maxs[1] - mins[1]
		max_range = max(lat_range, lon_range)
		factor = 0.2
		min_lat = max(-90, mins[0]-factor*max_range)
		min_lon = max(-90, mins[1]-factor*max_range)
		max_lat = min(180, maxs[0]+factor*max_range)
		max_lon = min(180, maxs[1]+factor*max_range)

		return (min_lon, max_lon, min_lat, max_lat)

	def add_latlon_lines(self, ax, lon_step=None, lat_step=None):
		if lon_step is not None:
			for i in range(-180, 180, lon_step):
				ax.axvline(x=i, ymin=-90, ymax=90, ls=':', linewidth=1, c='0.5')
		if lat_step is not None:
			for i in range(-90, 90, lat_step):
				ax.axhline(y=i, xmin=-180, xmax=180, ls=':', linewidth=1, c='0.5')


class traj_plot(abstract_plot):

	def world_plot(self, ax=None, set_extent=True, marker_spacing=1, labels=None):
		if ax == None:
			proj =  ccrs.PlateCarree()
			fig, ax = plt.subplots(figsize=(16,12),subplot_kw={'projection': proj})

			ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', 
	                                                edgecolor='black', facecolor="none"))

		if labels == None:
			labels = [None for i in range(len(self.trajectories))]

		ext = [float('inf'), -float('inf'), float('inf'), -float('inf')]
		for n, dat in enumerate(self.trajectories):
			color = ax.plot(dat['lon'], dat['lat'], transform=ccrs.PlateCarree())[0].get_color()
			ax.plot(dat['lon'], dat['lat'], 'o', markevery=marker_spacing, color=color, label=labels[n], transform=ccrs.PlateCarree())
			e = self.determine_extent(dat)
			ext[0] = min(e[0], ext[0])
			ext[1] = max(e[1], ext[1])
			ext[2] = min(e[2], ext[2])
			ext[3] = max(e[3], ext[3])

		if set_extent: ax.set_extent(ext, crs=ccrs.PlateCarree())

		ax.set_xlabel(r'Longitude($^{\circ}$E)')
		ax.set_ylabel(r'Latitude($^{\circ}$N)')

		return ax

	def altitude_plot(self, ax=None, marker_spacing=1):
		if ax == None:
			fig, ax = plt.subplots(figsize=(16,12))

		for dat in self.trajectories:
			color = ax.plot([-n for n in range(len(dat['altitude']))], dat['altitude'])[0].get_color()
			ax.plot([-n for n in range(len(dat['altitude']))], dat['altitude'], 'o', markevery=marker_spacing, color=color)

		return ax

	def add_grid(self, ax, x_range=None, y_range=None, x_step=None, y_step=None):
		if x_step is not None:
			for i in range(min(x_range), max(x_range), x_step):
				ax.axvline(x=i, ls=':', linewidth=1, c='0.5')
		if y_step is not None:
			for i in range(min(y_range), max(y_range), y_step):
				ax.axhline(y=i, ls=':', linewidth=1, c='0.5')

class freq_plot(abstract_plot):

	def _get_collisions(self, resolution, alt=False, pressure=False, alt_scale=None):
		points = [[],[]]
		recorded = []
		weights = []
		for traj in self.trajectories:
			for n in range(len(traj['lat'])):
				if alt:
					p1 = -n
					if alt_scale:
						p2 = self.better_round(alt_scale(traj['altitude'][n]), resolution)
					else:
						p2 = self.better_round(traj['altitude'][n], resolution)
				elif pressure:
					p1 = -n
					p2 = self.better_round(traj['pressure'][n], resolution)
				else:
					p1 = self.better_round(traj['lon'][n], resolution)
					p2 = self.better_round(traj['lat'][n], resolution)
				p = (p1, p2)
				if p not in recorded:
					recorded.append(p)
					points[0].append(p1)
					points[1].append(p2)
					weights.append(1)
				else:
					ind = recorded.index(p)
					weights[ind] += 1
		return np.array(points), np.array(weights)

	def _plot_kde(self, x, y, kde, resolution, x_range=None, y_range=None, resolution_y=None, pressure=False):
		if resolution_y is None:
			resolution_y = resolution
		if x_range is None:
			x_step = int((max(x)-min(x))/resolution)
			xgrid = np.linspace(min(x), max(x), x_step+1)
		else:
			x_step = int((x_range[1]-x_range[0])/resolution)
			xgrid = np.linspace(x_range[0], x_range[1], x_step+1)

		if y_range is not None:
			y_step = int((y_range[1]-y_range[0])/resolution_y)
			ygrid = np.linspace(y_range[0], y_range[1], y_step+1)
		else:
			y_step = int((max(y)-min(y))/resolution_y)
			ygrid = np.linspace(min(y), max(y), y_step+1)

			
		
		
		Xgrid, Ygrid = np.meshgrid(xgrid, ygrid)
		Z = kde.evaluate(np.vstack([Xgrid.ravel(), Ygrid.ravel()]))
		Zgrid = Z.reshape(Xgrid.shape)

		ext = ([min(x), max(x)] if x_range is None else x_range) + ([min(y), max(y)] if y_range is None else y_range)

		return Xgrid, Ygrid, Zgrid, Z, ext





	def contour_plot(self, ax=None, resolution=0.5, manual_bandwidth=None, scale=None, x_range=None, y_range=None, imshow_kwargs=None):
		if ax == None:
			proj =  ccrs.PlateCarree()
			fig, ax = plt.subplots(figsize=(16,12),subplot_kw={'projection': proj})

			ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', 
	                                                edgecolor='black', facecolor="none"))
		if imshow_kwargs == None:
			imshow_kwargs = {}
		collision_points, collision_weights = self._get_collisions(resolution)

		# print("Cont neff:", int(sum(collision_weights))**2 / sum([int(w)**2 for w in collision_weights]))

		if scale is not None:
			collision_weights = scale(collision_weights)

		if manual_bandwidth is None:
			kde = gaussian_kde(collision_points, weights=collision_weights)
		else:
			kde = gaussian_kde(collision_points, weights=collision_weights, bw_method = manual_bandwidth)

		x, y = collision_points

		vals = self._plot_kde(x, y, kde, resolution, x_range, y_range)

		Zgrid = vals[2]
		ext = vals[4]
		earth = plt.cm.gist_earth_r
		# print(ext)

		im = ax.imshow(Zgrid,
				origin='lower', aspect='auto',
				extent=ext,
				alpha=0.8,
				cmap=earth,
				transform=ccrs.PlateCarree(), **imshow_kwargs)

		ax.set_xlabel(r'Longitude($^{\circ}$E)')
		ax.set_ylabel(r'Latitude($^{\circ}$N)')

		ax.set_extent(ext, ccrs.PlateCarree())
		ax.set_aspect('equal')

		return ax, im

	def altitude_contour_plot(self, ax=None, resolution=100, manual_bandwidth=None, scale=None, x_range=None, y_range=None, 
		imshow_kwargs=None, yscale = None, yscale_inv = None, round_yaxis=None, yticks=None):
		if ax == None:
			fig, ax = plt.subplots(figsize=(16,12))

		if imshow_kwargs == None:
			imshow_kwargs = {}

		collision_points, collision_weights = self._get_collisions(resolution, alt=True, alt_scale=yscale)
		if scale is not None:
			collision_weights = scale(collision_weights)

		print("Alt neff:", int(sum(collision_weights))**2 / sum([int(w)**2 for w in collision_weights]))

		if manual_bandwidth is None:
			kde = gaussian_kde(collision_points, weights=collision_weights)
		else:
			kde = gaussian_kde(collision_points, weights=collision_weights, bw_method = manual_bandwidth)

		x, y = collision_points

		vals = self._plot_kde(x, y, kde, 1, x_range, y_range, resolution)

		Zgrid = vals[2]
		ext = vals[4]
		earth = plt.cm.gist_earth_r

		im = ax.imshow(Zgrid,
				origin='lower', aspect='auto',
				extent=ext,
				alpha=0.8,
				cmap=earth, **imshow_kwargs)

		ax.set_xlabel(r'Time (h)')
		ax.set_ylabel(r'Altitude (m AGL)')


		if yscale: 
			ax.set_yscale('function', **{'functions':[yscale, yscale_inv]})
			if yticks == None:
				alt_ticks = self.better_round(yscale_inv(ax.get_yticks()[1:-1]), 100)
				wanted_ticks = list(yscale(alt_ticks))
			else:
				wanted_ticks = list(map(yscale, yticks))
			print(wanted_ticks)
			ax.set_yticks(wanted_ticks)


			if round_yaxis:
				ax.set_yticklabels(self.better_round(yscale_inv(ax.get_yticks()), round_yaxis))
			else:
				ax.set_yticklabels(yscale_inv(ax.get_yticks()))

		return ax, im

	def pressure_contour_plot(self, ax=None, resolution=10, manual_bandwidth=None, scale=None, x_range=None, y_range=None, 
		imshow_kwargs=None, alt_yaxis=False, round_yaxis=False, yticks=None):

		if ax == None:
			fig, ax = plt.subplots(figsize=(16,12))

		if imshow_kwargs == None:
			imshow_kwargs = {}

		collision_points, collision_weights = self._get_collisions(resolution, pressure=True)
		if scale is not None:
			collision_weights = scale(collision_weights)

		print("Pressure neff:", int(sum(collision_weights))**2 / sum([int(w)**2 for w in collision_weights]))

		if manual_bandwidth is None:
			kde = gaussian_kde(collision_points, weights=collision_weights)
		else:
			kde = gaussian_kde(collision_points, weights=collision_weights, bw_method = manual_bandwidth)

		x, y = collision_points

		vals = self._plot_kde(x, y, kde, 1, x_range, y_range, resolution)

		Zgrid = vals[2]
		ext = vals[4]
		earth = plt.cm.gist_earth_r

		im = ax.imshow(Zgrid,
				origin='lower', aspect='auto',
				extent=ext,
				alpha=0.8,
				cmap=earth, **imshow_kwargs)

		ax.set_xlabel(r'Time (h)')
		ax.set_ylabel(r'Altitude (m ASL)' if alt_yaxis else r'Pressure (hPa)')
		if alt_yaxis: 
			ax.set_yscale('function', **{'functions':[self.press_to_alt, self.alt_to_press]})
			if yticks == None:
				alt_ticks = self.better_round(self.press_to_alt(ax.get_yticks()[1:-1]), 100)
				wanted_ticks = list(self.alt_to_press(alt_ticks))
			else:
				wanted_ticks = list(map(self.alt_to_press, yticks))

			ax.set_yticks(wanted_ticks)


			if round_yaxis:
				ax.set_yticklabels(self.better_round(self.press_to_alt(ax.get_yticks()), round_yaxis))
			else:
				ax.set_yticklabels(self.press_to_alt(ax.get_yticks()))
		return ax, im

		

	def kml_contour_plot(self, fname='freq_plot.kml', resolution=0.5, log_density=False):
		collision_points, collision_weights = self._get_collisions(resolution)
		if log_density:
			collision_weights = np.log(collision_weights)

		kde = gaussian_kde(collision_points, weights=collision_weights)

		x, y = collision_points

		Xgrid, Ygrid, Zgrid, Z = self._plot_kde(x, y, kde, resolution)

		ext = [min(x), max(x), min(y), max(y)]
		earth = plt.cm.gist_earth_r

		cmap = earth
		my_cmap = cmap(np.arange(cmap.N))
		my_cmap[:,-1] = np.array([1-exp(-pi*(i/256)) for i in range(256)]) # This is just a way of tapering off the alpha, so insignificant values are less visible
		my_cmap = ListedColormap(my_cmap)

		kml = simplekml.Kml()

		count_min = float(abs(min(Z)))
		count_max = float(abs(max(Z)))


		for lats, lons, counts in zip(Xgrid, Ygrid, Zgrid):
			for i in range(len(lats)):
				lat = lats[i]
				lon = lons[i]
				count = counts[i]

				pol = kml.newpolygon()
				pol.outerboundaryis = [(lat, lon), (lat+resolution, lon), (lat+resolution, lon+resolution), (lat, lon+resolution), (lat,lon)]

				normalised_count = (count-count_min)/(count_max-count_min)	
				rgb_color = [float(i) for i in my_cmap(normalised_count)]

				color = simplekml.Color.rgb(round(rgb_color[0]*255), round(rgb_color[1]*255), round(rgb_color[2]*255), round(rgb_color[3]*255))

				pol.style.linestyle.width = 0
				pol.style.polystyle.color = color
		kml.save(fname)





if __name__ == '__main__':
	base = 'N:/coala/'
	files = os.listdir(base)
	# tp = traj_plot()
	# for i in files[:10]:
	# 	tp.append(base + i)

	# fig = plt.figure(figsize=(12, 12))
	# grid = plt.GridSpec(8, 1, hspace=0.1, wspace=0.2)
	# main_ax = fig.add_subplot(grid[:4,0], projection=ccrs.PlateCarree())
	# main_ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', 
	#                                                 edgecolor='black', facecolor="none"))
	# main_ax.axis('on')
	# main_ax.get_xaxis().set_visible(True)
	# main_ax.get_yaxis().set_visible(True)
	# main_ax.set_title("Trajectory plot")
	# alt_ax = fig.add_subplot(grid[5:,0])
	# alt_ax.set_title("Altitude plot")

	# tp.add_latlon_lines(main_ax, 20, 10)
	# tp.world_plot(main_ax, set_extent=True)
	# tp.altitude_plot(alt_ax)


	#plt.show()


	fp = freq_plot()
	for i in files[:10]:
		if "." not in i:
			fp.append(base + i)

	#projec = ccrs.SouthPolarStereo()
	#projec = ccrs.LambertConformal(central_longitude=146.316, central_latitude=-60, false_easting=0, false_northing=-90)
	#projec = ccrs.Robinson(130)
	projec = ccrs.NearsidePerspective(central_longitude=146.316, central_latitude=-60)
	fig = plt.figure(figsize=(12, 12))
	grid = plt.GridSpec(8, 1, hspace=0.1, wspace=0.2)
	main_ax = fig.add_subplot(grid[:4,0], projection=projec)
	main_ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', 
	                                                edgecolor='black', facecolor="none"))
	main_ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

	# ocean110 = cfeature.NaturalEarthFeature('physical', 'ocean', \
 #        scale='110m', edgecolor='none', facecolor=cfeature.COLORS['water'])
	# main_ax.add_feature(ocean110, zorder=-5)

	# land110 = cfeature.NaturalEarthFeature('physical', 'land', '110m', \
	#         edgecolor='black', facecolor="silver")
	#main_ax.add_feature(land110, zorder=5)
	main_ax.axis('on')
	#main_ax.get_xaxis().set_visible(True)
	#main_ax.get_yaxis().set_visible(True)
	main_ax.set_title("Contour plot")
	alt_ax = fig.add_subplot(grid[5:,0])
	alt_ax.set_title("Altitude contour plot")

	#fp.add_latlon_lines(main_ax, 20, 10)
	fp.contour_plot(main_ax)
	#main_ax.invert_yaxis()
	#main_ax.invert_xaxis()
	fp.altitude_contour_plot(alt_ax, resolution=10)#, yscale=np.log1p, yscale_inv = lambda x: np.exp(x) - 1)
	alt_ax.set_yscale('symlog', linthreshy=0.5)

	plt.show()
