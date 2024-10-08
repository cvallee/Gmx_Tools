'''
Binding events analysis
'''
from .xvg import XVG
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def binding_events(
	ligand_file: str,
	site_file: str,
	min_dist: float = 0.6,
	unit: str = 'nm',
	time: str = 'ps',
	plot: bool = False
) -> dict:
	'''Analyse binding events between ligand and specific site from GROMACS simulation.
	It requires the trajectories of the CoM of the ligand(s) and the CoM of the specific site and return a dictionary of the number of binding events with their retention time per ligand.
	
	ligand_file: name of the XVG file for ligand(s) CoM trajectories
	site_file: name of the XVG file for site CoM trajectories
	min_dist: minimum distance to consider a binding events (default 0.6 nm or 6 Ang)
	unit: unit used for metrics (either "nm" or "Ang", default "nm")
	time: unit used for time (either "ps" or "ns", default "ps")
	plot: boolean (default False)
	'''
	data={
	'Binding_events':[],
	'Retention_times':[]
	}
	# Get the trajectories from files
	lig_xvg=XVG(ligand_file)
	lig_traj=lig_xvg.get3Dcoord()
	site_xvg=XVG(site_file)
	site_traj=site_xvg.get3Dcoord()
	# Assert the trajectories have the same length
	assert len(lig_traj) == len(site_traj)
	# Loop over all the ligands if there are multiple ligands in the ligand XVG file
	for i in range(len(lig_traj[0])):
		# Get distances between ligand(s) and center of mass for each timestep
		dist=[np.linalg.norm(np.array(site)-np.array(lig_traj[j][i])) for j, site in enumerate(site_traj)]
		if unit=='Ang':
			# Convert to Ang
			dist = [d*10 for d in dist]
		elif unit!='nm':
			raise ValueError('ERROR! Metrics unit not supported! Please use either nm or Ang.')
		# Get binding event and retention time for each ligands
		is_bound=False
		nb_binding_events=0
		retention_times=[]
		dt = lig_xvg.x_column[1]-lig_xvg.x_column[0]
		if time == 'ns':
			dt *= 0.001
		elif time != 'ps':
			raise ValueError('ERROR! Time unit not supported! Please use either ps or ns.')
		for d in dist:
			if d <= min_dist:
				if not is_bound:
					is_bound=True
					nb_binding_events += 1
					ret_time = dt
				else:
					ret_time += dt
			else:
				if is_bound:
					is_bound=False
					if ret_time:
						retention_times.append(ret_time)
		if is_bound:
			if ret_time:
				retention_times.append(ret_time)

		if nb_binding_events:
			data['Binding_events'].append(nb_binding_events)
			data['Retention_times'].append(retention_times)

		if plot:
			if time == 'ps':
				plt.plot([t for t in lig_xvg.x_column], dist)
			elif time == 'ns':
				plt.plot([t*0.001 for t in lig_xvg.x_column], dist)

	if plot:
		plt.title(f'Binding events\n({lig_xvg.filename} & {site_xvg.filename})')
		plt.ylim(0,min_dist*2)
		plt.xlabel(f'Time ({time})')
		plt.ylabel(f'Distance Ligand-Site ({unit})')
		plt.axhline(min_dist, linestyle=':', color='red', label='Binding threshold')
		plt.legend()
		plt.show()

	if not data['Retention_times']:
		data={'Binding_events':0,'Retention_times':None}
	
	return data