'''
Mixed Solvent MD tools
'''
from .xvg import XVG
import numpy as np

class MixSolvMD:
    '''
    MixSolvMD class to read and analyse mixed solvent MD run.
    It requires ligands trajectories in macromolecule and in bulk systems.

    macromol_file: name of the XVG file for ligand(s) trajectories in macromolecule system
    bulk_file: name of the XVG file for ligand(s) trajectories in bulk system
    system_name: name of the system for PDB title
    box: size of the simulation box (in nm)
    membrane: boolean (default False)
    temperature: temperature in which the trajectories have been sampled (in Kelvin, default 310.15)
    '''
    def __init__(
        self,
        macromol_file: str,
        bulk_file: str,
        system_name:str | None = None,
        box: list = [1.0,1.0,1.0],
        membrane: bool = False,
        temperature: float = 310.15
    ):
        self.hotspots = []
        self.system_name = system_name
        self.box = box
        self.avogadro=6.02214076e23         # Avogadro constant in mol-1
        self.R=1.985e-3                     # Gas constant in kcal.mol-1
        self.temperature = temperature
        self.macromol_traj=XVG(macromol_file).get3Dcoord()
        self.bulk_traj=XVG(bulk_file).get3Dcoord()
        if self.macromol_traj and self.bulk_traj:
            if len(self.macromol_traj) != len(self.bulk_traj):
                print("WARNING! Macromolecule and bulk trajectories don't have the same length! They probably don't have the same sampling time.")
            assert len(self.macromol_traj[0][0]) == 3 and len(self.bulk_traj[0][0]) == 3, 'ERROR! Coordinates provided are not 3D coordinates!'
        else:
            raise ValueError('ERROR! No macromolecule and/or bulk trajectories provided!')
        if membrane:
            memb_height = 4
            self.macromol_concentration = (len(self.macromol_traj[0])/(self.avogadro * (self.box[0]*self.box[1]*(self.box[2]-memb_height))*1e-24))
            self.bulk_concentration = (len(self.bulk_traj[0])/(self.avogadro * (self.box[0]*self.box[1]*(self.box[2]-memb_height))*1e-24))
        else:
            self.macromol_concentration = (len(self.macromol_traj[0])/(self.avogadro * (self.box[0]*self.box[1]*self.box[2])*1e-24))
            self.bulk_concentration = (len(self.bulk_traj[0])/(self.avogadro * (self.box[0]*self.box[1]*self.box[2])*1e-24))
        if self.macromol_concentration != self.bulk_concentration:
            print("WARNING! Macromolecule and bulk systems don't have the same ligand concentrations")

    def get_hotspots(
        self,
        step: float = 0.2,
        limits: list = []
    ):
        '''
        Function to generate the hotspots from the trajectories.

        step: voxel size (in nm, default = 0.2)
        limits: list of x, y, and z limits to reduce the hotspots to a subregion of the system ([(x_low_limit, x_high_limit), (y_low_limit, y_high_limit), (z_low_limit, z_high_limit)])
        '''
        def cell_center(cell: tuple):
            # cell is a tuple of (xmin, xmax, ymin, ymax, zmin, zmax) boundaries
            xmin, xmax, ymin, ymax, zmin, zmax = cell
            
            # get the coordinate of the center of the cell
            center = ((xmin+xmax)/2 , (ymin+ymax)/2 , (zmin+zmax)/2)

            return center

        # Flatten coordinates from trajectories of all ligands
        macromol_coord = np.asarray([coord for lig_traj in self.macromol_traj for coord in lig_traj])
        bulk_coord = np.asarray([coord for lig_traj in self.bulk_traj for coord in lig_traj])

        # Run the count in grid with np.histogramdd
        if limits:
            H_macro, edges_macro = np.histogramdd(macromol_coord, bins=(int(self.box[0]/step),int(self.box[1]/step),int(self.box[2]/step)), range=((limits[0][0], limits[0][1]), (limits[1][0], limits[1][1]), (limits[2][0], limits[2][1])))
            H_bulk, edges_bulk = np.histogramdd(bulk_coord, bins=(int(self.box[0]/step),int(self.box[1]/step),int(self.box[2]/step)), range=((limits[0][0], limits[0][1]), (limits[1][0], limits[1][1]), (limits[2][0], limits[2][1])))
        else:
            H_macro, edges_macro = np.histogramdd(macromol_coord, bins=(int(self.box[0]/step),int(self.box[1]/step),int(self.box[2]/step)), range=((0,self.box[0]), (0, self.box[1]), (0, self.box[2])))
            H_bulk, edges_bulk = np.histogramdd(bulk_coord, bins=(int(self.box[0]/step),int(self.box[1]/step),int(self.box[2]/step)), range=((0,self.box[0]), (0, self.box[1]), (0, self.box[2])))
        # Calculate occupancies and free energies
        self.macromol_occup = H_macro / len(self.macromol_traj)
        self.bulk_occup = H_bulk / len(self.bulk_traj)
        self.free_energy = -self.R*self.temperature*np.log(self.macromol_occup/self.bulk_occup)
        # Generate hotspots
        for i, ox in enumerate(self.macromol_occup):
            x_edges = (edges_macro[0][i], edges_macro[0][i+1])
            for j, oy in enumerate(ox):
                y_edges = (edges_macro[1][j], edges_macro[1][j+1])
                for k, oz in enumerate(oy):
                    z_edges = (edges_macro[2][k], edges_macro[2][k+1])
                    cell = (x_edges[0], x_edges[1], y_edges[0], y_edges[1], z_edges[0], z_edges[1])
                    if np.isnan(self.free_energy[i][j][k]) or np.isinf(self.free_energy[i][j][k]) or np.isneginf(self.free_energy[i][j][k]):
                        self.hotspots.append([cell_center(cell), oz, np.nan])
                    else:
                        self.hotspots.append([cell_center(cell), oz, self.free_energy[i][j][k]])

    def get_highest_occupancy(
        self,
        n: int = 5
    ):
        '''
        Returns N hotspots with the highest occupancy ranked from highest to lowest
        '''
        if not self.get_hotspots:
            raise ValueError('ERROR! Hotspots not generated')
        sorted_occup = sorted(enumerate([(hs[1], hs[2]) for hs in self.hotspots], start=1), key=lambda x: x[1], reverse=True)
        top = sorted_occup[:n]

        return top

    def get_lowest_free_energy(
        self,
        n: int = 5
    ):
        '''
        Returns N hotspots with the highest free energy ranked from highest to lowest
        '''
        if not self.get_hotspots:
            raise ValueError('ERROR! Hotspots not generated')
        # Need to convert nan to 0 for the sorted function
        sorted_fe = sorted(enumerate([(hs[2], hs[1]) if not np.isnan(hs[2]) else (0, hs[1]) for hs in self.hotspots], start=1), key=lambda x: x[1], reverse=False)
        top = sorted_fe[:n]

        return top

    def write_pdb(
        self,
        output: str = 'file.pdb',
        title: str | None = None
    ):
        '''
        Write the coordinates of the hotspots with their corresponding occupancy and free energy in PDB format.

        output: name of the PDB file
        title: title of the PDB file (default: self.system_name or "Mixed Solvent MD hotspots")
        '''
        from datetime import date
        if not self.get_hotspots:
            raise ValueError('ERROR! Hotspots not generated')
        if not title:
            if not self.system_name:
                title='Mixed Solvent MD hotspots'
            else:
                title=self.system_name
        with open(output, 'w') as f:
            f.writelines(f'HEADER    MIXED-SOLVENT MD    OCCUP=PROT_OCCUP; BETA=FREE_ENERGY_(KCAL/MOL)    {date.today()}\n')
            f.writelines(f'TITLE     {title}\n')
            f.writelines('')
            for i,d in enumerate(self.hotspots):
                if i < 99999:
                    j = i+1
                    f.writelines(f"ATOM  {'%5s' % j}  D   DUM A   1     {'%7.3f' % (d[0][0]*10)} {'%7.3f' % (d[0][1]*10)} {'%7.3f' % (d[0][2]*10)}  {'%.2f' % d[1]}  {'%.2f' % d[2]}           D\n")
                else:
                    if str(i)[-5:] == '99999':
                        j = 0
                    else:
                        j += 1
                    f.writelines(f"ATOM  {'%5s' % j}  D   DUM A   1     {'%7.3f' % (d[0][0]*10)} {'%7.3f' % (d[0][1]*10)} {'%7.3f' % (d[0][2]*10)}  {'%.2f' % d[1]}  {'%.2f' % d[2]}           D\n")

    def dist_to_probe(
        self,
        site_file: str,
        unit: str = 'nm'
    ):
        '''
        Retunrs the distances between the probe/ligand and a specific site of interest

        site_file: name of the XVG file for site trajectories in macromolecule system
        unit: metrics unit in which the distances will be returned (either "nm" or "Ang", default "nm")
        '''
        site_traj=XVG(site_file).get3Dcoord()
        distances = []
        
        for site, probe in zip(site_traj, self.macromol_traj):
            d = []
            for c in probe:
                dist = np.linalg.norm(np.array(site)-np.array(c))
                if unit=='nm':
                    d.append(dist)
                elif unit=='Ang':
                    d.append(dist*10)
                else:
                    raise ValueError('ERROR! Unit not specified or not valid')
            distances.append(min(d))
            
        return distances

    
