'''
Mixed Solvent MD tools
'''
import numpy as np

class MixSolvMD:
    '''
    MixSolvMD class to read and analyse mixed solvent MD runs (requires ligand trajectories in macromolecule and in bulk systems)
    '''
    def __init__(self, macromol_traj, bulk_traj, system_name='', box = [1.0,1.0,1.0], membrane=False, unit = 'Ang', temperature=310.15):
        self.grid = []
        self.hotspots = []
        self.system_name = system_name
        self.box = box
        self.unit = unit
        if self.unit == 'nm':
            vol_conv = 1e-24
            memb_height = 4
        elif self.unit == 'Ang':
            vol_conv = 1e-27
            memb_height = 40
        else:
            raise ValueError('ERROR! Metric unit is not valid!')
        self.avogadro=6.02214076e23         # Avogadro constant in mol-1
        self.R=1.985e-3                     # Gas constant in kcal.mol-1
        self.temperature = temperature
        if macromol_traj and bulk_traj:
            if len(macromol_traj) != len(bulk_traj):
                print("WARNING! Macromolecule and bulk trajectories don't have the same length! They probably don't have the same sampling time.")
            if len(macromol_traj[0][0]) == 3 and len(bulk_traj[0][0]) == 3:
                self.macromol_traj = macromol_traj
                self.bulk_traj = bulk_traj
            else:
                raise ValueError('ERROR! Coordinates provided are not 3D coordinates!')
        else:
            raise ValueError('ERROR! No macromolecule and/or bulk trajectories provided!')
        if membrane:
            self.macromol_concentration = (len(macromol_traj[0])/(self.avogadro * (self.box[0]*self.box[1]*(self.box[2]-memb_height))*vol_conv))
            self.bulk_concentration = (len(bulk_traj[0])/(self.avogadro * (self.box[0]*self.box[1]*(self.box[2]-memb_height))*vol_conv))
        else:
            self.macromol_concentration = (len(macromol_traj[0])/(self.avogadro * (self.box[0]*self.box[1]*self.box[2])*vol_conv))
            self.bulk_concentration = (len(bulk_traj[0])/(self.avogadro * (self.box[0]*self.box[1]*self.box[2])*vol_conv))
        if self.macromol_concentration != self.bulk_concentration:
            print("WARNING! Macromolecule and bulk systems don't have the same ligand concentrations")
        
    def gen_grid(self, x_limit=[], y_limit=[], z_limit=[], step=0.5):
        all_x = []
        all_y = []
        all_z = []
        for macromol_coords, bulk_coords in zip(self.macromol_traj, self.bulk_traj):
            for macromol_c,bulk_c in zip(macromol_coords,bulk_coords):
                all_x.append(macromol_c[0])
                all_x.append(bulk_c[0])
                all_y.append(macromol_c[1])
                all_y.append(bulk_c[1])
                all_z.append(macromol_c[2])
                all_z.append(bulk_c[2])
        xmin = min(all_x)
        xmax = max(all_x)
        if not x_limit:
            x_limit=[0, self.box[0]]
        if xmin < (x_limit[0]):
            xmin = x_limit[0]
        if xmax > (x_limit[-1]):
            xmax = x_limit[-1]
        ymin = min(all_y)
        ymax = max(all_y)
        if not y_limit:
            y_limit=[0, self.box[1]]
        if ymin < (y_limit[0]):
            ymin = y_limit[0]
        if ymax > (y_limit[-1]):
            ymax = y_limit[-1]
        zmin = min(all_z)
        zmax = max(all_z)
        if not z_limit:
            z_limit = [0, self.box[2]]
        if zmin < (z_limit[0]):
            zmin = z_limit[0]
        if zmax > (z_limit[-1]):
            zmax = z_limit[-1]
        
        xv = np.arange(xmin, xmax, step)
        yv = np.arange(ymin, ymax, step)
        zv = np.arange(zmin, zmax, step)
        
        for x in xv:
            for y in yv:
                for z in zv:
                    self.grid.append((x, x+step, y, y+step, z, z+step))
        
    def get_hotspots(self, step=0.5, progress=True):
        # Generate grid
        if not self.grid:
            self.gen_grid(step=step)
        
        # Define functions
        def is_inside(coords, cell):
            # coords is a list of 3D coordinates
            # cell is a tuple of (xmin, xmax, ymin, ymax, zmin, zmax) boundaries from self.grid
            xmin, xmax, ymin, ymax, zmin, zmax = cell

            check = [(xmin <= c[0] < xmax and ymin <= c[1] < ymax and zmin <= c[2] < zmax) for c in coords]

            # return a list of boolean for the coordinates inside the cell or not
            return check
        
        def cell_center(cell):
            # cell is a tuple of (xmin, xmax, ymin, ymax, zmin, zmax) boundaries from self.grid
            xmin, xmax, ymin, ymax, zmin, zmax = cell
            
            # get the coordinate of the center of the cell
            center = ((xmin+xmax)/2 , (ymin+ymax)/2 , (zmin+zmax)/2)
            
            return center
        
        # Run the count in grid
        if progress:
            from tqdm import tqdm
            iteration = tqdm(self.grid)
        else:
            iteration = self.grid
            # loop through the grid cells
        for cell in iteration:
            # initialize the count for this cell to zero
            macromol_counts = 0
            bulk_counts = 0
            # loop through the coordinates with protein
            for m_coord, b_coord in zip(self.macromol_traj, self.bulk_traj):
                # check if the macromolecule coordinates are inside the cell
                check = is_inside(m_coord, cell)
                # increment the count for this cell by one
                macromol_counts += sum(check)
                # repeat for coordinates in bulk
                check = is_inside(b_coord, cell)
                # increment the count for this cell by one
                bulk_counts += sum(check)
                    
            # append coordinates and counts if count_c1 - count_c2 > 0
            diff = macromol_counts - bulk_counts
            macromol_occup = macromol_counts / len(self.macromol_traj)
            bulk_occup = bulk_counts / len(self.bulk_traj)
            if diff > 0:
                if bulk_occup == 0:
                    self.hotspots.append([cell_center(cell), macromol_occup, np.nan])
                else:
                    free_energy = -self.R*self.temperature*np.log(macromol_occup/bulk_occup)
                    self.hotspots.append([cell_center(cell), macromol_occup, free_energy])

    def get_hotspots_v2(self, step=0.5, limits=[]):

        def cell_center(cell):
            # cell is a tuple of (xmin, xmax, ymin, ymax, zmin, zmax) boundaries from self.grid
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

    def get_highest_occupancy(self, n=5):
        if not self.get_hotspots:
            raise ValueError('ERROR! Hotspots not generated')
        sorted_occup = sorted(enumerate([(hs[1], hs[2]) for hs in self.hotspots], start=1), key=lambda x: x[1], reverse=True)
        top = sorted_occup[:n]

        return top

    def get_lowest_free_energy(self, n=5):
        if not self.get_hotspots:
            raise ValueError('ERROR! Hotspots not generated')
        # Need to convert nan to 0 for the sorted function
        sorted_fe = sorted(enumerate([(hs[2], hs[1]) if not np.isnan(hs[2]) else (0, hs[1]) for hs in self.hotspots], start=1), key=lambda x: x[1], reverse=False)
        top = sorted_fe[:n]

        return top

    def write_pdb(self, output='file.pdb', title=''):
        from datetime import date
        if not self.get_hotspots:
            raise ValueError('ERROR! Hotspots not generated')
        if not title:
            if not self.system_name:
                title='Mixed Solvent MD hotspots'
            else:
                title=self.system_name
        if self.unit == 'nm':
            convert=10          # Convert nm to Ang
        elif self.unit == 'Ang':
            convert=1           # Don't convert
        with open(output, 'w') as f:
            f.writelines(f'HEADER    MIXED-SOLVENT MD    OCCUP=PROT_OCCUP; BETA=FREE_ENERGY_(KCAL/MOL)    {date.today()}\n')
            f.writelines(f'TITLE     {title}\n')
            f.writelines('')
            for i,d in enumerate(self.hotspots):
                if i < 99999:
                    j = i+1
                    f.writelines(f"ATOM  {'%5s' % j}  D   DUM A   1     {'%7.3f' % (d[0][0]*convert)} {'%7.3f' % (d[0][1]*convert)} {'%7.3f' % (d[0][2]*convert)}  {'%.2f' % d[1]}  {'%.2f' % d[2]}           D\n")
                else:
                    if str(i)[-5:] == '99999':
                        j = 0
                    else:
                        j += 1
                    f.writelines(f"ATOM  {'%5s' % j}  D   DUM A   1     {'%7.3f' % (d[0][0]*convert)} {'%7.3f' % (d[0][1]*convert)} {'%7.3f' % (d[0][2]*convert)}  {'%.2f' % d[1]}  {'%.2f' % d[2]}           D\n")

    def dist_to_probe(self, site_traj, unit='nm'):
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

    
