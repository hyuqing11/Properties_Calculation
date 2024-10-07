
import numpy as np
from Read_String_Data import Read_String_Data
from String_Visualization import String_Visualization
from ComputeDynamicProperties import ComputeDynamicProperties
from ComputeDynamicFacilitation import ComputeDynamicFacilitation
from WriteFile import WriteFile

from ComputePropensity import ComputePropensity
import logging


class Simulator:
    def __init__(self, parameters, pos, vel, latt, atom_pos_dict, atom_vel_dict, latt_matrix):
        self.parameters = parameters
        self.pos = pos
        self.vel = vel
        self.latt = latt
        self.atom_pos_dict = atom_pos_dict
        self.atom_vel_dict = atom_vel_dict
        self.latt_matrix = latt_matrix
        self.logger = logging.getLogger(__name__)

        # Mapping property types to their corresponding methods
        self.property_methods = {
            0: self._compute_vacf_and_pdos,
            1: self._compute_van_hove,
            2: self._compute_dynamic_structure,
            4: self._compute_propensity,
            6: self._visualize_string,
            7: self._compute_dynamic_facilitation
        }


    def simulator_choice(self, folder):
        """Select and execute the appropriate simulation method based on the property type."""
        property_type = self.parameters.get('property_type')

        if property_type in self.property_methods:
            self.property_methods[property_type](folder)
        else:
            raise ValueError(f"Unknown property type: {property_type}")


    def _visualize_string(self, folder):
        logging.info("Visualize")
        string_reader = Read_String_Data(folder, self.parameters)
        string_libraries = string_reader.read_array_string()
        visual = String_Visualization(self.pos, self.parameters, folder, self.latt_matrix, string_libraries)
        visual.string_visual()
        print("Visualization completed")

    def _compute_van_hove(self,folder):
        """Compute the van Hove self-correlation function and write results."""
        logging.info("Compute van Hove self-correlation")
        for j in self.parameters['compute_type']:
            if j > self.parameters['num_types'] + 1:
                raise ValueError(f"Unknown atom type: {j}")
            if j == self.parameters['num_types'] + 1:
                van_hov_cal = ComputeDynamicProperties(self.pos, self.vel, self.parameters,
                                                       self.latt)
            else:
                van_hov_cal = ComputeDynamicProperties(self.atom_pos_dict[j], self.atom_vel_dict[j], self.parameters,
                                                       self.latt)
            Gr_mean, shells, r = van_hov_cal.calculate_van_hove_function()
            self._write_van_hove_results(folder, j, Gr_mean, r)

    def _compute_vacf_and_pdos(self,folder):
        """Compute the Velocity Auto-Correlation Function (VACF) and Phonon Density of States (PDOS)."""

        logging.info('Compute VACF and Pdos')
        omega = np.arange(0, self.parameters['max_omega'], self.parameters['d_omega'])
        nu = omega / (2 * np.pi)
        dt = self.parameters['dt']
        t = np.arange(self.parameters['Nc']) * dt

        for j in self.parameters['compute_type']:
            if j == self.parameters['num_types'] + 1:
                pdos_cal = ComputeDynamicProperties(self.pos, self.vel, self.parameters, self.latt)
                vacf_non, vacf_output, pdos = pdos_cal.pdos(omega)
            else:
                pdos_cal = ComputeDynamicProperties(self.atom_pos_dict[j], self.atom_vel_dict[j], self.parameters,
                                                    self.latt)
                vacf_non, vacf_output, pdos = pdos_cal.pdos(omega)
            self._write_vacf_and_pdos(j, nu, pdos, t, vacf_non, vacf_output,folder)

    def _compute_dynamic_structure(self,folder):
        logging.info("Compute dynamics structure")
        wr1, wr2, wr3 = self.parameters['write_parameters']
        omega = np.arange(0, self.parameters['max_omega'], self.parameters['d_omega'])
        nu = omega / (2 * np.pi)
        dt = self.parameters['dt']
        t = np.arange(self.parameters['Nc']) * dt
        for j in self.parameters['compute_type']:
            if j > self.parameters['num_types'] + 1:
                raise ValueError(f"Unknown atom type: {j}")
            if j == self.parameters['num_types'] + 1:
                dynamic_cal = ComputeDynamicProperties(self.pos, self.vel, self.parameters,
                                                       self.latt)
            else:
                dynamic_cal = ComputeDynamicProperties(self.atom_pos_dict[j], self.atom_vel_dict[j], self.parameters,self.latt)
            fd_scale = dynamic_cal.calculate_intermediate_scattering()
            if wr1:
                self._write_intermediate_scattering(folder, j, t, fd_scale)

            Sv = dynamic_cal.calculate_dynamic_structure(omega, fd_scale)
            if wr2:
                self._write_dynamic_structure(folder, j, nu, Sv)
            S_intgr = dynamic_cal.Integrate_dynamic_structure(Sv)
            if wr3:
                self._write_integration_dynamic_structure(folder, j, nu, S_intgr)

    def _compute_dynamic_facilitation(self,folder):
        logging.info("Compute dynamic facilitation")

        for j in self.parameters['compute_type']:
            if j > self.parameters['num_types'] + 1:
                raise ValueError(f"Unknown atom type: {j}")
            if j == self.parameters['num_types'] + 1:
                facilitation_cal = ComputeDynamicFacilitation(self.pos, self.vel, self.parameters,
                                                       self.latt)
            else:
                facilitation_cal = ComputeDynamicFacilitation(self.atom_pos_dict[j], self.atom_vel_dict[j], self.parameters,
                                                       self.latt)

        r, mean_Shells_record, mean_Shells1_record = facilitation_cal.compute_probability_distribution()


        self._write_dynamic_facilitation(folder,r,mean_Shells_record,mean_Shells1_record)

    def _compute_propensity(self, folder):
        """Compute spatial correlation of propensity and write results."""
        logging.info("Compute spatial correlation of propensity")
        for j in self.parameters['compute_type']:
            if j >= self.parameters['num_types'] + 1:
                raise ValueError(f"Unknown atom type: {j}")
            else:
                propensity_cal = ComputePropensity(folder, self.parameters)

        Cd = propensity_cal.compute_spatial_corr()
        self._write_Cd_results(folder, Cd)

    def _write_Cd_results(self,folder:str,Cd:list):
        sz_cd = np.shape(Cd)[0]
        time =self.parameters['dt'] * np.array(range(sz_cd))
        Cd_write = WriteFile(folder,self.parameters)
        Cd_write.write_cd(Cd, folder,time)

    ### Writing Methods (Reusable) ###

    def _write_results(self, folder, filename, x, y):
        """Helper function for writing results."""
        writer = WriteFile(folder, self.parameters)
        writer.write_results(filename, x, y)

    def _write_van_hove_results(self, folder: str, atom_type: int, Gr_mean: np.ndarray, r: np.ndarray):
        sz_time = np.size(self.parameters['time_series'])
        for i in range(sz_time):
            output_filename = f'atom_Gs_{atom_type}_{self.parameters["time_series"][i]}.txt'
            self._write_results(folder, output_filename, r, Gr_mean[i])

    def _write_intermediate_scattering(self, folder, atom_type, t, fd_scale):
        """Write results of intermediate scattering."""
        for i in range(self.parameters['vectors']):
            output_filename = f'Intermediate_scattering_{atom_type}_{i}.txt'
            self._write_results(folder, output_filename, t, fd_scale[i])

    def _write_dynamic_structure(self, folder, atom_type, nu, Sv):
        """Write results for dynamic structure."""
        for i in range(self.parameters['vectors']):
            output_filename = f'dynamic_structure_{atom_type}_{i}.txt'
            self._write_results(folder, output_filename, nu, Sv[i])

    def _write_integration_dynamic_structure(self, folder, atom_type, nu, S_intgr):
        """Write results for integrated dynamic structure."""
        for i in range(len(self.parameters['integration_list'])):
            output_filename = f'Integration_dynamic_structure_{atom_type}_{i}.txt'
            self._write_results(folder, output_filename, nu, S_intgr[i])

    def _write_vacf_and_pdos(self, atom_type, nu, pdos, t, vacf_non, vacf_output, folder):
        """Write the results of VACF and PDOS computations."""
        self._write_results(folder, f'atom{atom_type}_pdos_{self.parameters["Nc"]}.txt', nu, pdos)
        self._write_results(folder, f'atom{atom_type}_vacf_non_{self.parameters["Nc"]}.txt', t, vacf_non)
        self._write_results(folder, f'atom{atom_type}_vacf_scale_{self.parameters["Nc"]}.txt', t, vacf_output)


    def _write_dynamic_facilitation(self,folder,r,mean_Shells_record,mean_Shells1_record):
        dyna_fact_write= WriteFile(folder,self.parameters)
        dyna_fact_write.write_facilitation(r,mean_Shells_record,mean_Shells1_record)




