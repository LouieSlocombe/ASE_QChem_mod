import numpy as np
from ase.calculators.calculator import FileIOCalculator
from ase.calculators.calculator import SCFError
import ase.units


def make_neo_basis(neo_basis, neo_idx, neo_exp, atoms):
    neo_basis_list = neo_basis.split("_")
    ele = atoms.get_chemical_symbols()
    out_str = "$neo_basis \n"
    for idx in neo_idx:
        # write the atom name and index
        out_str += ele[idx] + " " + str(idx + 1) + "\n"
        for i in range(len(neo_basis_list)):
            for j in range(int(neo_basis_list[i][0])):
                out_str += neo_basis_list[i][-1].upper() + " 1 1.0\n"
                out_str += "  " + str(neo_exp[j]) + " 1.0\n"
        out_str += "****\n"
    out_str += "$end"
    return out_str


def make_neo_basis_presets(neo_basis_name, neo_idx, atoms):
    neo_basis_dict = {"PB4-D": "4s_3p_2d",
                      "PB4-F1": "4s_3p_2d_1f",
                      "PB4-F2": "4s_3p_2d_2f",
                      "PB5-D": "5s_4p_3d",
                      "PB5-F": "5s_4p_3d_2f",
                      "PB5-G": "5s_4p_3d_2f_1g",
                      "PB6-D": "6s_5p_4d",
                      "PB6-F": "6s_5p_4d_3f",
                      "PB6-G": "6s_5p_4d_3f_2g",
                      "PB6-H": "6s_5p_4d_3f_2g_1h",
                      }

    neo_basis_terms_dict = {"PB4-D": [1.957, 8.734, 16.010, 31.997, 9.438, 13.795, 24.028, 10.524, 19.016],
                            "PB4-F1": [5.973, 10.645, 17.943, 28.950, 7.604, 14.701, 23.308, 9.011, 19.787, 10.914],
                            "PB4-F2": [5.973, 10.645, 17.943, 28.950, 7.604, 14.701, 23.308, 9.011, 19.787, 10.914,
                                       20.985],
                            "PB5-D": [1.908, 9.051, 15.051, 29.766, 40.135, 4.907, 10.088, 15.893, 25.774, 10.352,
                                      18.358,
                                      36.366],
                            "PB5-F": [4.189, 6.231, 14.624, 20.481, 48.509, 2.349, 7.597, 18.521, 30.596, 8.971, 17.956,
                                      21.299, 10.321, 26.910],
                            "PB5-G": [3.283, 7.613, 16.077, 20.457, 47.011, 3.167, 7.785, 18.985, 29.425, 10.239,
                                      17.298,
                                      20.682, 10.930, 25.972, 10.513],
                            "PB6-D": [2.513, 4.840, 9.088, 16.231, 31.110, 38.754, 2.174, 5.897, 9.996, 16.347, 25.716,
                                      2.930, 11.015, 17.256, 25.541],
                            "PB6-F": [2.812, 5.703, 10.202, 15.450, 29.135, 39.753, 1.417, 6.884, 10.450, 16.187,
                                      26.447,
                                      4.483, 9.870, 18.313, 25.347, 6.001, 10.322, 24.409],
                            "PB6-G": [2.016, 3.049, 6.350, 8.910, 19.500, 29.296, 2.638, 3.716, 5.978, 16.252, 25.114,
                                      2.205, 3.654, 9.104, 24.039, 3.845, 9.921, 22.570, 10.062, 24.267],
                            "PB6-H": [1.386, 3.249, 9.562, 12.485, 21.429, 36.931, 1.399, 4.240, 9.930, 18.376, 24.119,
                                      2.607, 5.141, 7.750, 20.768, 0.509, 9.129, 26.408, 9.445, 28.407, 10.193],
                            }

    # Get the orbital types
    neo_basis_list = neo_basis_dict[neo_basis_name].split("_")
    # Get the elements
    ele = atoms.get_chemical_symbols()
    out_str = "$neo_basis \n"
    for idx in neo_idx:
        # write the atom name and index
        out_str += ele[idx] + " " + str(idx + 1) + "\n"
        # Get the exponents
        iter_exp = iter(neo_basis_terms_dict[neo_basis_name])
        for i in range(len(neo_basis_list)):
            for j in range(int(neo_basis_list[i][0])):
                out_str += neo_basis_list[i][-1].upper() + " 1 1.0\n"
                out_str += "  " + str(next(iter_exp)) + " 1.0\n"
        out_str += "****\n"
    out_str += "$end"
    return out_str


class QChem(FileIOCalculator):
    """
    QChem calculator
    """
    name = 'QChem'

    implemented_properties = ['energy', 'forces']
    command = 'qchem PREFIX.inp PREFIX.out'

    # Following the minimal requirements given in
    # http://www.q-chem.com/qchem-website/manual/qchem43_manual/sect-METHOD.html
    default_parameters = {'method': 'hf',
                          'basis': '6-31G*',
                          'jobtype': None,
                          'charge': 0}

    def __init__(self, restart=None,
                 ignore_bad_restart_file=FileIOCalculator._deprecated,
                 label='qchem', scratch=None, np=1, nt=1, pbs=False,
                 basisfile=None, ecpfile=None, atoms=None, neo_idx=[0],
                 neo_preset=None, neo_exp=[4.0, 8.0], neo_basis="1s_1p", solv_extra=None, **kwargs):
        """
        The scratch directory, number of processor and threads as well as a few
        other command line options can be set using the arguments explained
        below. The remaining kwargs are copied as options to the input file.
        The calculator will convert these options to upper case
        (Q-Chem standard) when writing the input file.

        scratch: str
            path of the scratch directory
        np: int
            number of processors for the -np command line flag
        nt: int
            number of threads for the -nt command line flag
        pbs: boolean
            command line flag for pbs scheduler (see Q-Chem manual)
        basisfile: str
            path to file containing the basis. Use in combination with
            basis='gen' keyword argument.
        ecpfile: str
            path to file containing the effective core potential. Use in
            combination with ecp='gen' keyword argument.
        """

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

        # Augment the command by various flags
        if pbs:
            self.command = 'qchem -pbs '
        else:
            self.command = 'qchem '
        if np != 1:
            self.command += '-np %d ' % np
        if nt != 1:
            self.command += '-nt %d ' % nt
        self.command += 'PREFIX.inp PREFIX.out'
        if scratch is not None:
            self.command += ' %s' % scratch

        # Redirect output to PREFIX.log
        self.command += " >> PREFIX.log"

        self.basisfile = basisfile
        self.ecpfile = ecpfile
        self.neo_idx = neo_idx
        self.neo_exp = neo_exp
        self.neo_basis = neo_basis
        self.neo_preset = neo_preset
        self.solv_extra = solv_extra

    def read(self, label):
        raise NotImplementedError

    def read_results(self):
        filename = self.label + '.out'

        with open(filename, 'r') as fileobj:
            lineiter = iter(fileobj)
            N_atoms = self.atoms.get_global_number_of_atoms()
            e_conv = ase.units.Hartree
            f_conv = ase.units.Hartree / ase.units.Bohr
            for line in lineiter:
                if 'SCF failed to converge' in line:
                    raise SCFError()
                elif 'ERROR: alpha_min' in line:
                    # Even though it is not technically a SCFError:
                    raise SCFError()
                elif ' Total energy in the final basis set =' in line:
                    self.results['energy'] = float(line.split()[8]) * e_conv
                elif ' Total energy =' in line:
                    self.results['energy'] = float(line.split()[3]) * e_conv
                elif ' E(NEO-SCF) =' in line:
                    # Adjust the energy from the NEO-SCF calculation
                    self.results['energy'] = float(line.split("=")[1]) * e_conv
                elif ' Gradient of SCF Energy' in line:
                    # Read gradient as 3 by N array and transpose at the end
                    gradient = [[] for _ in range(3)]
                    # Skip first line containing atom numbering
                    next(lineiter)
                    while True:
                        # Loop over the three Cartesian coordinates
                        for i in range(3):
                            # Cut off the component numbering and remove
                            # trailing characters ('\n' and stuff)
                            line = next(lineiter)[5:].rstrip()
                            # Cut in chunks of 12 symbols and convert into
                            # strings. This is preferred over string.split() as
                            # the fields may overlap for large gradients
                            gradient[i].extend(list(map(
                                float, [line[i:i + 12]
                                        for i in range(0, len(line), 12)])))

                        # After three force components we expect either a
                        # separator line, which we want to skip, or the end of
                        # the gradient matrix which is characterized by the
                        # line ' Max gradient component'.
                        # Maybe change stopping criterion to be independent of
                        # next line. Eg. if not lineiter.next().startswith(' ')
                        if ' Max gradient component' in next(lineiter):
                            # Minus to convert from gradient to force
                            self.results['forces'] = -np.array(gradient).T * f_conv
                            break
                elif ' NEO-SCF Analytic gradient:' in line:
                    # Read gradient as N array by 3
                    gradient = np.zeros((N_atoms, 3))
                    # Skip first line containing atom numbering
                    next(lineiter)
                    # Loop over the number of atoms
                    for i in range(N_atoms):
                        gradient[i, :] = np.fromstring(next(lineiter), dtype=float, sep=' ')[1:]
                    # Minus to convert from gradient to force
                    self.results['forces'] = -np.array(gradient) * f_conv

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        filename = self.label + '.inp'

        with open(filename, 'w') as fileobj:
            fileobj.write('$comment\n   ASE generated input file\n$end\n\n')
            # set main job type
            fileobj.write('$rem\n')
            if self.parameters['jobtype'] is None:
                if 'forces' in properties:
                    fileobj.write('   %-25s   %s\n' % ('JOBTYPE', 'FORCE'))
                else:
                    fileobj.write('   %-25s   %s\n' % ('JOBTYPE', 'SP'))

            # Write all parameters
            for prm in self.parameters:
                if prm not in ['charge', 'multiplicity']:
                    if self.parameters[prm] is not None:
                        # if "neo_".upper() not in prm.upper():
                        fileobj.write('   %-25s   %s\n' % (
                            prm.upper(), self.parameters[prm].upper()))

            # Not even a parameters as this is an absolute necessity
            fileobj.write('   %-25s   %s\n' % ('SYM_IGNORE', 'TRUE'))
            fileobj.write('$end\n\n')

            # Write the molecule block
            fileobj.write('$molecule\n')
            # Following the example set by the gaussian calculator
            if ('multiplicity' not in self.parameters):
                tot_magmom = atoms.get_initial_magnetic_moments().sum()
                mult = tot_magmom + 1
            else:
                mult = self.parameters['multiplicity']
            # Default charge of 0 is defined in default_parameters
            fileobj.write('   %d %d\n' % (self.parameters['charge'], mult))
            for a in atoms:
                fileobj.write('   %s  %f  %f  %f\n' % (a.symbol,
                                                       a.x, a.y, a.z))
            fileobj.write('$end\n\n')

            if self.basisfile is not None:
                with open(self.basisfile, 'r') as f_in:
                    basis = f_in.readlines()
                fileobj.write('$basis\n')
                fileobj.writelines(basis)
                fileobj.write('$end\n\n')

            if self.ecpfile is not None:
                with open(self.ecpfile, 'r') as f_in:
                    ecp = f_in.readlines()
                fileobj.write('$ecp\n')
                fileobj.writelines(ecp)
                fileobj.write('$end\n\n')

            # Write the solvent block
            if 'solvent_method' in self.parameters:
                if self.solv_extra is not None:
                    fileobj.write("\n")
                    fileobj.write(self.solv_extra)
                    fileobj.write("\n")

            # Write the neo basis block
            if 'neo' in self.parameters:
                if self.parameters['neo'].upper() == 'TRUE':
                    if self.neo_preset is None:
                        fileobj.write("\n")
                        out_str = make_neo_basis(self.neo_basis, self.neo_idx, self.neo_exp, atoms)
                        fileobj.write(out_str)
                        fileobj.write("\n")
                    else:
                        fileobj.write("\n")
                        out_str = make_neo_basis_presets(self.neo_preset, self.neo_idx, atoms)
                        fileobj.write(out_str)
                        fileobj.write("\n")
