'''
    ================================================================================================
    Core objects:  --- :mod: 'bilana.core.energy
        - Energy
            - run_calculation                     Reruns the MD simulations for specified resids to
                                                  extract the respective interaction energies
                                                  to lipid neighbors
            - run_lip_leaflet_interaction         Reruns MD simulations to extract lipid
                                                  leaflet energy
        - check_exist_xvgs                        Checks wether all xvg files from mdrun -rerun's
                                                  have correct length and exist
        - write_energyfile                        Gather all interaction energies from xvg files
                                                  into one single file
        - write_selfinteractionfile               Same as write_energyfile but for lipid
                                                  selfinteraction
        - create_lipid_water_interaction_file     Write file comprising all lipid water interaction
        - create_lipid_leaflet_interaction_file   Same as above for lipid leaflet interaction
        - create_indexfile                        Creates an indexfile used for mdrun -rerun
    ================================================================================================
'''
import os
import logging

import pandas as pd

from .systeminfo import Systeminfo
from ..lib.gromacswrapper import GMXNAME, exec_gromacs, write_log

LOGGER = logging.getLogger("bilana2.core.energy")


class Energy(Systeminfo):
    '''
        This class stores all information that is needed to automate the calculation of
        lipid interaction energies
        The main function is
            run_calculation
        For more information see its docstring
    '''

    __NLIPIDS_MAX_PER_CALCULATION = 40

    def __init__(self,
                 part,
                 neiblist,
                 inputfilepath="inputfile",
                 resindex_all='resindex_all.ndx',
                 energyname="all_energies.dat",
                 overwrite=True,
                 ):

        super().__init__(inputfilepath)
        knownparts = ['complete', 'head-tail']

        self.path.create_folders()

        if part not in knownparts:
            raise ValueError("Part keyword specified is not known.")
        if not os.path.isfile(resindex_all):
            raise FileNotFoundError("Missing resindex_all file. Create it using create_indexfile")

        self.neiblist      = neiblist
        self.resindex_all  = resindex_all
        self.overwrite     = overwrite
        self.groupblocks   = None  # Will be set in run_calculation function
        self.part          = part

        if part == 'complete':
            self.molparts            = ["resid_"]   # This is used as a descriptor for the index grp
            self.part                = ''           #
            self.max_lipids_per_calc = self.__NLIPIDS_MAX_PER_CALCULATION
            self.molparts_short      = [""]
            self.all_energies        = energyname
        elif part == 'head-tail':
            self.molparts            = ["resid_h_", "resid_t_"]
            self.max_lipids_per_calc = int( self.__NLIPIDS_MAX_PER_CALCULATION / 2 )
            self.molparts_short      = ["h_", "t_"]
            self.all_energies        = energyname

        LOGGER.info('\n Calculating for energygroups: %s', self.molparts)

    def run_calculation(self, resids):
        ''' Runs an energy calculation with settings from Energy() instance.
            For each residue the energy to all neighbors seen during MD is calculated
            and written to .edr files.
            Procedure is as follows:
            1. The neighbors are divided into fragments ("groupfragments")
            2. For each fragment:
                an mdp file is created (create_MDP)
                a tpr file is generated (create_TPR)
            3. The actual mdrun -rerun is performed (do_Energyrun)
            4. .xvg tables are generate from .edr files

        '''
        LOGGER.info('Rerunning MD for energyfiles')
        assert all( isinstance(x, int) for x in resids ), "resids have to be int"
        for res in resids:
            LOGGER.info('Working on lipid %s ...', res)

            ### Get all neighbors of res ###
            all_neibs_of_res = list(set([neibs for t in self.neiblist.keys()
                    for neibs in self.neiblist[t][res]]))
            nneibs = len(all_neibs_of_res)

            ### Split neighborlist such that max_lipids_per_calc is not exceeded ###
            if nneibs % self.max_lipids_per_calc == 0:
                number_of_groupfragments = ( nneibs // self.max_lipids_per_calc )
            else:
                number_of_groupfragments = ( nneibs // self.max_lipids_per_calc ) + 1
            LOGGER.info("Needing %s energy run(s)", number_of_groupfragments)

            for groupfragment in range(number_of_groupfragments):
                LOGGER.info("... on fragment %s ...", groupfragment)

                groupblockstart  = groupfragment * self.max_lipids_per_calc
                groupblockend    = (groupfragment+1) * self.max_lipids_per_calc
                self.groupblocks = (groupblockstart, groupblockend)

                # File in-/outputs
                groupfragment   = str(groupfragment)

                g_energy_output = '{}/xvgtables/energies_residue{}_{}{}.xvg'.format(
                                    self.path.energy, str(res), str(groupfragment), self.part)

                mdpout          = "{}/mdpfiles/energy_mdp_recalc_resid{}_{}{}.mdp".format(
                                    self.path.energy, str(res), groupfragment, self.part)

                tprout          = "{}/tprfiles/mdrerun_resid{}_{}{}.tpr".format(
                                    self.path.energy, str(res), groupfragment, self.part)

                energyf_output  = "{}/edrfiles/energyfile_resid{}_{}{}.edr".format(
                                    self.path.energy, str(res), groupfragment, self.part)

                xvg_out         = "{}/xvgtables/energies_residue{}_{}{}.xvg".format(
                                    self.path.energy, str(res), groupfragment, self.part)

                energygroups    = self._gather_energygroups(res, all_neibs_of_res)
                relev_energies  = self._get_relev_energies(res, all_neibs_of_res)

                # Self interactions
                relev_energies_self = self._get_relev_self_interaction(res)
                tprout_self         = "{}/tprfiles/mdrerun_resid{}_0{}.tpr"\
                                        .format(self.path.energy, str(res), self.part)
                energyf_output_self = "{}/edrfiles/energyfile_resid{}_0{}.edr"\
                                        .format(self.path.energy, str(res), self.part)
                xvg_out_self        = "{}/xvgtables/energies_residue{}_selfinteraction{}.xvg"\
                                        .format(self.path.energy, str(res), self.part)

                # Run functions
                self.create_MDP(mdpout, energygroups)
                self.create_TPR(mdpout, tprout)
                if os.path.isfile(energyf_output) and not self.overwrite:
                    LOGGER.info("Edrfile for lipid %s part %s already exists."
                                "Will skip this calculation.", res, groupfragment)
                else:
                    self.do_Energyrun(res, groupfragment, tprout, energyf_output)
                if os.path.isfile(g_energy_output) and not self.overwrite:
                    LOGGER.info("Xvgtable for lipid %s part %s already exists."
                                "Will skip this calculation.", res, groupfragment)
                else:
                    self.write_XVG(energyf_output, tprout, relev_energies, xvg_out)
                    self.write_XVG(energyf_output_self,
                        tprout_self, relev_energies_self, xvg_out_self)
        return 1

    def run_lip_leaflet_interaction(self, resids):
        ''' Run the energy calculation between each lipid and the opposing leaflet  '''
        LOGGER.info('Rerunning MD for energyfiles')

        for res in resids:
            LOGGER.info('Working on lipid %s ...', res)
            leaflet = self.convert.resid_to_leaflet[res]

            ### Get all leaflets from other leaflets ###
            ### ! This takes the leaflet assignment from frame 1, though  ! ###
            ### ! it should be dynamic, thus for future purpose this loop ! ###
            ### ! should remain inside the loop over resids               ! ###
            res_other_leaflet = []
            res_same_leaflet = []
            for nres in self.lipid_resids:
                leaf = self.convert.resid_to_leaflet[nres]
                if leaf != leaflet:
                    res_other_leaflet.append(nres)
                else:
                    if nres == res:
                        continue
                    res_same_leaflet.append(nres)

            ### Create an index file containing the reference molecule and ###
            ### indices of all mols in the opposing leaflet                ###
            outputndx = "{}/energy_leaflet_res{}.ndx".format(self.path.tmp, res)
            self.resindex_all = outputndx
            selectionstr = 'System=all;'\
                           'interleaflet=resname {resnames} and resid {interlipids};'\
                           'leaflet=resname {resnames} and resid {leafletlipids};'\
                           'resid_{hostid}=resid {hostid};'\
                           'resid_{hostid}; interleaflet; leaflet; System;'\
                           .format(resnames=' '.join(self.molecules),
                                interlipids=' '.join([str(i) for i in res_other_leaflet]),
                                leafletlipids=' '.join([str(i) for i in res_same_leaflet]),
                                hostid=res,
                           )

            cmd = [
                "-f", self.path.gro,
                "-s", self.path.tpr, "-select", selectionstr,
                "-on", outputndx,
                ]
            out, err = exec_gromacs(GMXNAME, "select", cmd)
            write_log("gmx_select", out, err, path=self.path.log)

            ### Defining input and output file paths ###
            g_energy_output = '{}/xvgtables/energies_residue{}_leaflet.xvg'\
                                .format(self.path.energy, str(res))
            mdpout          = '{}/mdpfiles/energy_mdp_recalc_resid{}_leaflet.mdp'\
                                .format(self.path.energy, str(res))
            tprout          = '{}/tprfiles/mdrerun_resid{}_leaflet.tpr'\
                                .format(self.path.energy, str(res))
            energyf_output  = '{}/edrfiles/energyfile_resid{}_leaflet.edr'\
                                .format(self.path.energy, str(res))
            xvg_out         = '{}/xvgtables/energies_residue{}_leaflet.xvg'\
                                .format(self.path.energy, str(res))
            energygroups    = "resid_{} leaflet interleaflet".format(res)
            relev_energies  = '\n'.join(
                ["Coul-SR:resid_{}-leaflet".format(res), "LJ-SR:resid_{}-leaflet".format(res),
                 "Coul-SR:resid_{}-interleaflet".format(res), "LJ-SR:resid_{}-interleaflet".format(res),
                '\n']
                )

            # Run functions
            self.create_MDP(mdpout, energygroups)
            self.create_TPR(mdpout, tprout)
            if os.path.isfile(energyf_output) and not self.overwrite:
                LOGGER.info("Edrfile for lipid %s part %s already exists. "
                            "Will skip this calculation.", res, self.part)
            else:
                self.do_Energyrun(res, 0, tprout, energyf_output)
            if os.path.isfile(g_energy_output) and not self.overwrite:
                LOGGER.info("Xvgtable for lipid %s part %s already exists. "
                            "Will skip this calculation.", res, self.part)
            else:
                self.write_XVG(energyf_output, tprout, relev_energies, xvg_out)
        return 1

    def create_MDP(self, mdpout: str, energygroups: str):
        ''' Create mdpfile '''
        os.makedirs(self.path.energy + '/mdpfiles', exist_ok=True)
        with open(mdpout, "w") as mdpfile_rerun:
            raw_mdp = [x.strip() for x in '''
            integrator              = md
            dt                      = 0.002
            nsteps                  =
            nstlog                  = 100000
            nstxout                 = 0
            nstvout                 = 0
            nstfout                 = 0
            nstcalcenergy           = 1000
            nstenergy               = 100
            cutoff-scheme           = Verlet
            nstlist                 = 20
            rlist                   = 1.2
            coulombtype             = pme
            rcoulomb                = 1.2
            vdwtype                 = Cut-off
            vdw-modifier            = Force-switch
            rvdw_switch             = 1.0
            rvdw                    = 1.2
            tcoupl                  = Nose-Hoover
            tau_t                   = 1.0
            tc-grps                 = System
            pcoupl                  = Parrinello-Rahman
            pcoupltype              = semiisotropic
            tau_p                   = 5.0
            compressibility         = 4.5e-5  4.5e-5
            ref_p                   = 1.0     1.0
            constraints             = h-bonds
            constraint_algorithm    = LINCS
            continuation            = yes
            nstcomm                 = 100
            comm_mode               = linear
            refcoord_scaling        = com
            '''.split('\n')]
            raw_mdp.append( 'ref_t = ' + str(self.temperature) )
            energygrpline = ''.join(['energygrps\t\t\t=', energygroups, '\n'])
            raw_mdp.append(energygrpline)
            mdpfile_rerun.write('\n'.join(raw_mdp)+'\n')

    def create_TPR(self, mdpoutfile: str, tprout: str):
        ''' Create TPRFILE with GROMPP '''
        os.makedirs(self.path.energy+'/tprfiles', exist_ok=True)

        grompp_arglist = ['-f', mdpoutfile, '-p',
                        self.path.top, '-c', self.path.gro, '-o', tprout,
                        '-n', self.resindex_all, '-po', mdpoutfile
                        ]
        out, err = exec_gromacs(GMXNAME, "grompp", grompp_arglist)

        write_log("gmx_grompp", out, err, path=self.path.log)

    def do_Energyrun(self, res, groupfragment, tprrerun_in, energyf_out):
        ''' Create .edr ENERGYFILE with mdrun -rerun '''
        LOGGER.info('...Rerunning trajectory for energy calculation...')

        os.makedirs(self.path.energy+'/edrfiles', exist_ok=True)
        os.makedirs(self.path.energy+'/logfiles', exist_ok=True)

        logoutput_file = "{}/logfiles/mdrerun_resid{}{}frag{}.log"\
                        .format( self.path.energy, str(res), self.part, str(groupfragment) )
        trajout        = 'DUMMY.trr' # As specified in mdpfile, !NO! .trr-file should be written

        mdrun_arglist = ['-s', tprrerun_in, '-rerun', self.path.trj,
                         '-e', energyf_out, '-o', trajout,'-g', logoutput_file,
                        ]
        out, err = exec_gromacs(GMXNAME, "mdrun", mdrun_arglist)

        write_log("gmx_mdrun", out, err, path=self.path.log)

    def write_XVG(self, energyf_in, tprrerun_in, all_relev_energies, xvg_out):
        ''' Create XVG-TABLE with all relevant energies '''
        os.makedirs(self.path.energy+'/xvgtables', exist_ok=True)

        g_energy_arglist = ['-f', energyf_in,
                            '-s', tprrerun_in,'-o', xvg_out,
                           ]
        inp_str=all_relev_energies.encode()

        out, err = exec_gromacs(GMXNAME, "energy", g_energy_arglist, inp_str)

        write_log("gmx_energy", out, err, path=self.path.log)

    def _gather_energygroups(self, res, all_neibs_of_res):
        ''' Set which part of molecule should be considered '''
        energygroup_resids = [res] + all_neibs_of_res[ self.groupblocks[0]:self.groupblocks[1] ]
        energygroup_list = []

        for resid in energygroup_resids:
            resname = self.convert.resid_to_resname[resid]
            if not ( self.ff.is_sterol(resname) or self.ff.is_protein(resname) ):
                energygroup_list.append( ''.join( ["resid_", str(resid)] ) )
            else:
                for part in self.molparts:
                    energygroup_list.append( ''.join( [part, str(resid)] ) )

        # Add water groups for molecule - solvent interaction
        energygroup_list += ["solv"]

        energygroup_string = ' '.join(energygroup_list)
        return energygroup_string

    def _get_relev_energies(self, res, all_neibs_of_res):
        '''
            Returns a string that contains all entries to be extracted from the energy file
            using gmx energy.
            This version is for lipid-lipid interaction
            for self interaction search function "get_relev_self_interaction"
        '''
        Etypes=["Coul-SR:", "LJ-SR:"]
        energyselection=[]
        resname = self.convert.resid_to_resname[res]
        host_is_pl = not ( self.ff.is_sterol(resname) or self.ff.is_protein(resname) )

        for interaction in Etypes:

            for neib in all_neibs_of_res[ self.groupblocks[0]:self.groupblocks[1] ]:
                resn_n = self.convert.resid_to_resname[neib]
                neib_is_pl = not ( self.ff.is_sterol(resn_n) or self.ff.is_protein(resn_n) )

                ### The following if clause is needed as a distinction between PL and non PL ###
                ### molecules, as the separation in different parts makes only sense for PL  ###
                ### Thus for sterols always the whole lipid part is used                     ###
                if host_is_pl and neib_is_pl:
                    for parthost in self.molparts:
                        for partneib in self.molparts:
                            energyselection.append( "{}{}{}-{}{}".format(
                                    interaction, parthost, str(res), partneib, str(neib) ) )
                elif host_is_pl and not neib_is_pl:
                    partneib = "resid_"
                    for parthost in self.molparts:
                        energyselection.append( "{}{}{}-{}{}".format(
                                interaction, parthost, str(res), partneib, str(neib) ) )
                elif not host_is_pl and neib_is_pl:
                    parthost = "resid_"
                    for partneib in self.molparts:
                        energyselection.append( "{}{}{}-{}{}".format(
                                interaction, parthost, str(res), partneib, str(neib) ) )
                elif not host_is_pl and not neib_is_pl:
                    parthost = "resid_"
                    partneib = "resid_"
                    energyselection.append( "{}{}{}-{}{}".format(
                            interaction, parthost, str(res), partneib, str(neib) ) )

        res_solv_interaction = [ "{}resid_{}-solv".format(etype, res) for etype in Etypes ]

        all_relev_energies = '\n'.join( res_solv_interaction + energyselection + ['\n'] )
        return all_relev_energies

    def _get_relev_self_interaction(self, res):
        ''' Returns string that describes all entries
            needed to be extracted from energy file using gmx energy
            This version is for lipid self interaction
        '''

        Etypes=["Coul-SR:", "LJ-SR:", "Coul-14:", "LJ-14:"]

        energyselection=[]
        for interaction in Etypes:
            for parthost in self.molparts:
                energyselection.append("{0}{1}{2}-{1}{2}".format(interaction, parthost, str(res)))
        all_relev_energies='\n'.join( energyselection + ['\n'] )
        return all_relev_energies

def create_lipid_leaflet_interaction_file(sysinfo, outputfilename="resid_leaflet_interaction.dat"):
    ''' Create a file with entries of
        interaction of resid at time to leaflet0 and leaflet1
        <time> <resid> <resname> <host_leaflet> <Etot> <Evdw> <Ecoul>
    '''
    energyoutput = open(outputfilename, "w")
    print('{: <10}{: <10}{: <10}{: <10}{: <20}{: <20}{: <20}'\
        .format("time", "resid", "resname", "leaflet_h", "Etot", "Evdw", "Ecoul"),
            file=energyoutput)
    for resid in sysinfo.lipid_resids:
        resname = sysinfo.convert.resid_to_resname[resid]
        leaflet = sysinfo.convert.resid_to_leaflet[resid]
        xvgfilename = "{}/xvgtables/energies_residue{}_leaflet.xvg"\
                .format( sysinfo.path.energy, str(resid) )

        with open(xvgfilename,"r") as xvgfile:
            res_to_rowindex = {}

            ### folderlayout is: <time> <Coul_resHost_resNeib> <LJ_resHost_resNeib> ... ###
            for energyline in xvgfile:
                energyline_cols = energyline.split()

                ### creating a dict to know which column(energies) belong to which residue ###
                if '@ s' in energyline:
                    rowindex  = int(energyline_cols[1][1:])+1 # time is at row 0 !
                    host = energyline_cols[3].split("resid_")[1].split("-")[0]
                    energytype = energyline_cols[3].split("-")[0][1:]
                    ### leaflet type can be either "leaflet" or "interleaflet" ###
                    leaflet_type = energyline_cols[3].split("-")[2][:-1]
                    LOGGER.debug("Hostid: %s:", host)

                    res_to_rowindex[(energytype, leaflet_type, host)] = rowindex
                    LOGGER.debug("Adding to dict: Etype %s, leaflettype: %s, host %s", energytype, leaflet_type, host)

                ### pick correct energies from energyfile and print ###
                elif '@' not in energyline and '#' not in energyline:
                    time = float(energyline_cols[0])
                    if time % sysinfo.dt != 0:
                        continue

                    vdw  = float( energyline_cols[ res_to_rowindex[ ( 'LJ', "leaflet", str(resid) ) ] ] )
                    coul = float( energyline_cols[ res_to_rowindex[ ( 'Coul', "leaflet", str(resid) ) ] ] )
                    Etot = vdw + coul
                    outpline = '{: <10}{: <10}{: <10}{: <10}{: <20.5f}{: <20.5f}{: <20.5f}'\
                        .format(time, resid, resname, leaflet, Etot, vdw, coul,)
                    print(outpline, file=energyoutput)
    energyoutput.close()


def create_lipid_interleaflet_interaction_file(sysinfo, outputfilename="resid_interleaflet_interaction.dat"):
    ''' Create a file with entries of
        interaction of resid at time to leaflet0 and leaflet1
        <time> <resid> <resname> <host_leaflet> <Etot> <Evdw> <Ecoul>
    '''
    energyoutput = open(outputfilename, "w")
    print('{: <10}{: <10}{: <10}{: <10}{: <20}{: <20}{: <20}'\
        .format("time", "resid", "resname", "leaflet_h", "Etot", "Evdw", "Ecoul"),
            file=energyoutput)
    for resid in sysinfo.lipid_resids:
        resname = sysinfo.convert.resid_to_resname[resid]
        leaflet = sysinfo.convert.resid_to_leaflet[resid]
        xvgfilename = "{}/xvgtables/energies_residue{}_leaflet.xvg"\
                .format( sysinfo.path.energy, str(resid) )

        with open(xvgfilename,"r") as xvgfile:
            res_to_rowindex = {}

            ### folderlayout is: <time> <Coul_resHost_resNeib> <LJ_resHost_resNeib> ... ###
            for energyline in xvgfile:
                energyline_cols = energyline.split()

                ### creating a dict to know which column(energies) belong to which residue ###
                if '@ s' in energyline:
                    rowindex  = int(energyline_cols[1][1:])+1 # time is at row 0 !
                    host = energyline_cols[3].split("resid_")[1].split("-")[0]
                    energytype = energyline_cols[3].split("-")[0][1:]
                    ### leaflet type can be either "leaflet" or "interleaflet" ###
                    leaflet_type = energyline_cols[3].split("-")[2][:-1]
                    LOGGER.debug("Hostid: %s:", host)

                    res_to_rowindex[(energytype, leaflet_type, host)] = rowindex
                    LOGGER.debug("Adding to dict: Etype %s, leaflettype: %s, host %s", energytype, leaflet_type, host)

                ### pick correct energies from energyfile and print ###
                elif '@' not in energyline and '#' not in energyline:
                    time = float(energyline_cols[0])
                    if time % sysinfo.dt != 0:
                        continue

                    vdw  = float( energyline_cols[ res_to_rowindex[ ( 'LJ', "interleaflet", str(resid) ) ] ] )
                    coul = float( energyline_cols[ res_to_rowindex[ ( 'Coul', "interleaflet", str(resid) ) ] ] )
                    Etot = vdw + coul
                    outpline = '{: <10}{: <10}{: <10}{: <10}{: <20.5f}{: <20.5f}{: <20.5f}'\
                        .format(time, resid, resname, leaflet, Etot, vdw, coul,)
                    print(outpline, file=energyoutput)
    energyoutput.close()

def write_selfinteractionfile(energy):
    ''' Extracts all self interaction energy entries from xvg files
        and writes them to "selfinteractions.dat"
    '''
    with open("selfinteractions.dat", "w") as energyoutput:
        print(\
              '{: <10}{: <10}{: <10}'
              '{: <20}{: <20}{: <20}{: <20}{: <20}{: <20}{: <20}'\
              .format("time", "resid", "resname",
                      "Etot", "VdWSR", "CoulSR", "VdW14", "Coul14", "VdWtot", "Coultot", ),
              file=energyoutput)

        for resid in energy.lipid_resids:
            xvg_out = '{}/xvgtables/energies_residue{}_selfinteraction{}.xvg'\
                .format(energy.path.energy, str(resid), energy.part)
            resname = energy.convert.resid_to_resname[resid]

            with open(xvg_out,"r") as xvgfile:
                res_to_row = {}

                ### folderlayout is: <time> <Coul_resHost_resNeib> <LJ_resHost_resNeib> ... ###
                for energyline in xvgfile:
                    energyline_cols = energyline.split()

                    ### creating a dict to know which column(energies) belong to which residue ###
                    if '@ s' in energyline:
                        ### time is at row 0 ! ###
                        row         = int(energyline_cols[1][1:]) + 1
                        host        = energyline_cols[3].split("resid_")[1][:-1]
                        energytype  = energyline_cols[3].split(":")[0][1:]
                        res_to_row.update( { (energytype, host):row } )

                    ### pick correct energies from energyfile and print ###
                    elif '@' not in energyline and '#' not in energyline:
                        time = float(energyline_cols[0])
                        if time % energy.dt != 0:
                            continue

                        try:
                            vdw_sr   = energyline_cols[res_to_row[('LJ-SR', str(host))]]
                            vdw_14   = energyline_cols[res_to_row[('LJ-14', str(host))]]
                            coul_sr  = energyline_cols[res_to_row[('Coul-SR', str(host))]]
                            coul_14  = energyline_cols[res_to_row[('Coul-14', str(host))]]
                            vdw_tot  = float(vdw_14) + float(vdw_sr)
                            coul_tot = float(coul_14) + float(coul_sr)
                        except KeyError:
                            continue

                        Etot = float(vdw_sr) + float(coul_sr) + float(vdw_14) + float(coul_14)
                        print(
                              '{: <10}{: <10}{: <10}{: <20.5f}'
                              '{: <20}{: <20}{: <20}{: <20}{: <20.5f}{: <20.5f}'
                              .format(time, resid, resname, Etot,
                                      vdw_sr, coul_sr, vdw_14, coul_14, vdw_tot, coul_tot,),
                              file=energyoutput)

def check_exist_xvgs(energy,
                     check_len=False,
                     missing_file_log="missing_xvg.log"):
    ''' Checks if all .xvg-files containing lipid interaction exist

        check_len can be set to simulation length and all .xvg files that differ from that
        length are included to missing_xvgfile.info
    '''
    def read_lastline_only(fname):
        ''' Reads last line of file without loop. Attention: Crashes if file is empty'''
        maxcount = 100000
        cnt = 0
        with open(fname, "rb") as f:
            f.seek(-2, os.SEEK_END)     # Jump to the second last byte.
            while f.read(1) != b"\n" and f.tell() != 1:   # Until EOL is found... if f.tell() == 0 means cursor is at the beginning of file
                f.seek(-2, os.SEEK_CUR) # ...jump back the read byte plus one more.
                cnt += 1
                if cnt > maxcount:
                    raise RuntimeError("Reach max byte count in file {}."\
                            "Is it corrupted?".format(fname))

            if f.tell() == 1:
                raise RuntimeError("File {} has no EOL character. Is it corrupted?".format(fname))

            last = f.readline()         # Read last line.
        return last

    missing_files = []
    reasons       = []

    for resid in energy.lipid_resids:
        all_neibs_of_res = list(set([neibs for t in energy.neiblist.keys()\
            for neibs in energy.neiblist[t][resid]]))
        n_neibs = len(all_neibs_of_res)

        if n_neibs % energy.max_lipids_per_calc == 0:
            number_of_groupfragments = ( n_neibs // energy.max_lipids_per_calc )
        else:
            number_of_groupfragments = ( n_neibs // energy.max_lipids_per_calc ) + 1

        for part in range(number_of_groupfragments):
            xvgfilename = "{}/xvgtables/energies_residue{}_{}{}.xvg"\
                .format(energy.path.energy, str(resid), str(part), energy.part)

            if not os.path.isfile(xvgfilename):
                # Check wether file exists
                missing_files.append(xvgfilename)
                reasons.append("missing")

            elif os.stat(xvgfilename).st_size == 0:
                # Check wether file is empty
                missing_files.append(xvgfilename)
                reasons.append("empty")

            elif check_len:
                # Check wether whole trajectory was used
                lastline = read_lastline_only(xvgfilename)
                time = float(lastline.decode().split()[0])
                if time < check_len:
                    missing_files.append(xvgfilename)
                    reasons.append("incomplete")

    if missing_files:
        LOGGER.error(
            "There were missing .xvg files! Check %s/%s for further information",
            energy.path.energy, missing_file_log)
        dat = pd.DataFrame({"file":missing_files, "reason":reasons})
        dat.to_csv(missing_file_log, index=False)
        return False
    return True

def write_energyfile(energy):
    ''' Creates files: "all_energies_<interaction>.dat
        NOTE: This function is too long. It should be separated into smaller parts.
    '''
    missing_energydata = []
    outp_fstr = '{: <10}{: <10}{: <10}{: <20}{: <20}{: <20}{: <20.5f}'
    LOGGER.info('Create energy file')
    with open(energy.all_energies, "w") as energyoutput:
        print(
              '{: <10}{: <10}{: <10}{: <20}'
              '{: <20}{: <20}{: <20}'\
              .format("time", "host", "neighbor", "molparts",\
                                       "VdW", "Coul", "Etot"),\
              file=energyoutput)

        for resid in energy.lipid_resids:
            LOGGER.info("Working on residue %s ...", resid)

            # Get neighborhood of resid
            hosttype         = energy.convert.resid_to_resname[resid]
            host_is_pl       = not (energy.ff.is_sterol(hosttype) or energy.ff.is_protein(hosttype))
            all_neibs_of_res = list(set([neibs for t in energy.neiblist.keys()\
                    for neibs in energy.neiblist[t][resid]]))
            n_neibs          = len(all_neibs_of_res)
            LOGGER.debug("All neibs of res %s are %s", resid, all_neibs_of_res)

            # Get number of fragments (how many runs per residue)
            if n_neibs % energy.max_lipids_per_calc == 0:
                number_of_groupfragments = ( n_neibs // energy.max_lipids_per_calc )
            else:
                number_of_groupfragments = ( n_neibs // energy.max_lipids_per_calc ) + 1
            LOGGER.debug("Nneibs: %s Nfrags: %s", n_neibs, number_of_groupfragments)

            ### Following line cuts the list of neighbors into sublists that corresponds to ###
            ### the respective mdrun -rerun's                                               ###
            all_neibs_of_res = [ all_neibs_of_res[ ( i*energy.max_lipids_per_calc ):( (i+1)*energy.max_lipids_per_calc ) ] for i in range(number_of_groupfragments) ]

            for frag in range(number_of_groupfragments):
                LOGGER.info("At fragment %s", frag)
                xvgfilename = "{}/xvgtables/energies_residue{}_{}{}.xvg"\
                    .format(energy.path.energy, str(resid), str(frag), energy.part)

                with open(xvgfilename,"r") as xvgfile:
                    res_to_rowindex = {}

                    for energyline in xvgfile:
                        energyline_cols = energyline.split()

                        ### folderlayout is: <time> <Coul_resHost_resNeib> <LJ_resHost_resNeib> ...
                        if '@ s' in energyline:
                            rowindex    = int(energyline_cols[1][1:])+1 # time is at row 0 !
                            interaction_str = energyline_cols[3]
                            host  = interaction_str.split("-")[1].replace("SR:", "")
                            neib  = interaction_str.split("-")[2][:-1]

                            #host        = left.split("_")[-1]
                            #if not "solv" in right:
                            #    neib    =  right.split("_")[-1]
                            #else:
                            #    neib    = "solv"

                            energytype  = energyline_cols[3].split("-")[0][1:]
                            LOGGER.debug("Hostid: %s, Neibid: %s", host, neib)

                            res_to_rowindex[(energytype, host, neib)] = rowindex
                            LOGGER.debug("Adding to dict: Etype %s, host %s, neib %s",
                                    energytype, host, neib)

                        ### pick correct energies from energyfile and print ###
                        elif '@' not in energyline and '#' not in energyline:
                            time = float(energyline_cols[0])
                            if time % energy.dt != 0:
                                continue

                            for neib in all_neibs_of_res[frag]:

                                neibtype = energy.convert.resid_to_resname[neib]
                                neib_is_pl = not (energy.ff.is_sterol(neibtype)\
                                        or energy.ff.is_protein(neibtype))

                                ### In the following a distinction is made between PL and  ###
                                ### non-PL molecules as separation in head/tail parts is   ###
                                ### only reasonable for PL molecules                       ###
                                if host_is_pl and neib_is_pl:

                                    for parthost in energy.molparts:

                                        ### remove "resid_" from parthost string ###
                                        parthost_short = parthost[7:]
                                        if parthost_short == '':
                                            parthost_short = "w"

                                        for partneib in energy.molparts:

                                            partneib_short = parthost[7:]
                                            if partneib_short == '':
                                                partneib_short = "w"

                                            inter = "{}_{}".format(parthost_short, partneib_short)
                                            key_LJ = ('LJ', parthost+str(resid), partneib+str(neib) )
                                            key_C  = ('Coul', parthost+str(resid), partneib+str(neib) )
                                            rowindex_LJ = res_to_rowindex[ key_LJ ]
                                            rowindex_C  = res_to_rowindex[ key_C ]

                                            ### Get energies from dict energyline_cols ###
                                            try:
                                                vdw  = energyline_cols[ rowindex_LJ ]
                                                coul = energyline_cols[ rowindex_C ]
                                            except KeyError:
                                                LOGGER.warning("Data not found for %s - %s",
                                                    parthost+str(resid), partneib+str(neib) )
                                                if resid not in missing_energydata:
                                                    missing_energydata.append(resid)

                                            Etot = float(vdw) + float(coul)
                                            outp_line = outp_fstr.format(
                                                time, resid, neib, inter, vdw, coul, Etot)
                                            print(outp_line, file=energyoutput)

                                elif host_is_pl and not neib_is_pl:

                                    for parthost in energy.molparts:

                                        ### remove "resid_" from parthost string ###
                                        parthost_short = parthost[7:]
                                        if parthost_short == '':
                                            parthost_short = "w"
                                        partneib       = "resid_"
                                        partneib_short = "w"

                                        inter = "{}_{}".format(parthost_short, partneib_short)
                                        key_LJ = ('LJ', parthost+str(resid), partneib+str(neib) )
                                        key_C  = ('Coul', parthost+str(resid), partneib+str(neib) )
                                        rowindex_LJ = res_to_rowindex[ key_LJ ]
                                        rowindex_C  = res_to_rowindex[ key_C ]

                                        # Get energies from dict energyline_cols
                                        try:
                                            vdw  = energyline_cols[ rowindex_LJ ]
                                            coul = energyline_cols[ rowindex_C ]
                                        except KeyError:
                                            LOGGER.warning("Data not found for %s - %s",
                                                parthost+str(resid), partneib+str(neib) )
                                            if resid not in missing_energydata:
                                                missing_energydata.append(resid)

                                        Etot = float(vdw) + float(coul)
                                        outp_line = outp_fstr.format(
                                            time, resid, neib, inter, vdw, coul, Etot)
                                        print(outp_line, file=energyoutput)

                                elif not host_is_pl and neib_is_pl:

                                    parthost       = "resid_"
                                    parthost_short = "w"

                                    for partneib in energy.molparts:

                                        partneib_short = parthost[7:]
                                        if partneib_short == '':
                                            partneib_short = "w"

                                        inter = "{}_{}".format(parthost_short, partneib_short)
                                        key_LJ = ('LJ', parthost+str(resid), partneib+str(neib) )
                                        key_C  = ('Coul', parthost+str(resid), partneib+str(neib) )
                                        rowindex_LJ = res_to_rowindex[ key_LJ ]
                                        rowindex_C  = res_to_rowindex[ key_C ]

                                        # Get energies from dict energyline_cols
                                        try:
                                            vdw  = energyline_cols[ rowindex_LJ ]
                                            coul = energyline_cols[ rowindex_C ]
                                        except KeyError:
                                            LOGGER.warning("Data not found for %s - %s",
                                                parthost+str(resid), partneib+str(neib) )
                                            if resid not in missing_energydata:
                                                missing_energydata.append(resid)

                                        Etot = float(vdw) + float(coul)
                                        outp_line = outp_fstr.format(
                                            time, resid, neib, inter, vdw, coul, Etot)
                                        print(outp_line, file=energyoutput)

                                elif not host_is_pl and  not neib_is_pl:

                                    parthost       = "resid_"
                                    parthost_short = "w"
                                    partneib       = "resid_"
                                    partneib_short = "w"

                                    key_LJ = ('LJ', parthost+str(resid), partneib+str(neib) )
                                    key_C  = ('Coul', parthost+str(resid), partneib+str(neib) )
                                    rowindex_LJ = res_to_rowindex[ key_LJ ]
                                    rowindex_C  = res_to_rowindex[ key_C ]

                                    inter = "{}_{}".format(parthost_short, partneib_short)

                                    # Get energies from dict energyline_cols
                                    try:
                                        vdw  = energyline_cols[ rowindex_LJ ]
                                        coul = energyline_cols[ rowindex_C ]
                                    except KeyError:
                                        LOGGER.warning("Data not found for %s - %s",
                                            parthost+str(resid), partneib+str(neib) )
                                        if resid not in missing_energydata:
                                            missing_energydata.append(resid)

                                    Etot = float(vdw) + float(coul)
                                    outp_line = outp_fstr.format(
                                        time, resid, neib, inter, vdw, coul, Etot)
                                    print(outp_line, file=energyoutput)

    if missing_energydata:
        LOGGER.warning("Missing energydata: %s", missing_energydata)
        raise RuntimeError("There were inconsistencies in the data. See log files for further information.")

    LOGGER.info("File %s written successfully", energy.all_energies)


def create_lipid_water_interaction_file(energy, outputfilename="water_interaction.dat"):
    ''' Create a file with entries of
        interaction of resid at time to solvent
        <time> <resid> <resname> <Etot> <Evdw> <Ecoul>
    '''
    energyoutput = open(outputfilename, "w")
    print('{: <10}{: <10}{: <10}{: <20}{: <20}{: <20}'\
        .format("time", "resid", "resname", "Etot", "Evdw", "Ecoul"), file=energyoutput)
    for resid in energy.lipid_resids:
        resname = energy.convert.resid_to_resname[resid]
        xvgfilename = "{}/xvgtables/energies_residue{}_0.xvg"\
                .format(energy.path.energy, str(resid))

        with open(xvgfilename,"r") as xvgfile:
            res_to_rowindex = {}

            for energyline in xvgfile:
                energyline_cols = energyline.split()

                if '@ s' in energyline:

                    rowindex  = int(energyline_cols[1][1:])+1 # time is at row 0 !
                    host = energyline_cols[3].split("-")[1].split("resid_")[1]
                    energytype = energyline_cols[3].split("-")[0][1:]
                    neib = energyline_cols[3].split("-")[2][:-1]
                    if neib == "solv":
                        LOGGER.debug("Hostid: %s", host)
                        LOGGER.debug("Adding to dict: Etype %s, host %s", energytype, host)
                        res_to_rowindex[(energytype, host)] = rowindex

                elif '@' not in energyline and '#' not in energyline:

                    time = float(energyline_cols[0])
                    if time % energy.dt != 0:
                        continue

                    vdw  = float( energyline_cols[ res_to_rowindex[ ( 'LJ', str(resid) ) ] ] )
                    coul = float( energyline_cols[ res_to_rowindex[ ( 'Coul', str(resid) ) ] ] )
                    Etot = vdw + coul
                    outpline = '{: <10}{: <10}{: <10}{: <20.5f}{: <20.5f}{: <20.5f}'\
                            .format(time, resid, resname, Etot, vdw, coul,)
                    LOGGER.debug("%s", outpline)
                    print(outpline, file=energyoutput)
    energyoutput.close()


def create_indexfile(sysinfo, resindex_filename="resindex_all.ndx"):
    '''
        Creates an indexfile that is used for the energy calculation
        containing indices for each molecule in the system (named resid_X) along with an index
        group of the whole system.

        The selection is made using gmx select tool (see respective manual for further information)
    '''

    sysinfo.path.create_folders()

    with open(resindex_filename, "w") as resindex_all:
        LOGGER.info("\n_____Creating index file____\n")
        for mol in sysinfo.lipid_resids:
            resname = sysinfo.convert.resid_to_resname[mol]
            LOGGER.info("Working on residue %s", mol)

            ### Create selection file ###
            selectionfilename = sysinfo.path.tmp + '/tmp_selectionfile'
            with open(selectionfilename,"w") as sf:

                selectionlist = []
                selectionlist += [ 'resid_{0}=resid {0} and resname {1};\n'\
                                .format(str(mol), resname,) ]

                if not (sysinfo.ff.is_sterol(resname) or sysinfo.ff.is_protein(resname)):
                    tailatoms = [i for j in sysinfo.ff.tail_atoms_of(resname) for i in j]
                    tailatoms_str = ' '.join( tailatoms )
                    headatoms_str = ' '.join( sysinfo.ff.head_atoms_of(resname) )
                    selectionlist += ['resid_{0}_h=resid {0} and resname {1} and name {2};\n'\
                                    .format(str(mol), resname, headatoms_str)]
                    selectionlist += ['resid_{0}_t=resid {0} and resname {1} and name {2};\n'\
                                    .format(str(mol), resname, tailatoms_str)]
                    selectionlist += ['resid_{}_h;\n'.format(str(mol)) ]
                    selectionlist += ['resid_{}_t;\n'.format(str(mol)) ]

                selectionlist += ['resid_{};\n'.format(str(mol))]
                selectionlist = ''.join( selectionlist )

                sf.write( selectionlist )

            ### Run gmx select ###
            outputindex = sysinfo.path.index + "/resid_" + str(mol) + ".ndx"
            gmx_select_arglist = ['-s', sysinfo.path.tpr, '-sf',
                                  selectionfilename, '-on', outputindex,
                                  ]
            out, err = exec_gromacs(GMXNAME, "select", gmx_select_arglist)
            write_log("gmx_select", out, err, path=sysinfo.path.log)

            ### Write output of last gmx_select command to resindex file ###
            with open(outputindex, "r") as output_index:
                filecontent = output_index.readlines()
                resindex_all.write( ''.join(filecontent)+'\n\n' )

        ### Add index group containing the whole system ###
        make_ndx_output = sysinfo.path.tmp + '/make_ndx_system.ndx'
        gmx_make_ndx_arglist = [ '-f', sysinfo.path.gro, '-o', make_ndx_output ]
        inp_str = b'keep 0\nq\n'
        out, err = exec_gromacs(GMXNAME, "make_ndx", gmx_make_ndx_arglist, inp_str)
        write_log("gmx_make_ndx", out, err, path=sysinfo.path.log)

        with open(make_ndx_output, "r") as output:
            filecontent = output.readlines()
            resindex_all.write(''.join(filecontent)+"\n\n")

    add_water_groups_to_index(sysinfo, add_grp_to=resindex_filename)
    add_leaflet_groups_to_index(sysinfo, add_grp_to=resindex_filename)

def add_water_groups_to_index(sysinfo, add_grp_to="resindex_all.ndx"):
    ''' Make index file using gmx select and append index group to <add_grp_to> '''
    selectionstr = 'solv=resname {}; solv;'.format(sysinfo.water_resname)
    outputsel = sysinfo.path.tmp + "/tmp_leaflet.ndx"
    cmd = [
        "-f", sysinfo.path.gro,
        "-s", sysinfo.path.tpr, "-select", selectionstr,
        "-on", outputsel
        ]
    out, err = exec_gromacs(GMXNAME, "select", cmd)
    write_log("gmx_select", out, err, sysinfo.path.log)
    with open(outputsel, "r") as selectionf, open(add_grp_to, "a") as ndxf:
        for line in selectionf:
            ndxf.write(line)

def add_leaflet_groups_to_index(sysinfo, add_grp_to="resindex_all.ndx"):
    ''' Make index file using gmx select and append index group to <add_grp_to> '''
    #leafdat = pd.read_table(sysinfo., sep="\s+")
    resid_list = [[], []]
    outputsel = sysinfo.path.tmp + "/tmp_leaflet.ndx"
    for resid in sysinfo.lipid_resids:
        resid_list[sysinfo.convert.resid_to_leaflet[resid]].append(resid)
    selectionstr = 'leaflet0=resname {0} and resid {1}; leaflet1=resname {0} and resid {2};'\
                    'leaflet0; leaflet1;'.format(
            ' '.join(sysinfo.molecules),
            ' '.join( [ str(i) for i in resid_list[0] ] ),
            ' '.join( [ str(i) for i in resid_list[1] ] ),
        )
    cmd = [
        "-f", sysinfo.path.gro,
        "-s", sysinfo.path.tpr, "-select", selectionstr,
        "-on", outputsel
        ]
    out, err = exec_gromacs(GMXNAME, "select", cmd)
    write_log("gmx_select", out, err, path=sysinfo.path.log)
    with open(outputsel, "r") as selectionf, open(add_grp_to, "a") as ndxf:
        for line in selectionf:
            ndxf.write(line)

