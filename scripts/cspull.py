#! /usr/bin/env python
# ==========================================================================
# This script generates the pull distribution for all model parameters.
#
# Copyright (C) 2011-2014 Juergen Knoedlseder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==========================================================================
import gammalib
from ctools import obsutils
import sys
import csv
import math


# ============ #
# cspull class #
# ============ #
class cspull(gammalib.GApplication):
    """
    This class implements the pull distribution generation script. It derives
    from the GammaLib::GApplication class which provides support for parameter
    files, command line arguments, and logging. In that way the Python
    script behaves just as a regular ctool.
    """
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self.name    = "cspull"
        self.version = "0.3.0"
        
        # Initialise some members
        self.obs      = None
        self.model    = None
        self.m_srcmdl = None
        
        # Make sure that parfile exists
        file = self.parfile()

        # Initialise application
        if len(argv) == 0:
            gammalib.GApplication.__init__(self, self.name, self.version)
        elif len(argv) ==1:
            gammalib.GApplication.__init__(self, self.name, self.version, *argv)
        else:
            raise TypeError("Invalid number of arguments given.")

        # Set logger properties
        self.log_header()
        self.log.date(True)

        # Return
        return
    
    def __del__(self):
        """
        Destructor.
        """
        #  Write separator into logger
        if self.logTerse():
            self.log("\n")
        
        # Return
        return

    def parfile(self):
        """
        Check if parfile exists. If parfile does not exist then create a
        default parfile. This kluge avoids shipping the cscript with a parfile.
        """
        # Set parfile name
        parfile = self.name+".par"
        
        try:
            pars = gammalib.GApplicationPars(parfile)
        except:
            # Signal if parfile was not found
            sys.stdout.write("Parfile "+parfile+" not found. Create default parfile.\n")
            
            # Create default parfile
            pars = gammalib.GApplicationPars()
            pars.append(gammalib.GApplicationPar("srcmdl","f","a","$CTOOLS/share/models/crab.xml","","","Source model"))
            pars.append(gammalib.GApplicationPar("outfile","f","a","pull.dat","","","Output file name"))
            pars.append(gammalib.GApplicationPar("ntrials","i","a","10","","","Number of trials"))
            pars.append(gammalib.GApplicationPar("caldb","s","a","dummy","","","Calibration database"))
            pars.append(gammalib.GApplicationPar("irf","s","a","cta_dummy_irf","","","Instrument response function"))
            pars.append(gammalib.GApplicationPar("edisp","b","h","no","","","Apply energy dispersion?"))
            pars.append(gammalib.GApplicationPar("ra","r","a","83.6331","0","360","RA of pointing (deg)"))
            pars.append(gammalib.GApplicationPar("dec","r","a","22.0145","-90","90","Dec of pointing (deg)"))
            pars.append(gammalib.GApplicationPar("emin","r","a","0.1","0.0","","Lower energy limit (TeV)"))
            pars.append(gammalib.GApplicationPar("emax","r","a","100.0","0.0","","Upper energy limit (TeV)"))
            pars.append(gammalib.GApplicationPar("enumbins","i","a","0","","","Number of energy bins (0=unbinned)"))
            pars.append(gammalib.GApplicationPar("duration","r","a","1800.0","","","Effective exposure time (s)"))
            pars.append(gammalib.GApplicationPar("deadc","r","h","0.95","","","Deadtime correction factor"))
            pars.append(gammalib.GApplicationPar("rad","r","h","5.0","","","Radius of ROI (deg)"))
            pars.append(gammalib.GApplicationPar("pattern","s","h","single","","","Observation pattern (single/four)"))
            pars.append(gammalib.GApplicationPar("offset","r","h","1.5","","","Observation pattern offset (deg)"))
            pars.append(gammalib.GApplicationPar("npix","i","h","200","","","Number of pixels for binned"))
            pars.append(gammalib.GApplicationPar("binsz","r","h","0.05","","","Pixel size for binned (deg/pixel)"))
            pars.append_standard()
            pars.save(parfile)
        
        # Return
        return
        
    def get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        # Get parameters
        if self.model == None:
            self.m_srcmdl = self["srcmdl"].filename()
        self.m_outfile  = self["outfile"].filename()
        self.m_ntrials  = self["ntrials"].integer()
        self.m_caldb    = self["caldb"].string()
        self.m_irf      = self["irf"].string()
        self.m_edisp    = self["edisp"].boolean()
        self.m_ra       = self["ra"].real()
        self.m_dec      = self["dec"].real()
        self.m_emin     = self["emin"].real()
        self.m_emax     = self["emax"].real()
        self.m_enumbins = self["enumbins"].integer()
        self.m_duration = self["duration"].real()
        self.m_deadc    = self["deadc"].real()
        self.m_rad      = self["rad"].real()
        self.m_pattern  = self["pattern"].string()
        self.m_offset   = self["offset"].real()
        self.m_npix     = self["npix"].integer()
        self.m_binsz    = self["binsz"].real()

        # Set some fixed parameters
        self.m_log   = False # Logging in client tools
        self.m_debug = False # Debugging in client tools
        
        # Setup observations
        self.obs = self.set_obs()
        
        # Load source model
        if self.m_srcmdl != None:
            self.model = gammalib.GModels(self.m_srcmdl)
        
        # Append source model to observation
        self.obs.models(self.model)

        # Return
        return
    
    def models(self, models):
        """
        Set model.
        """
        # Copy models
        self.model = models.copy()
    
        # Return
        return
        
    def execute(self):
        """
        Execute the script.
        """
        # Run the script
        self.run()
        
        # Return
        return

    def run(self):
        """
        Run the script.
        """
        # Switch screen logging on in debug mode
        if self.logDebug():
            self.log.cout(True)

        # Get parameters
        self.get_parameters()
        
        #  Write input parameters into logger
        if self.logTerse():
            self.log_parameters()
            self.log("\n")
        
        # Write observation into logger
        if self.logTerse():
            self.log("\n")
            self.log.header1("Observation")
            self.log(str(self.obs))
            self.log("\n")

        # Write header
        if self.logTerse():
            self.log("\n")
            self.log.header1("Generate pull distribution")

        # Loop over trials
        for seed in range(self.m_ntrials):
        
            # Make a trial
            result = self.trial(seed)
            
            # Write out result immediately
            if seed == 0:
                file = open(self.m_outfile, 'w')
                writer = csv.DictWriter(file, result['colnames'])
                writer.writerow(dict((_,_) for _ in result['colnames']))
            else:
                file = open(self.m_outfile, 'a')
            writer = csv.DictWriter(file, result['colnames'])
            writer.writerow(result['values'])
            file.close()
        
        # Return
        return
    
    def set_obs(self):
        """
        Returns an observation container with a set of CTA observations.
        
        Keywords:
        """
        # Setup observation definition list
        obsdeflist = obsutils.set_obs_patterns(self.m_pattern, \
                                               ra=self.m_ra, dec=self.m_dec, \
                                               offset=self.m_offset)
        
        # Create list of observations
        obs = obsutils.set_obs_list(obsdeflist, \
                                    tstart=0.0, duration=self.m_duration, \
                                    deadc=self.m_deadc, \
                                    emin=self.m_emin, emax=self.m_emax, \
                                    rad=self.m_rad, \
                                    irf=self.m_irf, caldb=self.m_caldb)
    
        # Return observation container
        return obs

    def trial(self, seed):
        """
        Create the pull for a single trial.
        
        Parameters:
         seed - Random number generator seed
        """
        # Write header
        if self.logExplicit():
            self.log.header2("Trial "+str(seed+1))

        # Simulate events
        obs = obsutils.sim(self.obs, \
                           nbins=self.m_enumbins, \
                           seed=seed, \
                           binsz=self.m_binsz, \
                           npix=self.m_npix, \
                           edisp=self.m_edisp, \
                           log=self.m_log, debug=self.m_debug)

        # Determine number of events in simulation
        nevents = 0.0
        for run in obs:
            nevents += run.events().number()

        # Write simulation results
        if self.logExplicit():
            self.log.header3("Simulation")
            self.log.parformat("Number of simulated events")
            self.log(nevents)
            self.log("\n")

        # Fit model
        like = obsutils.fit(obs, edisp=self.m_edisp, \
                            log=self.m_log, debug=self.m_debug)

        # Store results
        logL   = like.opt().value()
        npred  = like.obs().npred()
        models = like.obs().models()

        # Write result header
        if self.logExplicit():
            self.log.header3("Pulls")
        
        # Gather results
        colnames = []
        values   = {}
        colnames.append("LogL")
        colnames.append("Sim_Events")
        colnames.append("Npred_Events")
        values["LogL"]         = logL
        values["Sim_Events"]   = nevents
        values["Npred_Events"] = npred
        for i in range(models.size()):
            model      = models[i]
            model_name = model.name()
            for k in range(model.size()):
                par = model[k]
                if par.is_free():
                
                    # Set parameter name
                    name = model_name+"_"+par.name()
                    
                    # Append parameter, Pull_parameter and Unc_parameter
                    colnames.append(name)
                    colnames.append("Pull_"+name)
                    colnames.append("Unc_"+name)
                
                    # Compute pull
                    fitted_value = par.value()
                    real_value   = self.model[i][k].value()
                    error        = par.error()
                    if error != 0.0:
                        pull = (fitted_value - real_value) / error
                    else:
                        pull = 99.0
                        
                    # Store results
                    values[name] = fitted_value
                    values["Pull_"+name] = pull
                    values["Unc_"+name] = error

                    # Write result
                    if self.logExplicit():
                        self.log.parformat(name)
                        self.log(pull)
                        self.log(" (")
                        self.log(fitted_value)
                        self.log(" +/- ")
                        self.log(error)
                        self.log(")\n")
        
        # Bundle together results
        result = {'colnames': colnames, 'values': values}
        
        # Return
        return result


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Generates pull distribution for a source model.
    """
    # Create instance of application
    app = cspull(sys.argv)
    
    # Open logfile
    app.logFileOpen()
    
    # Execute application
    app.execute()
    
