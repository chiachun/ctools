/***************************************************************************
 *                   ctbkgcube - CTA background cube tool                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Chia-Chun Lu                                     *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file ctbkgcube.cpp
 * @brief CTA background cube tool implementation
 * @author Chia-Chun Lu
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctbkgcube.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_GET_PARAMETERS                        "ctbkgcube::get_parameters()"
#define G_GET_OBS                                      "ctbkgcube::get_obs()"
#define G_GET_EBOUNDS                              "ctbkgcube::get_ebounds()"
#define G_SET_FROM_CNTMAP          "ctbkgcube::set_from_cntmap(std::string&)"

/* __ Debug definitions __________________________________________________ */

/* __ Coding definitions _________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
ctbkgcube::ctbkgcube(void) : GApplication(CTBKGCUBE_NAME, CTBKGCUBE_VERSION)
{
    // Initialise members
    init_members();

    // Write header into logger
    log_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Observations constructor
 *
 * param[in] obs Observation container.
 *
 * This method creates an instance of the class by copying an existing
 * observations container.
 ***************************************************************************/
ctbkgcube::ctbkgcube(const GObservations& obs) :
           GApplication(CTBKGCUBE_NAME, CTBKGCUBE_VERSION)
{
    // Initialise members
    init_members();

    // Set observations
    m_obs = obs;

    // Write header into logger
    log_header();

    // Return
    return;
}



/***********************************************************************//**
 * @brief Command line constructor
 *
 * @param[in] argc Number of arguments in command line.
 * @param[in] argv Array of command line arguments.
 ***************************************************************************/
ctbkgcube::ctbkgcube(int argc, char *argv[]) :
           GApplication(CTBKGCUBE_NAME, CTBKGCUBE_VERSION, argc, argv)
{
    // Initialise members
    init_members();

    // Write header into logger
    log_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] app ctbkgcube application.
 ***************************************************************************/
ctbkgcube::ctbkgcube(const ctbkgcube& app) : GApplication(app)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(app);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
ctbkgcube::~ctbkgcube(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] app ctbkgcube application.
 * @return Returns ctbkgcube application.
 ***************************************************************************/
ctbkgcube& ctbkgcube::operator=(const ctbkgcube& app)
{
    // Execute only if object is not identical
    if (this != &app) {

        // Copy base class members
        this->GApplication::operator=(app);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(app);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 ***************************************************************************/
void ctbkgcube::clear(void)
{
    // Free members
    free_members();
    this->GApplication::free_members();

    // Initialise members
    this->GApplication::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Execute application
 *
 * This is the main execution method of the ctbkgcube class. It is invoked
 * when the executable is called from command line. The method generates
 * the background cube and saves the result.
 ***************************************************************************/
void ctbkgcube::execute(void)
{
    // Signal that some parameters should be read ahead
    m_read_ahead = true;

    // Create the background cube
    run();

    // Save the background cube into the output FITS file
    save();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Generate the background cube(s).
 *
 * This method reads the task parameters from the parfile, sets up the
 * observation container, loops over all CTA observations in the container
 * and generates a background cube from the CTA observations.
 ***************************************************************************/
void ctbkgcube::run(void)
{
    // If we're in debug mode then all output is also dumped on the screen
    if (logDebug()) {
        log.cout(true);
    }

    // Get task parameters
    get_parameters();

    // Write parameters into logger
    if (logTerse()) {
        log_parameters();
        log << std::endl;
    }

    // Write observation(s) into logger
    if (logTerse()) {
        log << std::endl;
        if (m_obs.size() > 1) {
            log.header1("Observations");
        }
        else {
            log.header1("Observation");
        }
        log << m_obs << std::endl;
    }

    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Prepare model");
    }

    // Copy models from observation container
    m_bkgmdl = m_obs.models();

    // Remove all model that are not CTA background models from the
    // container
    int num = m_bkgmdl.size();
    for (int i = num-1; i >= 0; --i) {

        // Flag removal
        bool remove = true;

        // Do we have a CTA, HESS or VERITAS specific model?
        if (m_bkgmdl[i]->is_valid("CTA", "") ||
            m_bkgmdl[i]->is_valid("HESS", "") ||
            m_bkgmdl[i]->is_valid("VERITAS", "")) {

            // Do we have a background model?
            if (dynamic_cast<GModelData*>(m_bkgmdl[i]) != NULL) {

                // ... then keep model
                remove = false;

            } // endif: had a background model
        } // endif: had a CTA, HESS or VERITAS model

        // Log results
        if (logTerse()) {
            if (remove) {
                log << gammalib::parformat("Remove model");
            }
            else {
                log << gammalib::parformat("Keep model");
            }
            log << m_bkgmdl[i]->name();
            log << " ";
            log << m_bkgmdl[i]->type();
            log << "(";
            log << m_bkgmdl[i]->instruments();
            log << ")";
            log << std::endl;
        }

        // Remove model if requested
        if (remove) {
            m_bkgmdl.remove(i);
        }

    } // endfor: looped over all background models

    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Generate background cube");
    }

    // Initialise observation counter
    int n_observations = 0;

    // Loop over all observations in the container
    for (int i = 0; i < m_obs.size(); ++i) {
      
        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

        // Continue only if observation is a CTA observation
        if (obs != NULL) {

            // Write header for observation
            if (logTerse()) {
                if (obs->name().length() > 1) {
                    log.header3("Observation "+obs->name());
                }
                else {
                    log.header3("Observation "+gammalib::str(i));
                }
            }

            // Increment number of observations
            n_observations++;

            // Fill the cube
            fill_cube(obs);
	    	  
        } // endif: CTA observation found

    } // endfor: looped over observations

    // Scale the cube by ontime to get the unit right
    for (int i = 0; i < m_bkgcube.size(); ++i) {

            // Get event bin
            GCTAEventBin* bin = m_bkgcube[i];
	    
	    // Divide the content by ontime
            bin->counts(bin->counts()/m_ontime);
    }

    // Log results
    if (logTerse()) {
        log << gammalib::parformat("Number of observations");
        log << n_observations << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Generate background cube
 *
 * @param[in] obs Pointer to CTA observation.
 *
 * Fills
 ***************************************************************************/
void ctbkgcube::fill_cube(GCTAObservation* obs)
{
    // Continue only if observation pointer is valid
    if (obs != NULL) {

        // Initialise statistics
        double sum = 0.0;

	// Add ontime to total ontime
        m_ontime = m_ontime + obs->ontime();

        // Set GTI of actual observations as the GTI of the event cube
        m_bkgcube.gti(obs->events()->gti());

        // Loop over all bins in background cube
        for (int i = 0; i < m_bkgcube.size(); ++i) {

            // Get event bin
            GCTAEventBin* bin = m_bkgcube[i];
            
            // Compute model value for event bin
            double model = 
                   m_bkgmdl.eval(*(const_cast<const GCTAEventBin*>(bin)), *obs) *
	           bin->size();

            // Add existing number of counts
            model += bin->counts();

            // Store cumulated value
            bin->counts(model);

            // Sum all events
            sum += model;
        }

        // Log results
        if (logTerse()) {
            log << gammalib::parformat("Background events in cube");
            log << sum << std::endl;
        }

    } // endif: observation pointer was not valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save background cube
 *
 * Save the background cube into the file specified by the outfile parameter.
 ***************************************************************************/
void ctbkgcube::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Save background cube");
    }

    // Create energies container from energy boundaries
    GEnergies energies;
    for (int i = 0; i < m_ebounds.size(); ++i) {
        energies.append(m_ebounds.elogmean(i));
    }

    // Get output filename
    std::string filename = (*this)["outfile"].filename();

    // Create empty FITS file
    GFits fits;

    // Write background cube
    m_bkgcube.map().write(fits);

    // Write energies
    energies.write(fits);
    
    // Save FITS file
    fits.saveto(filename, clobber());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * @exception GException::invalid_value
 *            No background model definition XML file specified.
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user. The parameters are read in the correct order.
 ***************************************************************************/
void ctbkgcube::get_parameters(void)
{
    // If we do not have any observations in the container then load them
    // using the "infile" parameter
    if (m_obs.size() == 0) {
        get_obs();
    }

    // If there are no models associated with the observations then load now
    // the model definition from the XML file
    if (m_obs.models().size() == 0) {
        if ((*this)["bkgmdl"].is_undefined()) {
            std::string msg = "No background model definition XML file "
                              "specified. Please set the \"bkgmdl\" parameter "
                              "to the background model definition XML "
                              "filename.";
            throw GException::invalid_value(G_GET_PARAMETERS, msg);
        }
        std::string bkgmdl = (*this)["bkgmdl"].filename();
        GModels     models(bkgmdl);
        m_obs.models(models);
    }

    // If no counts map is specified then setup the background cube from
    // the user parameters
    std::string cntmap = (*this)["cntmap"].filename();
    if ((gammalib::toupper(cntmap) == "NONE") ||
        (gammalib::strip_whitespace(cntmap) == "")) {
    
        // Get user parameters for counts map definition
        std::string proj     = (*this)["proj"].string();
        std::string coordsys = (*this)["coordsys"].string();
        double      xref     = (*this)["xref"].real();
        double      yref     = (*this)["yref"].real();
        double      binsz    = (*this)["binsz"].real();
        int         nxpix    = (*this)["nxpix"].integer();
        int         nypix    = (*this)["nypix"].integer();

        // Get energy definition
        get_ebounds();

        // Define skymap for background cube
        GSkymap map(proj, coordsys, xref, yref,
                    -binsz, binsz, nxpix, nypix,
                    m_ebounds.size());

        // Allocate background cube using a dummy GTI
        GGti gti;
        gti.append(GTime(0.0), GTime(1.0)); // Dummy GTI
        m_bkgcube = GCTAEventCube(map, m_ebounds, gti);

    }

    // ... otherwise setup the background cube from the counts map
    else {
    
        // Set background cube from counts map
        set_from_cntmap(cntmap);
    
    }

    // Read output filename (if needed)
    if (m_read_ahead) {
        std::string filename = (*this)["outfile"].filename();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get observation definition
 *
 * @exception GException::invalid_value
 *            No event list, counts cube or observation definition file
 *            specified.
 *
 * Get observation definition from the user parameters.
 ***************************************************************************/
void ctbkgcube::get_obs(void)
{
    // Get input filename
    std::string filename = (*this)["infile"].filename();

    // Check whether infile parameter is undefined
    if (gammalib::toupper(filename) == "NONE") {
        std::string msg = "No event list, counts cube or observation definition "
                          "file specified. Please set the \"infile\" parameter "
                          "to either an event list filename, a counts cube "
                          "filename or an observation definition filename.";
        throw GException::invalid_value(G_GET_OBS, msg);
    }

    // Try first to open as FITS file
    try {

        // Allocate CTA observation
        GCTAObservation obs;
        
        // Load input file in CTA observation
        obs.load(filename);

        // Append CTA observation to container
        m_obs.append(obs);
            
    }
        
    // ... otherwise try to open as XML file
    catch (GException::fits_open_error &e) {

        // Load observations from XML file
        m_obs.load(filename);

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get energy boundaries from parameters
 *
 * @exception GException::invalid_value
 *            Energy bin definition file not defined or invalid extension
 *            name encountered.
 *
 * Get the energy boundaries from the user parameters.
 ***************************************************************************/
void ctbkgcube::get_ebounds(void)
{
    // Determine the energy binning alogrithm
    std::string ebinalg = (*this)["ebinalg"].string();

    // If we have the binning given by a file then try to get the boundaries
    // from that file
    if (ebinalg == "FILE") {

        // Get filename
        std::string filename = (*this)["ebinfile"].filename();

        // Check whether infile parameter is undefined
        if (gammalib::toupper(filename) == "NONE") {
            std::string msg = "No energy bin definition file specified. "
                              "Please set the \"ebinfile\" parameter to a "
                              "file containing the energy bin definition.";
            throw GException::invalid_value(G_GET_EBOUNDS, msg);
        }

        // Open fits file to check which extension is given
        GFits file(filename);

        // Check first for EBOUNDS extension
        if (file.contains("EBOUNDS")) {
            file.close();
            m_ebounds.load(filename, "EBOUNDS");
        }

        // ... then check for ENERGYBINS extension
        else if (file.contains("ENERGYBINS")) {
            file.close();
            m_ebounds.load(filename, "ENERGYBINS");
        }

        // ... otherwise throw an exception
        else {
            file.close();
            std::string msg = "No extension with name \"EBOUNDS\" or"
                              " \"ENERGYBINS\" found in FITS file"
                              " \""+filename+"\".\n"
                              "An \"EBOUNDS\" or \"ENERGYBINS\" extension"
                              " is required if the parameter \"ebinalg\""
                              " is set to \"FILE\".";
            throw GException::invalid_value(G_GET_EBOUNDS, msg);
        }
    }
    
    // ... otherwise read emin, emax and nebins
    else {

        // Get the relevant parameters
    	double emin     = (*this)["emin"].real();
    	double emax     = (*this)["emax"].real();
    	int    enumbins = (*this)["enumbins"].integer();
        bool   log      = ((*this)["ebinalg"].string() == "LIN") ? false : true;

        // Create energy boundaries
        m_ebounds = GEbounds(enumbins,
                             GEnergy(emin, "TeV"),
                             GEnergy(emax, "TeV"),
                             log);

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set background cube definition from counts map
 *
 * @exception GException::invalid_value
 *            Invalid counts map projection or invalid events encountered.
 *
 * Set background cube definition from counts map.
 ***************************************************************************/
void ctbkgcube::set_from_cntmap(const std::string& filename)
{
    // Allocate CTA observation
    GCTAObservation obs;
        
    // Load counts map in CTA observation
    obs.load(filename);

    // Set background cube from counts map
    const GCTAEventCube* cube = dynamic_cast<const GCTAEventCube*>(obs.events());

    // Continue only if cube is valid
    if (cube != NULL) {

        // Get sky map projection
        const GWcs* wcs = dynamic_cast<const GWcs*>(cube->map().projection());
        
        // Continue only if projection is valid
        if (wcs != NULL) {
            
            // Get user parameters for counts map definition
            std::string proj     = wcs->code();
            std::string coordsys = wcs->coordsys();
            double      xref     = wcs->crval(0);
            double      yref     = wcs->crval(1);
            double      dx       = wcs->cdelt(0);
            double      dy       = wcs->cdelt(1);
            int         nx       = cube->map().nx();
            int         ny       = cube->map().ny();

            // Get energy definition
            m_ebounds = cube->ebounds();

            // Define skymap for background cube
            GSkymap map(proj, coordsys, xref, yref,
                        dx, dy, nx, ny,
                        m_ebounds.size());

            // Allocate background cube
            m_bkgcube = GCTAEventCube(map, m_ebounds, obs.events()->gti());

        } // endif: WCS projection was valid

        // ... projection is not of WCS type
        else {
            std::string msg = "Counts map project is not of WCS type.";
            throw GException::invalid_value(G_SET_FROM_CNTMAP, msg);
        }

    } // endif: observation contained an events cube

    // ... there is not events cube
    else {
        std::string msg = "No events cube found in file \""
                          ""+filename+"\".";
        throw GException::invalid_value(G_SET_FROM_CNTMAP, msg);
    }

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void ctbkgcube::init_members(void)
{
    // Initialise members
    m_read_ahead = false;
    m_obs.clear();
    m_bkgcube.clear();
    m_bkgmdl.clear();
    m_ebounds.clear();

    // Set logger properties
    log.date(true);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Application.
 ***************************************************************************/
void ctbkgcube::copy_members(const ctbkgcube& app)
{
    // Copy members
    m_read_ahead = app.m_read_ahead;
    m_obs        = app.m_obs;
    m_bkgcube    = app.m_bkgcube;
    m_bkgmdl     = app.m_bkgmdl;
    m_ebounds    = app.m_ebounds;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctbkgcube::free_members(void)
{
    // Write separator into logger
    if (logTerse()) {
        log << std::endl;
    }

    // Return
    return;
}
