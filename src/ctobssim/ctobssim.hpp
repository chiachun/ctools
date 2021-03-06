/***************************************************************************
 *                ctobssim - CTA observation simulator tool                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2014 by Juergen Knoedlseder                         *
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
 * @file ctobssim.hpp
 * @brief CTA observation simulator tool definition
 * @author Juergen Knoedlseder
 */

#ifndef CTOBSSIM_HPP
#define CTOBSSIM_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCTALib.hpp"

/* __Definitions _________________________________________________________ */
#define CTOBSSIM_NAME    "ctobssim"
#define CTOBSSIM_VERSION "00-08-00"


/***********************************************************************//**
 * @class ctobssim
 *
 * @brief CTA observation simulator tool interface definition
 *
 * This class simulates CTA observation(s) using Monte Carlo sampling of the
 * source and background models. The class supports simulation of data of
 * multiple CTA observations in one shot. If multiple CTA observations are
 * processed and the save method is called, events FITS files will be written
 * for each observation.
 ***************************************************************************/
class ctobssim : public GApplication  {
public:
    // Constructors and destructors
    ctobssim(void);
    explicit ctobssim(GObservations obs);
    ctobssim(int argc, char *argv[]);
    ctobssim(const ctobssim& app);
    virtual ~ctobssim(void);

    // Operators
    ctobssim& operator=(const ctobssim& app);

    // Methods
    void           clear(void);
    void           execute(void);
    void           run(void);
    void           save(void);
    GObservations& obs(void) { return m_obs; }
    void           get_parameters(void);
    void           set_list(GCTAObservation* obs);
    void           simulate_source(GCTAObservation* obs, const GModels& models,GRan& ran, GLog* wrklog=NULL);
    void           simulate_background(GCTAObservation* obs, const GModels& models,GRan& ran, GLog* wrklog=NULL);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctobssim& app);
    void free_members(void);
    void save_fits(void);
    void save_xml(void);

    // User parameters
    std::string   m_infile;      //!< Input model
    std::string   m_outfile;     //!< Output events file
    std::string   m_prefix;      //!< Prefix for multiple event lists
    std::string   m_caldb;       //!< Calibration database repository
    std::string   m_irf;         //!< Instrument response function
    int           m_seed;        //!< Random number generator seed
    double        m_ra;          //!< RA of pointing direction
    double        m_dec;         //!< DEC of pointing direction
    double        m_rad;         //!< FOV radius (degrees)
    double        m_tmin;        //!< Start time (MET)
    double        m_tmax;        //!< Stop time (MET)
    double        m_emin;        //!< Lower energy (TeV)
    double        m_emax;        //!< Upper energy (TeV)
    bool          m_apply_edisp; //!< Apply energy dispersion?
    double        m_deadc;       //!< Average deadtime correction

    // Protected members
    double            m_area;        //!< Surface area for simulation (cm2)
    int               m_max_photons; //!< Maximum number of photons/slice
    std::vector<GRan> m_rans;        //!< Random number generators
    GObservations     m_obs;         //!< Observation container
    bool              m_use_xml;     //!< Use XML file instead of FITS file
    bool              m_read_ahead;  //!< Read ahead parameters
    GTimeReference    m_cta_ref;     //!< CTA time reference
    int               m_event_id;    //!< Event identifier
};

#endif /* CTOBSSIM_HPP */
