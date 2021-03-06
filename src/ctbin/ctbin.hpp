/***************************************************************************
 *                      ctbin - CTA data binning tool                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @file ctbin.hpp
 * @brief CTA data binning tool definition
 * @author Juergen Knoedlseder
 */

#ifndef CTBIN_HPP
#define CTBIN_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCTALib.hpp"

/* __Definitions _________________________________________________________ */
#define CTBIN_NAME    "ctbin"
#define CTBIN_VERSION "00-05-00"


/***********************************************************************//**
 * @class ctbin
 *
 * @brief CTA data binning tool interface defintion
 *
 * This class bins CTA event list(s) into a single counts cube. The class
 * can operate on predefined observation containers, on individual event list
 * FITS files, and on observation definition XML files.
 *
 * If multiple event lists are specified in the observation container or the
 * XML definition file, the class will merge these events into a single
 * counts cube.
 *
 * Results are stored in an observation container that can be written to disk
 * in form of a single FITS file.
 *
 * WARNING: Note that the pointing direction of the counts cube will be set
 * to the skymap centre used for the counts cube definition. If usepnt=yes
 * is used, the pointing direction will be extracted from the first
 * observation encountered in the list. Ultimately, pointing direction
 * information should not be removed from the counts cube and exposure and
 * PSF cubes should be used for response computation. This is however not
 * yet fully implemented.
 ***************************************************************************/
class ctbin : public GApplication  {
public:
    // Constructors and destructors
    ctbin(void);
    explicit ctbin(GObservations obs);
    ctbin(int argc, char *argv[]);
    ctbin(const ctbin& app);
    virtual ~ctbin(void);

    // Operators
    ctbin& operator=(const ctbin& app);

    // Methods
    void                 clear(void);
    void                 execute(void);
    void                 run(void);
    void                 save(void);
    const GObservations& obs(void) const;
    GCTAEventCube        cube(void) const;
    void                 get_parameters(void);
    void                 init_cube(void);
    void                 fill_cube(GCTAObservation* obs);
    void                 obs_cube(void);

protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const ctbin& app);
    void           free_members(void);
    void           get_ebounds(void);
    std::string    set_outfile_name(const int index) const;
    void           save_fits(void);
    void           save_xml(void);
    void           save_counts_map(const GCTAObservation* obs,
                                   const std::string&     outfile) const;

    // Obsolete methods (to be reomved later)
    void           bin_events(GCTAObservation* obs);

    // User parameters
    std::string              m_evfile;     //!< Input event list or XML file
    std::string              m_outfile;    //!< Output counts map or XML file
    std::string              m_prefix;     //!< Prefix for multiple counts maps
    bool                     m_usepnt;     //!< Use pointing instead of xref/yref parameters
    double                   m_emin;       //!< Lower energy
    double                   m_emax;       //!< Upper energy
    int                      m_enumbins;   //!< Number of energy bins
    std::string              m_proj;       //!< WCS projection
    std::string              m_coordsys;   //!< Coordinate system
    std::string              m_ebinalg;    //!< Algorithm for energy binning
    std::string              m_ebinfile;   //!< FITS-file containing energy binning
    double                   m_xref;       //!< Longitude reference coordinate
    double                   m_yref;       //!< Latitude reference coordinate
    double                   m_binsz;      //!< Pixel size
    int                      m_nxpix;      //!< Number of pixels in longitude
    int                      m_nypix;      //!< Number of pixels in latitude

    // Protected members
    GObservations            m_obs;        //!< Observation container
    std::vector<std::string> m_infiles;    //!< Input event filenames
    bool                     m_use_xml;    //!< Use XML file instead of FITS file
    bool                     m_read_ahead; //!< Read ahead parameters
    GSkymap                  m_cube;       //!< Event cube
    GEbounds                 m_ebounds;    //!< Energy boundaries
    GGti                     m_gti;        //!< Good time intervals
    double                   m_ontime;     //!< Total ontime
    double                   m_livetime;   //!< Total livetime
};

/***********************************************************************//**
 * @brief Return observation container
 *
 * @return Reference to observation container
 ***************************************************************************/
inline
const GObservations& ctbin::obs(void) const
{
    return m_obs;
}


/***********************************************************************//**
 * @brief Return event cube
 *
 * @return Event cube
 *
 * Returns a CTA event cube built from the ctbin member information.
 ***************************************************************************/
inline
GCTAEventCube ctbin::cube(void) const
{
    // Build event cube
    GCTAEventCube cube(m_cube, m_ebounds, m_gti);
    
    // Return cube
    return cube;
}

#endif /* CTBIN_HPP */
