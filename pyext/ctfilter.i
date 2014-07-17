/***************************************************************************
 *                    ctselect - CTA data selection tool                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @file ctfilter.i
 * @brief CTA data selection tool Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctfilter.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class ctfilter
 *
 * @brief CTA data selection tool Python interface
 ***************************************************************************/
class ctfilter : public GApplication  {
public:
    // Constructors and destructors
    ctfilter(void);
    explicit ctfilter(GObservations obs);
    ctfilter(int argc, char *argv[]);
    ctfilter(const ctfilter& app);
    virtual ~ctfilter(void);

    // Methods
    void           clear(void);
    void           execute(void);
    void           run(void);
    void           save(void);
    GObservations& obs(void) { return m_obs; }
    void           get_parameters(void);
    void           apply_filter(GCTAObservation* obs, const std::string& filename);
};


/***********************************************************************//**
 * @brief CTA data selection tool Python extension
 ***************************************************************************/
%extend ctfilter {
    ctfilter copy() {
        return (*self);
    }
}
