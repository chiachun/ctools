/***************************************************************************
 *                     ctmodel - CTA counts model tool                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file ctmodel/main.cpp
 * @brief CTA counts model tool main code
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctmodel.hpp"


/***********************************************************************//**
 * @brief Main entry point of application
 *
 * @param[in] argc Number of command line arguments.
 * @param[in] argv Command line arguments.
 *
 * This is the main entry point of the gtbin application. It allocates a
 * ctmodel object and runs the application.
 ***************************************************************************/
int main (int argc, char *argv[])
{
    // Initialise return code
    int rc = 1;

    // Create instance of application
    ctmodel application(argc, argv);

    // Run application
    try {
        // Execute application
        application.execute();

        // Signal success
        rc = 0;
    }
    catch (std::exception &e) {

        // Extract error message
        std::string message = e.what();
        std::string signal  = "*** ERROR encounterted in the execution of"
                              " ctmodel. Run aborted ...";

        // Write error in logger
        application.log << signal  << std::endl;
        application.log << message << std::endl;

        // Write error on standard output
        std::cout << signal  << std::endl;
        std::cout << message << std::endl;

    } // endcatch: catched any application error

    // Return
    return rc;
}
