# DGMETA (Version 1): Dissolved Gas Modeling and Environmental Tracer Analysis Computer Program v

DGMETA (Dissolved Gas Modeling and Environmental Tracer Analysis) is a Microsoft Excel-based computer program that is used for modeling
air-water equilibrium conditions from measurements of dissolved gases and for computing concentrations of environmental tracers that 
rely on air-water equilibrium model results. DGMETA can solve for the temperature, salinity, excess air, fractionation of gases, 
or pressure/elevation of water when it is equilibrated with the atmosphere. Models are calibrated inversely using one or more measurements 
of dissolved gases such as helium, neon, argon, krypton, xenon, and nitrogen. Excess nitrogen gas, originating from denitrification or 
other sources, also can be included as a fitted parameter or as a separate calculation from the dissolved gas modeling results. 
DGMETA uses the air-water equilibrium models to separate measured concentrations of gases and isotopes of gases into components that are 
used for tracing water in the environment. DGMETA calculates atmospheric dry-air mole fractions (mixing ratios) for transient atmospheric 
gas tracers such as chlorofluorocarbons, sulfur hexafluoride, and bromotrifluoromethane (Halon-1301); and concentrations of tritiogenic helium-3 
and radiogenic helium-4, which accumulate from the decay of tritium in water and the decay of uranium and thorium in rocks, respectively. 

## Installation and Requirements

Download the zipped DGMETA file and unzip the contents to a user directory. 
 - The zip file contains an empty DGMETA excel file, this read me, the publication document of the program, and a folder containing 
 example DGMETA files described in the program documentation.
 - The files in the example folder contain the data and results used in the program documentation.
 - To use the program, Microsoft Excel version 2007 or later must be installed on the machine
 - The program is written in Visual Basic for Applications (VBA) and macros must be enabled to use program.
 - In many cases the program may not function when first opened because it is not a trusted source. In this case, it is necessary to 
 save the program using a new name (DGMETA_Central_Valley for example), close it, then reopen it to enable the functionality.

## Disclaimer
This software has been approved for release by the U.S. Geological Survey (USGS). Although the software has been subjected to rigorous review, 
the USGS reserves the right to update the software as needed pursuant to further analysis and review. No warranty, expressed or implied, 
is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute 
any such warranty. Furthermore, the software is released on condition that neither the USGS nor the U.S. Government shall be held liable for any 
damages resulting from its authorized or unauthorized use.

## Authors and acknowledgment
Bryant Jurgens, Jon Karl Böhlke, Karl Haase, Eurybiades Busenberg, Andrew G. Hunt, and Jeffrey A. Hansen

## License
The DGMETA file is open source software and all code is contained within the program and can be modified freely. The program is made available by
U.S. Geological Survey (USGS) to be used in the public interest and in the advancement of science.  You may, without any fee or cost, use, copy, 
modify, or distribute this software, and any derivative works thereof, and its supporting documentation, subject to the following restrictions and 
understandings.

If you distribute copies or modifications of the software and related material, make sure the recipients receive a copy of this notice and receive 
or can get a copy of the original distribution.  If the software and (or) related material are modified and distributed, it must be made clear 
that the recipients do not have the original and they must be informed of the extent of the modifications.  For example, modified files must 
include a prominent notice stating the modifications made, the author of the modifications, and the date the modifications were made. This restriction 
is necessary to guard against problems introduced in the software by others.

## Project status
DGMETA is actively managed. Users are encouraged to report any bugs, issues, or questions to Bryant Jurgens at bjurgens@usgs.gov.
Please see the Version History for known issues of past versions.

## IPDS number
IP-156811