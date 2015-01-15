To compile the code, [ROOT](http://root.cern.ch/ "ROOT") has to be installed. For building use the following command:

``g++ -std=gnu++11 -O3 -o simulate simulate.cpp `root-config --cflags --glibs` -lSpectrum``

Usage: ``./simulate particle_name [output_file]``

Possible particles are:
- photon
- proton
- electron (positron)
- muon (antimuon)
