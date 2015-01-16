// compile with:
// g++ -std=gnu++11 -O3 -o simulate simulate.cpp `root-config --cflags --glibs` -lSpectrum

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string>

#include <TROOT.h>
#include <TFile.h>
#include <TNtupleD.h>
#include <TRandom3.h>

int main(int argc, char **argv)
{
	enum Particle_id {photon = 1, positron, electron, antimuon = 5, muon, proton = 14};
	const char* name;  // output file name
	Particle_id id_part;

	if (argc < 2) {
		printf("[WARNING] No arguments given, simulate photon\n");
		id_part = photon;
		name = "output.root";
	} else if (argc == 2 || argc == 3) {
		if (strstr(argv[1], "photon"))
			id_part = photon;
		else if (strstr(argv[1], "proton"))
			id_part = proton;
		else if (strstr(argv[1], "electron"))
			id_part = electron;
		else if (strstr(argv[1], "positron"))
			id_part = positron;
		else if (strstr(argv[1], "antimu"))
			id_part = antimuon;
		else if (strstr(argv[1], "muon"))
			id_part = muon;
		else {
			fprintf(stderr, "[ERROR] Unknown particle \"%s\", will exit\n", argv[1]);
			return 1;
		}

		if (argc == 3) {
			std::string tmp = argv[2];
			if (!strstr(argv[2], ".root"))
				tmp += ".root";
			name = tmp.c_str();
		} else
			name = "output.root";

		printf("[INFO] Simulate %s (id: %d), output will be written to %s\n", argv[1], id_part, name);
	} else {
		fprintf(stderr, "[ERROR] Too many arguments given, will exit\n");
		printf("   Usage: %s particle_name [output_file]\n", argv[0]);
		return 1;
	}

	static const double PI = 3.141592653589793;
	static const double MASS_PROTON = 938.272;
	static const double MASS_ELECTRON = .5109989;
	static const double MASS_MUON = 105.65837;

	/*
	 * Important!
	 * Values for the simulated ranges and statistics which should be changed!
	 */
	double start_energy = /*.8*/1.4, end_energy = 1.604, step_energy = .0005;  // GeV
	double start_theta = /*0.*/PI-.2, end_theta = PI, step_theta = PI/360.;  //radians, half a degree step size
	double start_phi = /*0.*/PI-.2, end_phi = 2*PI, step_phi = PI/360.;
	const int unsigned count = 1;  // number of particles per step
	// should the z vertex be randomized regarding to the target length? if yes, change target length accordingly [cm]
	bool vtx = true;
	const double target_length = 10.;
	// print more information
	bool dbg = false;
	/*
	 * End of simulation specific user modifications
	 */

	// set particle mass [GeV]
	double m;
	if (id_part == proton)
		m = MASS_PROTON/1000.;
	else if (id_part == electron || id_part == positron)
		m = MASS_ELECTRON/1000.;
	else if (id_part == muon || id_part == antimuon)
		m = MASS_MUON/1000.;
	else
		m = 0.;

	int n_part = 1;  // number of particle(s)
	char var_names[128];

	// vertex and beam information
	const double x_vtx = 0., y_vtx = 0., px_bm = 0., py_bm = 0., pz_bm = 1.;
	double z_vtx = 0., pt_bm = .1, en_bm = .1;  // 100 MeV beam, just an arbitrary value [GeV]

	// names for the branches
	sprintf(var_names, "X_vtx:Y_vtx:Z_vtx:Px_bm:Py_bm:Pz_bm:Pt_bm:En_bm:Px_l%02d%02d:Py_l%02d%02d:Pz_l%02d%02d:Pt_l%02d%02d:En_l%02d%02d",
		n_part, id_part, n_part, id_part, n_part, id_part, n_part, id_part, n_part, id_part);
	TFile f(name, "RECREATE");
	if (!f.IsOpen()) {
		fprintf(stderr, "[ERROR] Can't create file %s: %s\n", name, strerror(errno));
		exit(1);
	}
	TNtupleD tpl("h1", "mkin MC file", var_names);
	if (dbg) {
		tpl.Print();
		std::cout << "#args: " << tpl.GetNvar() << std::endl;
	}

	// random number generator
	TRandom3 rand(0);

	double st, sp, ct, cp;
	double px, py, pz, pt, en;
	//Double_t buffer[8 + 5*n_part];  // 8 parameters for vertex (3) and beam (5) + 5 parameters per particle (px, py, pz, pt, e)
	Double_t* buffer;
	buffer = (Double_t*)malloc((8 + 5*n_part)*sizeof(Double_t));
	//memset(buffer, 0, (8 + 5*n_part)*sizeof(Double_t));
	unsigned long long int n_events = 0;
	// write vertex and beam information (fixed values) to array
	/**buffer++ = x_vtx;
	*buffer++ = y_vtx;
	*buffer++ = z_vtx;
	*buffer++ = px_bm;
	*buffer++ = py_bm;
	*buffer++ = pz_bm;
	*buffer++ = pt_bm;
	*buffer++ = en_bm;*/
	buffer[0] = x_vtx;
	buffer[1] = y_vtx;
	buffer[2] = z_vtx;
	buffer[3] = px_bm;
	buffer[4] = py_bm;
	buffer[5] = pz_bm;
	buffer[6] = pt_bm;
	buffer[7] = en_bm;

	for (double e = start_energy; e < end_energy; e += step_energy) {
		for (double t = start_theta; t <= end_theta; t += step_theta) {
			for (double p = start_phi; p <= end_phi; p += step_phi) {
			//TODO: Phi symmetric -> randomize uniformly?
			// Or use amount of particles per step for phi loop to have all possible phi values, even for low count values:
			// step_phi = (end_phi - start_phi)/count; for (as above)
				st = sin(t);
				sp = sin(p);
				ct = cos(t);
				cp = cos(p);
				px = st*cp;
				py = st*sp;
				pz = ct;
				if (id_part == photon)
					en = pt = e;
				else {
					en = e;  // e total energy, kinetic + m
					pt = sqrt(e*e - m*m);
				}

				buffer[8] = px;
				buffer[9] = py;
				buffer[10] = pz;
				buffer[11] = pt;
				buffer[12] = en;

				for (unsigned int i = 0; i < count; i++) {
					// change Z vertex randomly based on a uniform distribution over the target length
					if (vtx)
						buffer[2] = target_length * rand.Uniform(-.5, .5);
					tpl.Fill(buffer);
					n_events++;
					if (n_events % 10000000 == 0 || (n_events < 10000000 && n_events % 1000000 == 0))
						std::cout << "[INFO] Created " << n_events/1000000 << "M events" << std::endl;
				}
			}
		}
	}

	// write the event ntuple to the output file
	printf("\n[INFO] Writing to %s . . .\n", f.GetName());
	printf("[INFO]     =>  %llu events\n", n_events);
	tpl.Write();
	// close file and free memory
	f.Close();
	free(buffer);

	printf("[INFO] Done!\n");

	return 0;
}
