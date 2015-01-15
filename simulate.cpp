// compile with:
// g++ -std=gnu++11 -O3 -o simulate simulate.cpp `root-config --cflags --glibs` -lSpectrum

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <map>
#include <errno.h>
#include <string>

#include <TROOT.h>
#include <TFile.h>
#include <TNtuple.h>
//#include <TLorentzVector.h>

typedef std::map<const int, double> IDMap;
typedef std::pair<const int, const double> IDPair;
typedef std::map<const int, const double>::iterator IDIter;

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
		else if (strstr(argv[1], "positron"))
			id_part = positron;
		else if (strstr(argv[1], "electron"))
			id_part = electron;
		else if (strstr(argv[1], "antimu"))
			id_part = antimuon;
		else if (strstr(argv[1], "muon"))
			id_part = muon;
		else if (strstr(argv[1], "proton"))
			id_part = proton;
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

	bool dbg = true;

	static const double PI = 3.141592653589793;
	static const double MASS_PROTON = 938.272;
	static const double MASS_ELECTRON = .5109989;
	static const double MASS_MUON = 105.65837;

	IDMap mass;
	mass.insert(IDPair(photon, 0.));
	mass.insert(IDPair(positron, MASS_ELECTRON));
	mass.insert(IDPair(electron, MASS_ELECTRON));
	mass.insert(IDPair(muon, MASS_MUON));
	mass.insert(IDPair(antimuon, MASS_MUON));
	mass.insert(IDPair(proton, MASS_PROTON));

	int n_part = 1;  // number of particle(s)
	char var_names[128];

	double start_energy = /*.8*/1.4, end_energy = 1.604, step_energy = .0005;  // GeV
	double start_theta = /*0.*/PI-.2, end_theta = PI, step_theta = PI/360.;  //radians, half a degree step size
	double start_phi = /*0.*/PI-.2, end_phi = 2*PI, step_phi = PI/360.;
	const int unsigned count = 1;  // number of particles per step

	const double x_vtx = 0., y_vtx = 0., px_bm = 0., py_bm = 0., pz_bm = 1.;
	double z_vtx = 0., pt_bm = .1, en_bm = .1;  // 100 MeV beam, just an arbitrary value [GeV]

	sprintf(var_names, "X_vtx:Y_vtx:Z_vtx:Px_bm:Py_bm:Pz_bm:Pt_bm:En_bm:Px_l%02d%02d:Py_l%02d%02d:Pz_l%02d%02d:Pt_l%02d%02d:En_l%02d%02d",
		n_part, id_part, n_part, id_part, n_part, id_part, n_part, id_part, n_part, id_part);
	TFile f(name, "RECREATE");
	if (!f.IsOpen()) {
		fprintf(stderr, "[ERROR] Can't create file %s: %s\n", name, strerror(errno));
		exit(1);
	}
	TNtuple tpl("h1", "mkin MC file", var_names);
	if (dbg) {
		tpl.Print();
		std::cout << "#args: " << tpl.GetNvar() << std::endl;
	}

	//TLorentzVector p;
	double st, sp, ct, cp;
	double px, py, pz, pt;
	//Float_t buffer[8 + 5*n_part];  // 8 parameters for vertex (3) and beam (5) + 5 parameters per particle (px, py, pz, pt, e)
	Float_t* buffer;
	buffer = (Float_t*)malloc((8 + 5*n_part)*sizeof(Float_t));
	//memset(buffer, 0, (8 + 5*n_part)*sizeof(Float_t));
	unsigned int n_events = 0;
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
				st = sin(t);
				sp = sin(p);
				ct = cos(t);
				cp = cos(p);
				px = e*st*cp;
				py = e*st*sp;
				pz = e*ct;
				pt = e;

				buffer[8] = px;
				buffer[9] = py;
				buffer[10] = pz;
				buffer[11] = pt;
				buffer[12] = e;
//TODO: phi unnötig? (symmetrisch); z vertex uniform verteilen -> theta resolution; andere teilchen hinzufügen, masse etc berücksichtigen (proton, elektron, ...)
				for (unsigned int i = 0; i < count; i++) {
					tpl.Fill(buffer);
					n_events++;
					if (n_events % 1000000 == 0)
						std::cout << "[INFO] Created " << n_events/1000000 << "M events" << std::endl;
				}
			}
		}
	}

	// write the event ntuple to the output file
	printf("\n\n[INFO] Writing to %s . . .\n", f.GetName());
	printf("[INFO]    => %d events\n", n_events);
	tpl.Write();
	// close file and free memory
	f.Close();
	free(buffer);

	printf("[INFO] Done!\n");

	return 0;
}
