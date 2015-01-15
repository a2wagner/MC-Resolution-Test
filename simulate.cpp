//{
//	// define the reaction as usual
//	PReaction my_reaction(1.5,"g","p","p g","test",1,0,0,0);

//	my_reaction.Print();  // some infos about the reaction
//	// start simulation
//	my_reaction.Loop(3);  //570 für etwas mehr als 200 Ereignisse nach Geant-Simulation

//	// generate files for accessing the root-file
///*	TFile *f = new TFile("test_sample.root", "OPEN");
//	TTree *t = (TTree*)f->Get("data");
//	t->MakeClass("analysis");*/
//}

// compile with:
// g++ -std=gnu++11 -o simulate simulate.C `root-config --cflags --glibs` -lSpectrum

#include <iostream>
#include <stdio.h>
#include <math.h>

#include <TROOT.h>
#include <TFile.h>
#include <TNtuple.h>
//#include <TLorentzVector.h>

int main(int argc, char **argv)
{
	static const double PI = 3.141592653589793;

	int n_part = 1;  // number of particles
	int id_part = 1;  // id of particle, only one
	const int id_photon = 1;
	const char* name = "test_output.root";
	char var_names[128];

	double start_energy = /*.8*/1.4, end_energy = 1.604, step_energy = .0005;  // GeV
	double start_theta = /*0.*/PI-.2, end_theta = PI, step_theta = PI/360.;  //radians, half a degree step size
	double start_phi = /*0.*/PI-.2, end_phi = 2*PI, step_phi = PI/360.;
	const int unsigned count = 1;  // number of particles per step

	const double x_vtx = 0., y_vtx = 0., z_vtx = 0., px_bm = 0., py_bm = 0.;
	const double pz_bm = 1.;
	double pt_bm = .1, en_bm = .1;  // 100 MeV beam, just an arbitrary value [GeV]

	sprintf(var_names, "X_vtx:Y_vtx:Z_vtx:Px_bm:Py_bm:Pz_bm:Pt_bm:En_bm:Px_l%02d%02d:Py_l%02d%02d:Pz_l%02d%02d:Pt_l%02d%02d:En_l%02d%02d",
		id_part, id_photon, id_part, id_photon, id_part, id_photon, id_part, id_photon, id_part, id_photon);
	TFile f(name, "RECREATE");
	TNtuple tpl("h1", "mkin MC file", var_names);
	tpl.Print();
	std::cout << "args: " << tpl.GetNvar() << std::endl;

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

				for (unsigned int i = 0; i < count; i++) {
					//printf("%d\n", n_events);
					//tpl->Print();
					tpl.Fill(buffer);
					//tpl->Print();
					n_events++;
					if (n_events % 1000000 == 0)
						std::cout << "Created " << n_events/1000000 << "M events" << std::endl;
				}
			}
		}
	}

	// write the event ntuple to the output file
	printf("\n\nWriting to %s . . .\n", f.GetName());
	printf("   => %d events\n", n_events);
	tpl.Write();
	// close file and free memory
	f.Close();
	free(buffer);

	printf("Done!\n");

	return 0;
}
