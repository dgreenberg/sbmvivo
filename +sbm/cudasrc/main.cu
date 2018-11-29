#include <stdlib.h>
#include <stdint.h>
#include <string>
#include "cuda5s.h"

int main(int argc, char **argv) {
	if (argc < 2) return 0;
	printf("\nopening data file...");
	FILE * pFile = fopen(argv[1],"r");
	if (pFile == NULL) {
		printf("failed to open file: %s\n", argv[1]);
		perror("Error");
		return -1;
	}
	printf("done.\n");
	printf("Reading 5s pf data and allocating GPU space...");
	gpu5s_problem * g = init_gpu5sproblem_fromraw(pFile);
	fclose(pFile);

	//FIXME: add an option to initialize everything on the GPU, so we can debug that too

	if (g == NULL) {
		printf("read failed!\n");
		return -2;
	}
	printf("done.\n");
	if (argc > 2) {
		unsigned long long seed = (unsigned long long) atoi(argv[2]);
		printf("initializing with seed %d\n", (int) seed);
		cudaseedrng(g->options.nblocks_pf, g->d.rngstates, g->options.nblocks_pf * BLOCKSIZE_PF, seed);
	}
	pushparameterstodevice(g);

	printf("Running cuda pf code. ");
	int mainloopresult = runpfmainloop(g);
	float marglik = gpu5s_marglik(g);
	printf("Marginal likelihood: %f\n",marglik);

	std::string filename(argv[1]);
	std::string outname;
	size_t lastdot = filename.find_last_of(".");
	if (lastdot == std::string::npos) {
		outname = filename + "_results";
	} else {
		outname = filename.substr(0, lastdot) + "_results." + filename.substr(lastdot + 1);
	}
	pFile = fopen(outname.c_str(),"w");
	if (pFile == NULL) {
		printf("failed to write results!\n");
		return -1;
	}
	int saveresult = save_gpu5sresults_toraw(pFile, g);
	fclose(pFile);
	free(g->h.fobs); g->h.fobs = NULL;
	free(g->h.u); g->h.u = NULL;
	destroy_gpu5s_problem(g);
	printf("\n");
	return saveresult;
}
