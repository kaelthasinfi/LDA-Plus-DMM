/*
 * Copyright (C) 2007 by
 * 
 * 	Xuan-Hieu Phan
 *	hieuxuan@ecei.tohoku.ac.jp or pxhieu@gmail.com
 * 	Graduate School of Information Sciences
 * 	Tohoku University
 *
 * GibbsLDA++ is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation; either version 2 of the License,
 * or (at your option) any later version.
 *
 * GibbsLDA++ is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GibbsLDA++; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 */

/* 
 * References:
 * + The Java code of Gregor Heinrich (gregor@arbylon.net)
 *   http://www.arbylon.net/projects/LdaGibbsSampler.java
 * + "Parameter estimation for text analysis" by Gregor Heinrich
 *   http://www.arbylon.net/publications/text-est.pdf
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "constants.h"
#include "strtokenizer.h"
#include "utils.h"
#include "dataset.h"
#include "model.h"
#include <algorithm>
//#include "bigdouble.h"

using namespace std;

model::~model() {
    if (p) {
		delete p;
    }

    if (ptrndata) {
		delete ptrndata;
    }
    
    if (pnewdata) {
		delete pnewdata;
    }

    if (z) {
		for (int m = 0; m < M; m++) {
			if (z[m]) {
				delete z[m];
			}
		}
    }
    
    if (nw) {
		for (int w = 0; w < V; w++) {
			if (nw[w]) {
				delete nw[w];
			}
		}
    }

    if (nd) {
		for (int m = 0; m < M; m++) {
			if (nd[m]) {
				delete nd[m];
			}
		}
    } 
    
    if (nwsum) {
		delete nwsum;
    }   
    
    if (ndsum) {
		delete ndsum;
    }
    
    if (theta) {
		for (int m = 0; m < M; m++) {
			if (theta[m]) {
				delete theta[m];
			}
		}
    }
    
    if (phi) {
		for (int k = 0; k < K; k++) {
			if (phi[k]) {
				delete phi[k];
			}
		}
    }

    // only for inference
    if (newz) {
		for (int m = 0; m < newM; m++) {
			if (newz[m]) {
				delete newz[m];
			}
		}
    }
    
    if (newnw) {
		for (int w = 0; w < newV; w++) {
			if (newnw[w]) {
				delete newnw[w];
			}
		}
    }

    if (newnd) {
		for (int m = 0; m < newM; m++) {
			if (newnd[m]) {
				delete newnd[m];
			}
		}
    } 
    
    if (newnwsum) {
		delete newnwsum;
    }   
    
    if (newndsum) {
		delete newndsum;
    }
    
    if (newtheta) {
		for (int m = 0; m < newM; m++) {
			if (newtheta[m]) {
				delete newtheta[m];
			}
		}
    }
    
    if (newphi) {
		for (int k = 0; k < K; k++) {
			if (newphi[k]) {
				delete newphi[k];
			}
		}
    }
}

void model::set_default_values() {
    wordmapfile = "wordmap.txt";
    trainlogfile = "trainlog.txt";
    tassign_suffix = ".tassign";
    theta_suffix = ".theta";
    phi_suffix = ".phi";
    others_suffix = ".others";
    twords_suffix = ".twords";
    
    dir = "./";
    dfile = "trndocs.dat";
    model_name = "model-final";    
    model_status = MODEL_STATUS_UNKNOWN;
    
    ptrndata = NULL;
    pnewdata = NULL;
    
    M = 0;
    V = 0;
    K = 100;
    alpha = 50.0 / K;
    beta = 0.1;
    niters = 2000;
    liter = 0;
    savestep = 200;    
    twords = 0;
    withrawstrs = 0;
    
    p = NULL;
    z = NULL;
    nw = NULL;
    nd = NULL;
    nwsum = NULL;
    ndsum = NULL;
    theta = NULL;
    phi = NULL;
    
    newM = 0;
    newV = 0;
    newz = NULL;
    newnw = NULL;
    newnd = NULL;
    newnwsum = NULL;
    newndsum = NULL;
    newtheta = NULL;
    newphi = NULL;

	DMM_K = 40;
}

int model::parse_args(int argc, char ** argv) {
    return utils::parse_args(argc, argv, this);
}

int model::init(int argc, char ** argv) {
    // call parse_args
    if (parse_args(argc, argv)) {
		return 1;
    }
    
    if (model_status == MODEL_STATUS_EST) {
		// estimating the model from scratch
		if (init_est()) {
			return 1;
	}
	
    } else if (model_status == MODEL_STATUS_ESTC) {
		
			// estimating the model from a previously estimated one
		if (init_estc()) {
			return 1;
	}
	
    } else if (model_status == MODEL_STATUS_INF) {
		// do inference
		if (init_inf()) {
			return 1;
		}
    }
    
    return 0;
}

int model::load_model(string model_name) {
    int i, j;
    
    string filename = dir + model_name + tassign_suffix;
    FILE * fin = fopen(filename.c_str(), "r");
    if (!fin) {
		printf("Cannot open file %d to load model!\n", filename.c_str());
		return 1;
    }
    
    char buff[BUFF_SIZE_LONG];
    string line;

    // allocate memory for z and ptrndata
    z = new int*[M];
    ptrndata = new dataset(M);
    ptrndata->V = V;

    for (i = 0; i < M; i++) {
		char * pointer = fgets(buff, BUFF_SIZE_LONG, fin);
		if (!pointer) {
	   		printf("Invalid word-topic assignment file, check the number of docs!\n");
	    	return 1;
		}
	
		line = buff;
		strtokenizer strtok(line, " \t\r\n");
		int length = strtok.count_tokens();
	
		vector<int> words;
		vector<int> topics;
		for (j = 0; j < length; j++) {
	    	string token = strtok.token(j);
    
	    	strtokenizer tok(token, ":");
	    	if (tok.count_tokens() != 2) {
				printf("Invalid word-topic assignment line!\n");
				return 1;
	    	}
	    
	    	words.push_back(atoi(tok.token(0).c_str()));
	    	topics.push_back(atoi(tok.token(1).c_str()));
		}
	
		// allocate and add new document to the corpus
		document * pdoc = new document(words);
		ptrndata->add_doc(pdoc, i);
	
		// assign values for z
		z[i] = new int[topics.size()];
		for (j = 0; j < topics.size(); j++) {
	    	z[i][j] = topics[j];
		}
    }   
    
    fclose(fin);
    
    return 0;
}

int model::save_model(string model_name) {
    if (save_model_tassign(dir + model_name + tassign_suffix)) {
		return 1;
    }
    
    if (save_model_others(dir + model_name + others_suffix)) {
		return 1;
    }
    
    if (save_model_theta(dir + model_name + theta_suffix)) {
		return 1;
    }
    
    if (save_model_phi(dir + model_name + phi_suffix)) {
		return 1;
    }
    
    if (twords > 0) {
		if (save_model_twords(dir + model_name + twords_suffix)) {
			return 1;
		}
    }
    
    return 0;
}

int model::save_model_tassign(string filename) {
    int i, j;
    
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
		printf("Cannot open file %s to save!\n", filename.c_str());
		return 1;
    }

    // wirte docs with topic assignments for words
    for (i = 0; i < ptrndata->M; i++) {    
		for (j = 0; j < ptrndata->docs[i]->length; j++) {
			fprintf(fout, "%d:%d ", ptrndata->docs[i]->words[j], z[i][j]);
		}
		fprintf(fout, "\n");
    }

    fclose(fout);
    
    return 0;
}

int model::save_model_theta(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
		printf("Cannot open file %s to save!\n", filename.c_str());
		return 1;
    }
    
    for (int i = 0; i < M; i++) {
		for (int j = 0; j < K; j++) {
			fprintf(fout, "%f ", theta[i][j]);
		}
		fprintf(fout, "\n");
    }
    
    fclose(fout);
    
    return 0;
}

int model::save_model_phi(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
		printf("Cannot open file %s to save!\n", filename.c_str());
		return 1;
    }
    
    for (int i = 0; i < K; i++) {
		for (int j = 0; j < V; j++) {
			fprintf(fout, "%f ", phi[i][j]);
		}
		fprintf(fout, "\n");
    }
    
    fclose(fout);    
    
    return 0;
}

int model::save_model_others(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
		printf("Cannot open file %s to save!\n", filename.c_str());
		return 1;
    }

    fprintf(fout, "alpha=%f\n", alpha);
    fprintf(fout, "beta=%f\n", beta);
    fprintf(fout, "ntopics=%d\n", K);
    fprintf(fout, "ndocs=%d\n", M);
    fprintf(fout, "nwords=%d\n", V);
    fprintf(fout, "liter=%d\n", liter);
    
    fclose(fout);    
    
    return 0;
}

int model::save_model_twords(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
		printf("Cannot open file %s to save!\n", filename.c_str());
		return 1;
    }
    
    if (twords > V) {
		twords = V;
    }
    mapid2word::iterator it;
    
    for (int k = 0; k < K; k++) {
		vector<pair<int, double> > words_probs;
		pair<int, double> word_prob;
		for (int w = 0; w < V; w++) {
			word_prob.first = w;
			word_prob.second = phi[k][w];
			words_probs.push_back(word_prob);
		}
    
        // quick sort to sort word-topic probability
		utils::quicksort(words_probs, 0, words_probs.size() - 1);
	
		fprintf(fout, "Topic %dth:\n", k);
		for (int i = 0; i < twords; i++) {
			it = id2word.find(words_probs[i].first);
			if (it != id2word.end()) {
				fprintf(fout, "\t%s   %f\n", (it->second).c_str(), words_probs[i].second);
			}
		}
    }
    
    fclose(fout);    
    
    return 0;    
}

int model::save_inf_model(string model_name) {
    if (save_inf_model_tassign(dir + model_name + tassign_suffix)) {
		return 1;
    }
    
    if (save_inf_model_others(dir + model_name + others_suffix)) {
		return 1;
    }
    
    if (save_inf_model_newtheta(dir + model_name + theta_suffix)) {
		return 1;
    }
    
    if (save_inf_model_newphi(dir + model_name + phi_suffix)) {
		return 1;
    }

    if (twords > 0) {
		if (save_inf_model_twords(dir + model_name + twords_suffix)) {
			return 1;
		}
    }
    
    return 0;
}

int model::save_inf_model_tassign(string filename) {
    int i, j;
    
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
		printf("Cannot open file %s to save!\n", filename.c_str());
		return 1;
    }

    // wirte docs with topic assignments for words
    for (i = 0; i < pnewdata->M; i++) {    
		for (j = 0; j < pnewdata->docs[i]->length; j++) {
			fprintf(fout, "%d:%d ", pnewdata->docs[i]->words[j], newz[i][j]);
		}
		fprintf(fout, "\n");
    }

    fclose(fout);
    
    return 0;
}

int model::save_inf_model_newtheta(string filename) {
    int i, j;

    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
		printf("Cannot open file %s to save!\n", filename.c_str());
		return 1;
    }
    
    for (i = 0; i < newM; i++) {
		for (j = 0; j < K; j++) {
			fprintf(fout, "%f ", newtheta[i][j]);
		}
		fprintf(fout, "\n");
    }
    
    fclose(fout);
    
    return 0;
}

int model::save_inf_model_newphi(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
		printf("Cannot open file %s to save!\n", filename.c_str());
		return 1;
    }
    
    for (int i = 0; i < K; i++) {
		for (int j = 0; j < newV; j++) {
			fprintf(fout, "%f ", newphi[i][j]);
		}
		fprintf(fout, "\n");
    }
    
    fclose(fout);    
    
    return 0;
}

int model::save_inf_model_others(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
		printf("Cannot open file %s to save!\n", filename.c_str());
		return 1;
    }

    fprintf(fout, "alpha=%f\n", alpha);
    fprintf(fout, "beta=%f\n", beta);
    fprintf(fout, "ntopics=%d\n", K);
    fprintf(fout, "ndocs=%d\n", newM);
    fprintf(fout, "nwords=%d\n", newV);
    fprintf(fout, "liter=%d\n", inf_liter);
    
    fclose(fout);    
    
    return 0;
}

int model::save_inf_model_twords(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
		printf("Cannot open file %s to save!\n", filename.c_str());
		return 1;
    }
    
    if (twords > newV) {
		twords = newV;
    }
    mapid2word::iterator it;
    map<int, int>::iterator _it;
    
    for (int k = 0; k < K; k++) {
		vector<pair<int, double> > words_probs;
		pair<int, double> word_prob;
		for (int w = 0; w < newV; w++) {
			word_prob.first = w;
			word_prob.second = newphi[k][w];
			words_probs.push_back(word_prob);
		}
		
			// quick sort to sort word-topic probability
		utils::quicksort(words_probs, 0, words_probs.size() - 1);
		
		fprintf(fout, "Topic %dth:\n", k);
		for (int i = 0; i < twords; i++) {
				_it = pnewdata->_id2id.find(words_probs[i].first);
				if (_it == pnewdata->_id2id.end()) {
				continue;
			}
			it = id2word.find(_it->second);
			if (it != id2word.end()) {
				fprintf(fout, "\t%s   %f\n", (it->second).c_str(), words_probs[i].second);
			}
		}
    }
    
    fclose(fout);    
    
    return 0;    
}


int model::init_est() {
    int m, n, w, k;

    p = new double[K];

    // + read training data
    ptrndata = new dataset;
    if (ptrndata->read_trndata(dir + dfile, dir + wordmapfile)) {
        printf("Fail to read training data!\n");
        return 1;
    }
		
    // + allocate memory and assign values for variables
    M = ptrndata->M;
    V = ptrndata->V;
    // K: from command line or default value
    // alpha, beta: from command line or default values
    // niters, savestep: from command line or default values

    nw = new int*[V];
    for (w = 0; w < V; w++) {
        nw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    nw[w][k] = 0;
        }
    }
	
    nd = new int*[M];
    for (m = 0; m < M; m++) {
        nd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    nd[m][k] = 0;
        }
    }
	
    nwsum = new int[K];
		for (k = 0; k < K; k++) {
		nwsum[k] = 0;
    }
    
    ndsum = new int[M];
    for (m = 0; m < M; m++) {
		ndsum[m] = 0;
    }

    srandom(time(0)); // initialize for random number generation
    z = new int*[M];
    for (m = 0; m < ptrndata->M; m++) {
		int N = ptrndata->docs[m]->length;
		z[m] = new int[N];
	
        // initialize for z
        for (n = 0; n < N; n++) {
    	    int topic = (int)(((double)random() / RAND_MAX) * K);
    	    z[m][n] = topic;
    	    
    	    // number of instances of word i assigned to topic j
    	    nw[ptrndata->docs[m]->words[n]][topic] += 1;
    	    // number of words in document i assigned to topic j
    	    nd[m][topic] += 1;
    	    // total number of words assigned to topic j
    	    nwsum[topic] += 1;
        } 
        // total number of words in document i
        ndsum[m] = N;      
    }
    
    theta = new double*[M];
    for (m = 0; m < M; m++) {
        theta[m] = new double[K];
    }
	
    phi = new double*[K];
    for (k = 0; k < K; k++) {
        phi[k] = new double[V];
    }    
    
    return 0;
}

int model::init_estc() {
    // estimating the model from a previously estimated one
    int m, n, w, k;

    p = new double[K];

    // load moel, i.e., read z and ptrndata
    if (load_model(model_name)) {
		printf("Fail to load word-topic assignmetn file of the model!\n");
		return 1;
    }

    nw = new int*[V];
    for (w = 0; w < V; w++) {
        nw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    nw[w][k] = 0;
        }
    }
	
    nd = new int*[M];
    for (m = 0; m < M; m++) {
        nd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    nd[m][k] = 0;
        }
    }
	
    nwsum = new int[K];
    for (k = 0; k < K; k++) {
		nwsum[k] = 0;
    }
    
    ndsum = new int[M];
	for (m = 0; m < M; m++) {
			ndsum[m] = 0;
    }

    for (m = 0; m < ptrndata->M; m++) {
		int N = ptrndata->docs[m]->length;

		// assign values for nw, nd, nwsum, and ndsum	
        for (n = 0; n < N; n++) {
    	    int w = ptrndata->docs[m]->words[n];
    	    int topic = z[m][n];
    	    
    	    // number of instances of word i assigned to topic j
    	    nw[w][topic] += 1;
    	    // number of words in document i assigned to topic j
    	    nd[m][topic] += 1;
    	    // total number of words assigned to topic j
    	    nwsum[topic] += 1;
        } 
        // total number of words in document i
        ndsum[m] = N;      
    }
	
    theta = new double*[M];
    for (m = 0; m < M; m++) {
        theta[m] = new double[K];
    }
	
    phi = new double*[K];
    for (k = 0; k < K; k++) {
        phi[k] = new double[V];
    }    

    return 0;        
}

void model::estimate() {
    if (twords > 0) {
		// print out top words per topic
		dataset::read_wordmap(dir + wordmapfile, &id2word);
    }

    printf("Sampling %d iterations!\n", niters);

    int last_iter = liter;
	long long count_lda_time = 0;
    for (liter = last_iter + 1; liter <= niters + last_iter; liter++) {
		printf("Iteration %d ...\n", liter);
		
		// for all z_i
		auto start_time = utils::now();	
		for (int m = 0; m < M; m++) {
			for (int n = 0; n < ptrndata->docs[m]->length; n++) {
				// (z_i = z[m][n])
				// sample from p(z_i|z_-i, w)
				int topic = sampling(m, n);
				z[m][n] = topic;
			}
		}
		auto cost_time = utils::now()-start_time;
		count_lda_time += cost_time;
	
		if (savestep > 0) {
			if (liter % savestep == 0) {
				// saving the model
				printf("Saving the model at iteration %d ...\n", liter);
				compute_theta();
				compute_phi();
				save_model(utils::generate_model_name(liter));
			}
		}
    }
    
    printf("Gibbs sampling completed!\n");
	printf("Total time of Gibbs sampling:%llds\n", count_lda_time/1000);
    printf("Saving the final model!\n");
    compute_theta();
    compute_phi();
    liter--;
    save_model(utils::generate_model_name(-1));
}

int model::sampling(int m, int n) {
    // remove z_i from the count variables
    int topic = z[m][n];
    int w = ptrndata->docs[m]->words[n];
    nw[w][topic] -= 1;
    nd[m][topic] -= 1;
    nwsum[topic] -= 1;
    ndsum[m] -= 1;

    double Vbeta = V * beta;
    double Kalpha = K * alpha;    
    // do multinomial sampling via cumulative method
    for (int k = 0; k < K; k++) {
		p[k] = (nw[w][k] + beta) / (nwsum[k] + Vbeta) *
				(nd[m][k] + alpha) / (ndsum[m] + Kalpha);
    }
    // cumulate multinomial parameters
    for (int k = 1; k < K; k++) {
		p[k] += p[k - 1];
    }
    // scaled sample because of unnormalized p[]
    double u = ((double)random() / RAND_MAX) * p[K - 1];
    
    for (topic = 0; topic < K; topic++) {
		if (p[topic] > u) {
			break;
		}
    }
    
    // add newly estimated z_i to count variables
    nw[w][topic] += 1;
    nd[m][topic] += 1;
    nwsum[topic] += 1;
    ndsum[m] += 1;    
    
    return topic;
}

void model::compute_theta() {
    for (int m = 0; m < M; m++) {
		for (int k = 0; k < K; k++) {
			theta[m][k] = (nd[m][k] + alpha) / (ndsum[m] + K * alpha);
		}
    }
}

void model::compute_phi() {
    for (int k = 0; k < K; k++) {
		for (int w = 0; w < V; w++) {
			phi[k][w] = (nw[w][k] + beta) / (nwsum[k] + V * beta);
		}
    }
}

int model::init_inf() {
    // estimating the model from a previously estimated one
    int m, n, w, k;

    p = new double[K];

    // load moel, i.e., read z and ptrndata
    if (load_model(model_name)) {
		printf("Fail to load word-topic assignmetn file of the model!\n");
		return 1;
    }

    nw = new int*[V];
    for (w = 0; w < V; w++) {
        nw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    nw[w][k] = 0;
        }
    }
	
    nd = new int*[M];
    for (m = 0; m < M; m++) {
        nd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    nd[m][k] = 0;
        }
    }
	
    nwsum = new int[K];
    for (k = 0; k < K; k++) {
		nwsum[k] = 0;
    }
    
    ndsum = new int[M];
    for (m = 0; m < M; m++) {
		ndsum[m] = 0;
    }

    for (m = 0; m < ptrndata->M; m++) {
		int N = ptrndata->docs[m]->length;

		// assign values for nw, nd, nwsum, and ndsum	
        for (n = 0; n < N; n++) {
    	    int w = ptrndata->docs[m]->words[n];
    	    int topic = z[m][n];
    	    
    	    // number of instances of word i assigned to topic j
    	    nw[w][topic] += 1;
    	    // number of words in document i assigned to topic j
    	    nd[m][topic] += 1;
    	    // total number of words assigned to topic j
    	    nwsum[topic] += 1;
        } 
        // total number of words in document i
        ndsum[m] = N;      
    }
    
    // read new data for inference
    pnewdata = new dataset;
    if (withrawstrs) {
		if (pnewdata->read_newdata_withrawstrs(dir + dfile, dir + wordmapfile)) {
				printf("Fail to read new data!\n");
				return 1;
		}    
    } else {
		if (pnewdata->read_newdata(dir + dfile, dir + wordmapfile)) {
				printf("Fail to read new data!\n");
				return 1;
		}    
    }
    
    newM = pnewdata->M;
    newV = pnewdata->V;
    
    newnw = new int*[newV];
    for (w = 0; w < newV; w++) {
        newnw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    newnw[w][k] = 0;
        }
    }
	
    newnd = new int*[newM];
    for (m = 0; m < newM; m++) {
        newnd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    newnd[m][k] = 0;
        }
    }
	
    newnwsum = new int[K];
    for (k = 0; k < K; k++) {
		newnwsum[k] = 0;
    }
    
    newndsum = new int[newM];
    for (m = 0; m < newM; m++) {
		newndsum[m] = 0;
    }

    srandom(time(0)); // initialize for random number generation
    newz = new int*[newM];
    for (m = 0; m < pnewdata->M; m++) {
		int N = pnewdata->docs[m]->length;
		newz[m] = new int[N];

		// assign values for nw, nd, nwsum, and ndsum	
        for (n = 0; n < N; n++) {
    	    int w = pnewdata->docs[m]->words[n];
    	    int _w = pnewdata->_docs[m]->words[n];
    	    int topic = (int)(((double)random() / RAND_MAX) * K);
    	    newz[m][n] = topic;
    	    
    	    // number of instances of word i assigned to topic j
    	    newnw[_w][topic] += 1;
    	    // number of words in document i assigned to topic j
    	    newnd[m][topic] += 1;
    	    // total number of words assigned to topic j
    	    newnwsum[topic] += 1;
        } 
        // total number of words in document i
        newndsum[m] = N;      
    }    
    
    newtheta = new double*[newM];
    for (m = 0; m < newM; m++) {
        newtheta[m] = new double[K];
    }
	
    newphi = new double*[K];
    for (k = 0; k < K; k++) {
        newphi[k] = new double[newV];
    }    
    
    return 0;        
}

void model::inference() {
    if (twords > 0) {
		// print out top words per topic
		dataset::read_wordmap(dir + wordmapfile, &id2word);
    }

    printf("Sampling %d iterations for inference!\n", niters);
    
    for (inf_liter = 1; inf_liter <= niters; inf_liter++) {
		printf("Iteration %d ...\n", inf_liter);
		
		// for all newz_i
		for (int m = 0; m < newM; m++) {
			for (int n = 0; n < pnewdata->docs[m]->length; n++) {
				// (newz_i = newz[m][n])
				// sample from p(z_i|z_-i, w)
				int topic = inf_sampling(m, n);
				newz[m][n] = topic;
			}
		}
    }
    
    printf("Gibbs sampling for inference completed!\n");
    printf("Saving the inference outputs!\n");
    compute_newtheta();
    compute_newphi();
    inf_liter--;
    save_inf_model(dfile);
}

int model::inf_sampling(int m, int n) {
    // remove z_i from the count variables
    int topic = newz[m][n];
    int w = pnewdata->docs[m]->words[n];
    int _w = pnewdata->_docs[m]->words[n];
    newnw[_w][topic] -= 1;
    newnd[m][topic] -= 1;
    newnwsum[topic] -= 1;
    newndsum[m] -= 1;
    
    double Vbeta = V * beta;
    double Kalpha = K * alpha;
    // do multinomial sampling via cumulative method
    for (int k = 0; k < K; k++) {
		p[k] = (nw[w][k] + newnw[_w][k] + beta) / (nwsum[k] + newnwsum[k] + Vbeta) *
		    (newnd[m][k] + alpha) / (newndsum[m] + Kalpha);
    }
    // cumulate multinomial parameters
    for (int k = 1; k < K; k++) {
		p[k] += p[k - 1];
    }
    // scaled sample because of unnormalized p[]
    double u = ((double)random() / RAND_MAX) * p[K - 1];
    
    for (topic = 0; topic < K; topic++) {
		if (p[topic] > u) {
			break;
		}
    }
    
    // add newly estimated z_i to count variables
    newnw[_w][topic] += 1;
    newnd[m][topic] += 1;
    newnwsum[topic] += 1;
    newndsum[m] += 1;    
    
    return topic;
}

void model::compute_newtheta() {
    for (int m = 0; m < newM; m++) {
		for (int k = 0; k < K; k++) {
			newtheta[m][k] = (newnd[m][k] + alpha) / (newndsum[m] + K * alpha);
		}
    }
}

void model::compute_newphi() {
    map<int, int>::iterator it;
    for (int k = 0; k < K; k++) {
		for (int w = 0; w < newV; w++) {
			it = pnewdata->_id2id.find(w);
			if (it != pnewdata->_id2id.end()) {
				newphi[k][w] = (nw[it->second][k] + newnw[w][k] + beta) / (nwsum[k] + newnwsum[k] + V * beta);
			}
		}
    }
}

void compute_Index(int len, vector<BigDouble> &num, double basic)
{
	for (int i = 1; i < len; ++i)
	{
		num.push_back(num[i-1] * (basic+i-1));
	}
}

void model::build_Index(int len1, int len2)
{
	auto now_time = utils::now();
	Index_numerator.reserve(static_cast<size_t>(len1));
	Index_denominator.reserve(static_cast<size_t>(len2));
	double basic = V*DMM_beta;
	Index_denominator.emplace_back(1.0);
	compute_Index(len2, Index_denominator, basic);

	basic = DMM_beta;
	Index_numerator.emplace_back(1.0);
	compute_Index(len1, Index_numerator, basic);
	
	auto build_Index_time = (utils::now()-now_time);
	printf("time of building Index:%lldms\n",build_Index_time);
}

//LDA+DMM
int model::DMM_init(dataset *data)
{
	int k, w, m, n;
	DMM_alpha = 0.1;
	DMM_beta = 0.1;
	DMM_iter = 5000;
//	DMM_V = V;
	DMM_M = data->M;
	printf("enter DMM_INIT()\n");
	
	prob = new BigDouble[DMM_K];

	DMM_end_doc_topic = new int[M];
	for (m = 0; m < M; ++m)
		DMM_end_doc_topic[m] = 0;
	DMM_doc_topicsum = new int*[M];
	for (m = 0; m < M; ++m)
	{
		DMM_doc_topicsum[m] = new int[DMM_K];
		for (k = 0; k < DMM_K; ++k)
			DMM_doc_topicsum[m][k] = 0;
	}

	DMM_zd = new int[DMM_M];
	for (m = 0; m < DMM_M; ++m)
		DMM_zd[m] = 0;

	DMM_Mz = new int[DMM_K];
	for (k = 0; k < DMM_K; ++k)
		DMM_Mz[k] = 0;
	
/*	DMM_nwd = new int*[V];
	for (w = 0; w < V; w++)
	{
		DMM_nwd[w] = new int[M];
		for (m = 0; m < M; ++m)
			DMM_nwd[w][m] = 0;
	}*/

	DMM_nwz = new int*[V];
	for (w = 0; w < V; w++)
	{
		DMM_nwz[w] = new int[DMM_K];
		for (k = 0; k < DMM_K; ++k){
			DMM_nwz[w][k] = 0;
		}
	}

	DMM_nwsum = new int[DMM_K];
	for (k = 0; k < DMM_K; ++k)
	{
		DMM_nwsum[k] = 0;
	}

	DMM_ndsum = new int[DMM_M];
	for (m = 0; m < DMM_M; ++m)
	{
		DMM_ndsum[m] = data->docs[m]->length;
		//printf("%d\n", DMM_ndsum[m]);
	}
	doc_word_id = new vector<int>[DMM_M];
	doc_word_count = new vector<int>[DMM_M];

//WHY M, maybe for trans doc;
	if (M!=DMM_M)	
	{
		unique_doc_word_id = new vector<int>[M];
		unique_doc_word_count = new vector<int>[M];
		get_uninque_word();
	}
//	printf("NEW VECTOR\n");

	//printf("%d\n", K);
	for (m = 0; m < DMM_M; ++m)
	{
		int topic = (int)(((double)random() / RAND_MAX) * DMM_K);
		DMM_zd[m] = topic;
		for (n = 0; n < data->docs[m]->length; ++n)
		{
			int word_id = data->docs[m]->words[n];
			//printf("word_topic:%d\n", word_topic);
			//printf("sum_k[word_topic]:%d\n", sum_k[word_topic]);
			//printf("word_id:%d\n", word_id);

			int len_now = doc_word_id[m].size(), flag = 0;	
			for (int i = 0; i < len_now; ++i)
			{
				if (doc_word_id[m][i] == word_id)
				{
					doc_word_count[m][i]++;
					flag = 1;
					break;
				}
			}
			if (flag == 0) 
			{
				doc_word_id[m].push_back(word_id);	
				doc_word_count[m].push_back(1);
			}
		}
		DMM_Mz[topic]++;
		DMM_nwsum[topic] += data->docs[m]->length;
		int d_size = doc_word_id[m].size();
		printf("doc_%d:topic:%d d_size:%d\n", m, topic, d_size);
		for (int i = 0; i < d_size; ++i)
		{
			int word_id = doc_word_id[m][i];
			int word_count = doc_word_count[m][i];
			DMM_nwz[word_id][topic] += word_count;
		}
	}
	return 0;
}
void model::DMM_estimate()
{
	int it;
	printf("DMM_sampling %d iteration\n", DMM_iter);
	build_Index(DMM_INDEX_SIZE, DMM_INDEX_SIZE);
	//build index
	auto dmm_start_time = utils::now();
	for (it = 0; it < DMM_iter; ++it)
	{
		printf("DMM Iteration %d ...\n", it);
		for (int m = 0; m < DMM_M; ++m)
		{
			//auto now_time = utils::now();
			int topic = DMM_sampling(m);
			//auto single_iter_time = utils::now()-now_time;
			//printf("single iteration time of DMM:%lldms\n", single_iter_time);
			printf("doc %d sampling topic:%d\n", m, topic);
			DMM_zd[m] = topic;	
		}
	}
	printf("DMM_sampling Completed\n");
	auto dmm_total_time = (utils::now()-dmm_start_time)/1000;
	printf("Total DMM sampling time:%llds\n",dmm_total_time);
	DMM_assign_topic();
	
	printf("!!!\n");
	int topic_doc_sum[DMM_K];
	for (int i = 0; i < DMM_K; ++i)
		topic_doc_sum[i] = 0;
	for (int i = 0; i < M; ++i)
		topic_doc_sum[DMM_end_doc_topic[i]]++;
	sort(topic_doc_sum, topic_doc_sum+DMM_K);
	for (int i = 0; i < DMM_K; ++i)
		printf("%d ", topic_doc_sum[i]);
	cout << endl;


	printf("LDA+DMM-Perplexity:%lf\n", calPerplexity(unique_doc_word_id, unique_doc_word_count, clusterWordTotalFreq, clusterWordFreq));
}
int model::DMM_sampling(int m)
{
//	printf("DMM_sampling document %d\n", m);
	int topic = DMM_zd[m], doc_word_len = doc_word_id[m].size();
//	printf("document%d-doc_word_len:%d\n", m, doc_word_len);
	//double Index_denominator[doc_word_len], Index_numerator[doc_word_len];
	DMM_Mz[topic]--;
	DMM_nwsum[topic] -= DMM_ndsum[m];
	for (int i = 0; i < doc_word_len; i++)
	{
		int word_id = doc_word_id[m][i];
		//printf("word_id:%d\n", word_id);
		DMM_nwz[word_id][topic] -= doc_word_count[m][i];
	}
	//double Vbeta = V*DMM_beta;
	double Kalpha = DMM_K*DMM_alpha;
	/*Index_numerator[0] = beta;
	Index_denominator[0] = Vbeta;

	for (int i = 1; i < doc_word_len/10; ++i)
	{
		Index_numerator[i] = Index_numerator[i-1]*(beta+i);
		Index_denominator[i] = Index_denominator[i-1]*(Vbeta+i);
		printf("%.12lf %.12lf\n", Index_numerator[i], Index_denominator[i]);
	}*/
	for (int k = 0; k < DMM_K; ++k)
	{
		prob[k] = (DMM_Mz[topic] + alpha)/(M-1+Kalpha);
		BigDouble denominator = 1.0;
		denominator *= Index_denominator[DMM_nwsum[k]+DMM_ndsum[m]]/Index_denominator[DMM_nwsum[k]];
		//cout << "denominator:" << denominator.num << " " << denominator.exp << endl;

		double numv = 1;
		int nume = 0;

		for (int i = 0; i < doc_word_len; ++i)
		{
			int word_id = doc_word_id[m][i];
			int word_count = doc_word_count[m][i];
			
			int index_down = DMM_nwz[word_id][k];
			int index_up = DMM_nwz[word_id][k]+word_count;
			//cout << "Index_down:" << index_down << endl;
			//cout << "Index_up:" << index_up << endl;
			numv *= Index_numerator[index_up].num / Index_numerator[index_down].num;
			//cout << "Index_up_num:" << Index_numerator[index_up].num << endl << "Index_down_num:" << Index_numerator[index_down].num << endl;
			nume += Index_numerator[index_up].exp - Index_numerator[index_down].exp;
			while (numv > BigDouble::LARGE_DOUBLE)
			{
				numv /= BigDouble::LARGE_DOUBLE;
				nume++;
			}
			while (numv < BigDouble::SMALL_DOUBLE)
			{
				numv *= BigDouble::LARGE_DOUBLE;
				nume--;
			}
			//numerator *= (Index_numerator[DMM_nwz[word_id][k]+word_count-1]/Index_numerator[DMM_nwz[word_id][k]-1]);
		}
		BigDouble numerator{numv, nume};	
		//printf("numv:%lf nume:%d\n", numv, nume);
		prob[k] *= (numerator/denominator);
	}
	vector<double> probs = normalizeProb(prob);
	/*cout << "probs:";
	for (int i = 0; i < probs.size(); ++i)
		cout << probs[i] << " ";
	cout << endl;*/
	//printf("!!!!\n");	
	/*for (int k = 1; k < DMM_K; ++k)
	{
		//printf("%d\n", k);
		prob[k] += prob[k-1];
	}*/
	topic = chooseK(probs);	
	//BigDouble u = ((BigDouble)random()/ RAND_MAX * prob[DMM_K-1]);

	//for (topic = 0; topic < DMM_K; topic++)	
		//if (prob[topic] > u) break;	
	
	DMM_Mz[topic]++;
	DMM_nwsum[topic] += DMM_ndsum[m];
	for (int i = 0; i < doc_word_len; i++)
	{
		int word_id = doc_word_id[m][i];
		DMM_nwz[word_id][topic] += doc_word_count[m][i];
	}
	return topic;
}

void model::trans_into_new_doc()
{
	int m, n;
	newdoc_forDMM = new dataset(K);
	temp_point = new document*[K];
	newdoc_forDMM->V = V;
	
	temp_doc = new vector<int>[K];
	word_belong = new vector<int>[K];
	for (m = 0; m < M; ++m)
	{
		for (n = 0; n < ptrndata->docs[m]->length; ++n)
		{
			int word_topic = z[m][n];
			//printf("word_topic:%d\n", word_topic);
			int word_id = ptrndata->docs[m]->words[n];
			temp_doc[word_topic].push_back(word_id);
			word_belong[word_topic].push_back(m);
		}
	}
	//for (int k = 0; k < K; ++k)
//		printf("temp_doc_length:%lld\n", temp_doc[k].size());
	for (int k = 0; k < K; ++k)
	{
		temp_point[k] = new document(temp_doc[k]);
		newdoc_forDMM->add_doc(temp_point[k],k);
//		printf("add_doc_point:%lld\n", temp_point);
	}	
	printf("new_doc_large:%d\n", newdoc_forDMM->M);
	for (int i = 0; i < K; ++i)
	{
		printf("new_doc_%d_Length:%d\n", i+1, newdoc_forDMM->docs[i]->length);
	}
}

//after DMM, assign topic to the original doc.
void model::DMM_assign_topic()
{
	//for (int d = 0; d < DMM_M; ++d)
	//	printf("doc %d cluster %d\n", d, DMM_zd[d]);
	for (int d = 0; d < DMM_M; ++d)	
	{
		int len = word_belong[d].size();
		int topic = DMM_zd[d];
		for (int w = 0; w < len; ++w)
		{
			int doc = word_belong[d][w];
			DMM_doc_topicsum[doc][topic]++;
		}	
	}
	for (int d = 0; d < M; ++d)
	{
		int max = 0, max_z = 0;
		for (int t = 0; t < DMM_K; ++t)
			if (DMM_doc_topicsum[d][t] > max)
			{
				max = DMM_doc_topicsum[d][t];
				max_z = t;
			}
		DMM_end_doc_topic[d] = max_z;
	}
	//recalculate the frequence of word in cluster;
	clusterWordFreq = new int*[V];
	for (int w = 0; w < V; ++w)
	{
		clusterWordFreq[w] = new int[DMM_K];
		for (int k = 0; k < DMM_K; ++k)
			clusterWordFreq[w][k] = 0;
	}
	
	clusterWordTotalFreq = new int[DMM_K];
	for (int k = 0; k < DMM_K; ++k)
		clusterWordTotalFreq[k] = 0;
	//printf("!!!\n");
	for (int d = 0; d < DMM_M; ++d)
	{
		int doc_len = temp_doc[d].size();
		for (int i = 0; i < doc_len; ++i)
		{
			int w = temp_doc[d][i];
			int doc = word_belong[d][i];
			int topic = DMM_end_doc_topic[doc];
			clusterWordFreq[w][topic]++;
			clusterWordTotalFreq[topic]++;
		}
	}	
}
void model::get_uninque_word()
{
	for (int m = 0; m < M; ++m)
	{
		for (int n = 0; n < ptrndata->docs[m]->length; ++n)
		{
			int word_id = ptrndata->docs[m]->words[n];
			int len_now = unique_doc_word_count[m].size();
			int flag = 0;
			for (int i = 0; i < len_now; ++i)
				if (unique_doc_word_id[m][i] == word_id)
				{
					unique_doc_word_count[m][i]++;
					flag = 1;
				}
			if (!flag)
			{
				unique_doc_word_id[m].push_back(word_id);
				unique_doc_word_count[m].push_back(1);
			}
		}
	}
}
//perplexity

double model::calPerplexity(vector<int> *the_doc_word_id, vector<int> *the_doc_word_count, int *ClusterWordTotalFreq, int **ClusterWordFreq)
{
	auto acc_count = 0;
	auto acc = 0.0;

	//for (int d = 0; d < M; ++d)
	//	printf("doc %d belongs to cluster %d\n", d, DMM_end_doc_topic[d]);
	for (int d = 0; d < M; ++d)
	{
		int cluster;
		if (DMM_M != M) 
			cluster = DMM_end_doc_topic[d];
		else 
			cluster = DMM_zd[d];
		int doc_len = the_doc_word_id[d].size();
		for (int i = 0; i < doc_len; ++i)
		{
			int w = the_doc_word_id[d][i];
			int w_count = the_doc_word_count[d][i];
			auto w_gen_prob = ((double)ClusterWordFreq[w][cluster]+DMM_beta)/((double)ClusterWordTotalFreq[cluster]+DMM_beta*V);
			acc += (std::log(1 / w_gen_prob) * w_count);
			acc_count += w_count;
		}
	}
	return std::exp(acc / (double) acc_count);
}

int model::chooseK(vector<double> probs)
{
	double x = utils::random_double();	
	for (int i = 0; i < DMM_K; ++i)
	{
		x -= probs[i];
		if (x < 0)
			return i;
	}
	return DMM_K-1; 
}

vector<double> model::normalizeProb(BigDouble *probs)
{
	BigDouble sum = 0;
	for (int d = 0; d < DMM_K; ++d) sum += probs[d];
	vector<double> result(DMM_K);
	for (int i = 0; i < DMM_K; ++i)
		result[i] = (double) (probs[i] / sum);
/*	printf("probs:");
	for (int i = 0; i < DMM_K; ++i)
	{
		probs[i].print();cout << " ";
	}
	cout << endl;
	double all = 0;
	printf("result:");
	for (int i = 0; i < DMM_K; ++i)
		printf("%lf ", result[i]);
	cout << endl;
	for (int i = 0; i < DMM_K; ++i)
		all += result[i];
	printf("all:%lf\n", all);*/
	return result;
}

void model::Baseline_GSDMM()
{
	DMM_init(ptrndata);	
	printf("Baseline:GSDMM Sampling %d iteration\n", DMM_iter);
	build_Index(DMM_INDEX_SIZE, DMM_INDEX_SIZE);
	//build index
	auto dmm_start_time = utils::now();
	printf("DMM_M:%d\n", DMM_M);
	for (int it = 0; it < DMM_iter; ++it)
	{
		printf("Baseline_GSDMM Iteration %d ...\n", it);
		for (int m = 0; m < DMM_M; ++m)
		{
			int topic = DMM_sampling(m);
			//printf("doc %d sampling topic:%d\n", m, topic);
			DMM_zd[m] = topic;
		}
	}
	printf("Baseline_GSDMM_sampling Completed\n");
	auto Baseline_GSDMM_total_time = (utils::now()-dmm_start_time)/1000;
	printf("Total GSDMM sampling time:%llds\n", Baseline_GSDMM_total_time);
	vector<int> topic_res[DMM_K];
	int topic_doc_sum[DMM_K];
	for (int i = 0; i < DMM_K; ++i)
		topic_doc_sum[i] = 0;
	for (int i = 0; i < DMM_M; ++i)
	{
		topic_res[DMM_zd[i]].push_back(i);
		topic_doc_sum[DMM_zd[i]]++;
	}
	for (int i = 0; i < DMM_K; ++i)
	{
		printf("topic %d:", i);
		for (int j = 0; j < topic_res[i].size(); ++j)
			printf("%d ", topic_res[i][j]);
		cout << endl;
	}
	sort(topic_doc_sum, topic_doc_sum+DMM_K);
	for (int i = 0; i < DMM_K; ++i)
		printf("%d ",topic_doc_sum[i]);
	cout << endl;



	printf("Baseline_GSDMM-Perplexity:%lf\n", calPerplexity(doc_word_id, doc_word_count, DMM_nwsum, DMM_nwz));
}
