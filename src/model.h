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

#ifndef	_MODEL_H
#define	_MODEL_H

#include "constants.h"
#include "dataset.h"
#include "bigdouble.h"

using namespace std;

// LDA model
class model {
public:
    // fixed options
    string wordmapfile;		// file that contains word map [string -> integer id]
    string trainlogfile;	// training log file
    string tassign_suffix;	// suffix for topic assignment file
    string theta_suffix;	// suffix for theta file
    string phi_suffix;		// suffix for phi file
    string others_suffix;	// suffix for file containing other parameters
    string twords_suffix;	// suffix for file containing words-per-topics

    string dir;			// model directory
    string dfile;		// data file    
    string model_name;		// model name
    int model_status;		// model status:
				// MODEL_STATUS_UNKNOWN: unknown status
				// MODEL_STATUS_EST: estimating from scratch
				// MODEL_STATUS_ESTC: continue to estimate the model from a previous one
				// MODEL_STATUS_INF: do inference

    dataset * ptrndata;	// pointer to training dataset object
    dataset * pnewdata; // pointer to new dataset object

    mapid2word id2word; // word map [int => string]
    
    // --- model parameters and variables ---    
    int M; // dataset size (i.e., number of docs)
    int V; // vocabulary size
    int K; // number of topics
    double alpha, beta; // LDA hyperparameters 
    int niters; // number of Gibbs sampling iterations
    int liter; // the iteration at which the model was saved
    int savestep; // saving period
    int twords; // print out top words per each topic
    int withrawstrs;

    double * p; // temp variable for sampling
    int ** z; // topic assignments for words, size M x doc.size()
    int ** nw; // cwt[i][j]: number of instances of word/term i assigned to topic j, size V x K
    int ** nd; // na[i][j]: number of words in document i assigned to topic j, size M x K
    int * nwsum; // nwsum[j]: total number of words assigned to topic j, size K
    int * ndsum; // nasum[i]: total number of words in document i, size M
    double ** theta; // theta: document-topic distributions, size M x K
    double ** phi; // phi: topic-word distributions, size K x V
    
    // for inference only
    int inf_liter;
    int newM;
    int newV;
    int ** newz;
    int ** newnw;
    int ** newnd;
    int * newnwsum;
    int * newndsum;
    double ** newtheta;
    double ** newphi;
    // --------------------------------------
	// for DMM
	int DMM_K, DMM_V, DMM_M;
	int * DMM_zd; // the topic of document i:DMM_zd[i];
	int * DMM_Mz; // number of documents in cluster z;
	//int ** DMM_nwd; // number of occurrences of word w in document d;
	int ** DMM_nwz; // number of occurrences of word w in cluster z;
	int * DMM_nwsum; // number of words in cluster z; 
	int * DMM_ndsum; // number of words in document d;
	int DMM_iter; // number of iterations for DMM;
	int ** DMM_doc_topicsum;
	int * DMM_end_doc_topic;// the final topic of doc;
	double DMM_alpha;
	double DMM_beta;
 	vector<int> *doc_word_id;
	vector<int> *doc_word_count;// about word w in trans_document d
	vector<int> *unique_doc_word_id;
	vector<int> *unique_doc_word_count;
	vector<BigDouble> Index_numerator;
	vector<BigDouble> Index_denominator;
	const int DMM_INDEX_SIZE = 2000000;
	BigDouble* prob;

	dataset *newdoc_forDMM;
	document** temp_point;
	vector<int> *temp_doc;
	vector<int> *word_belong;//mark the word belong to which doc
	int ** clusterWordFreq;
	int * clusterWordTotalFreq;
	// --------------------------------------
    
    model() {
		set_default_values();
    }
          
    ~model();
    
    // set default values for variables
    void set_default_values();   

    // parse command line to get options
    int parse_args(int argc, char ** argv);
    
    // initialize the model
    int init(int argc, char ** argv);
    
    // load LDA model to continue estimating or to do inference
    int load_model(string model_name);
    
    // save LDA model to files
    // model_name.tassign: topic assignments for words in docs
    // model_name.theta: document-topic distributions
    // model_name.phi: topic-word distributions
    // model_name.others: containing other parameters of the model (alpha, beta, M, V, K)
    int save_model(string model_name);
    int save_model_tassign(string filename);
    int save_model_theta(string filename);
    int save_model_phi(string filename);
    int save_model_others(string filename);
    int save_model_twords(string filename);
    
    // saving inference outputs
    int save_inf_model(string model_name);
    int save_inf_model_tassign(string filename);
    int save_inf_model_newtheta(string filename);
    int save_inf_model_newphi(string filename);
    int save_inf_model_others(string filename);
    int save_inf_model_twords(string filename);
    
    // init for estimation
    int init_est();
    int init_estc();
	
    // estimate LDA model using Gibbs sampling
    void estimate();
    int sampling(int m, int n);
    void compute_theta();
    void compute_phi();
    
    // init for inference
    int init_inf();
    // inference for new (unseen) data based on the estimated LDA model
    void inference();
    int inf_sampling(int m, int n);
    void compute_newtheta();
    void compute_newphi();

	//try to use LDA+DMM	
	void trans_into_new_doc();
	int DMM_init(dataset *doc);
	void DMM_estimate();
	int DMM_sampling(int m);
	void build_Index(int len1, int len2);
	void get_uninque_word();
	double calPerplexity(vector<int> *doc_word_id, vector<int> * doc_word_count, int *ClusterWordTotalFreq, int **ClusterWordFreq);
	void DMM_assign_topic();
	int chooseK(vector<double> probs);
	vector<double> normalizeProb(BigDouble* prob);
	void Baseline_GSDMM();
};

#endif

