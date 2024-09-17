// g++ -o a a.cpp
// g++ -fopenmp -o a a.cpp

//#include <bits/stdc++.h>
#include <vector>
#include <iostream>
#include <random>
#include <fstream>
#include <cassert>
//#define assertm(exp, msg) assert(((void)msg, exp))

// Fixed params
#define n_sbox 4 // Number of sboxes
#define n_bits 4 // Number of bits for each sbox

// Inspection params
#define n_keys 32768 // Number of key generated
#define min_rounds 2 // Minimum number of rounds
#define max_rounds 17 // Maximum number of rounds
#define n_mix 1 // Total number of testable mixing layers

// DEBUG
#define x_MAX -17.0 // Target score for the xor 
#define c_MAX -17.0 // Target score for the xor 

using namespace std;

// map<pair<int,int>, int> full_cache; // Already computed sums for full vectors;
vector<vector<int> > block_cache((1<<n_bits), vector<int>((1<<n_bits)));
// map<pair<int,int>, int> block_cache; // Already computed sums for blocks;

// Sbox S25
short unsigned int Sbox[16] = {0,14,11,1,7,12,9,6,13,3,4,15,2,8,10,5};
// Sbox S26
// short unsigned int Sbox[16] = {0,12,15,1,5,10,9,4,11,7,2,13,6,8,14,3};

int attack(int N_ROUNDS, int N_KEYS, int MIX_ID, int ciph_size, vector<int> L, vector<vector<int> > M, vector<int> in_diff);
int cipher (int x, vector<int>& L, vector<int>& key, int N_ROUNDS);
int sum(int x, int y, vector<vector<int> >& M);
int circ(int x, int y, vector<vector<int> >& M);

int main(){

	int ciph_size = n_sbox * n_bits; // Total size of cipher
	
	vector<double> scoreXor(n_mix, 0); // For each mixing layer, the best score for xor diff
	vector<double> scoreCirc(n_mix, 0); // For each mixing layer, the best score for circ diff

	
	// Set of all possible differences activating ONLY ONE S-BOX
	vector<int> in_diff;

	for(int i=1;i<(1<<n_bits);i++){
	 	for(int j=0;j<n_sbox;j++){
	 		int diff = i << (j*n_bits);
	 		in_diff.push_back(diff);
	 	}
	}

	// ALL POSSIBLE DIFFERENCES
	// vector<int> in_diff;
	
	// for(int i=1;i<(1<<ciph_size);i++){
	//	in_diff.push_back(i);
	// }


	// |=========================|
	// |     Reading circ sum    |
	// |=========================|
	vector<vector<int> > M(n_bits, vector<int>(n_bits));
	ifstream matrices_file("data/circle.txt");
	for(int i=0;i<n_bits;i++){
		for(int j=0;j<n_bits;j++){
			matrices_file >> M[i][j];
		}
	}
	matrices_file.close();

	// |========================|
	// |      Reading Mixing    |
	// |========================|
	// Every number is a column of the layer - we can test multiple ones
	vector<vector<int> > L(n_mix, vector<int>(ciph_size));
	ifstream mixlayers_file("data/mixing_layers.txt");

	for(int layer=0;layer<n_mix;layer++){
		for(int i=0;i<ciph_size;i++){
			mixlayers_file >> L[layer][i];
		}
	}
	mixlayers_file.close();

	if(n_mix>1){
		printf("[MAIN] Loaded %d mixing layers\n", n_mix);
		for(int layer=0;layer<n_mix;layer++){
			printf("[MIXING] ML %d:", layer);
			for(int i=0;i<ciph_size;i++){
				printf(" %d", L[layer][i]);
			}
			printf("\n");
		}
	}
	
	// |========================|
	// |      Precomputation    |
	// |========================|
	for(int i=0;i<(1<<n_bits);i++){
		for(int j=0;j<(1<<n_bits);j++){
			block_cache[i][j] = sum(i, j, M);
		}
	}

	// DEBUG
	// for (auto it=block_cache.begin(); it!=block_cache.end(); ++it) {
	// 	printf("[MAIN] 0x%.1x o 0x%.1x = 0x%.1x \n", it->first.first, it->first.second, it->second);
	// }

	// |======================|
	// |        Attack        | 
	// |======================|

	for(int mix_id=0;mix_id<n_mix;mix_id++){
		for(int n_rounds=min_rounds;n_rounds<=max_rounds;n_rounds++){
			attack(n_rounds, n_keys, mix_id, ciph_size, L[mix_id], M, in_diff);
		}
	}

}

// Compare xor and circ on single mixing
int attack(int N_ROUNDS, int N_KEYS, int MIX_ID, int ciph_size, vector<int> L, vector<vector<int> > M, vector<int> in_diff){

	// |========================|
	// |      Key generation    |
	// |========================|
	
	cout << "[COMPARE] Starting key generation" << endl;
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> distr(0, (1<<ciph_size));

	vector<vector<int> > keys(N_KEYS, vector<int>(N_ROUNDS));
	for (int i=0;i<N_KEYS;i++){
		for(int j=0;j<N_ROUNDS;j++){
			keys[i][j] = distr(gen);
		}
	}
	cout << "[COMPARE] Generated " << N_KEYS << " keys for " << N_ROUNDS << " rounds" << endl;

	// |======================|
	// |     XOR analysis     |
	// |======================|
	
	cout << "[COMPARE] Starting XOR" << endl;
	// Every mixing layer
	cout << "[XOR m" << MIX_ID+1 << " r" << N_ROUNDS << "] Mixing layer " << MIX_ID+1 << "/" << n_mix << " - " << N_ROUNDS << " rounds" << endl;

	double x_max = -1000000;

	// vector<int> in_diff = {0x7000}; // DEBUG

	// Check every possible input difference activating a single sbox
	#pragma omp parallel for reduction(max:x_max) // Compute the max for each thread and then the global one
	for(int e = 0;e < in_diff.size();e++){
		int D = in_diff[e];
		// printf("[XOR m%d r%d] In diff: 0x%.4x (%d/%d)\n",MIX_ID+1,N_ROUNDS,D,e+1,(int)(in_diff.size()));
		
		// Possible outcomes
		vector<int> cnt_ddt((1<<ciph_size), 0.0);
		
		for(int k=0;k < n_keys;k++){ 
			for(int x=0;x < (1<<ciph_size);x++){
				int out_D = (cipher(x,L,keys[k],N_ROUNDS))^(cipher(x^D,L,keys[k],N_ROUNDS));
				if(out_D == 0){
					printf("x: 0x%.4x - in_diff: 0x%.4x - N_ROUNDS: %d\n", x, D, N_ROUNDS);
					printf("cip(x): 0x%.4x - cip(x^D): 0x%.4x\n", cipher(x,L,keys[k],N_ROUNDS), cipher(x^D,L,keys[k],N_ROUNDS));
					for(int i=0;i<keys[k].size();i++){
						printf("k[%d]: 0x%.4x\n",i,keys[k][i]);
					}
				}
				assert(out_D != 0);
				cnt_ddt[out_D] += 1;

			}
		}

		// Averaging
		double tmp_best = -1000000;
		for(int i=0;i<(1<<ciph_size);i++){
			double avg = (double)cnt_ddt[i] / n_keys;
			double tmp_score = (double)log2((double)avg/(1<<ciph_size));
			tmp_best = max(tmp_best, tmp_score);
			x_max = max(x_max, tmp_score);
			// if (x_max == tmp_score){ // DEBUG
				// printf("[XOR m%d r%d - NEW MAX] 0x%.4x --> 0x%.4x %f p=2^(%f) \n",MIX_ID+1,N_ROUNDS,D,i,avg,tmp_score);
			// }
			if(tmp_score>x_MAX){
				printf("[XOR m%d r%d - HIGH] 0x%.4x --> 0x%.4x %f p=2^(%f) \n",MIX_ID+1,N_ROUNDS,D,i,avg,tmp_score);
			}

		}
		printf("[XOR m%d r%d] In diff: 0x%.4x (%d/%d) - score: %f \n",MIX_ID+1,N_ROUNDS,D,e+1,(int)(in_diff.size()), tmp_best);
		// cout << "[XOR m" << MIX_ID+1 << " r" << N_ROUNDS << "] xmax[" << e+1 << "]: " << x_max << endl;
	}

	printf("[SCORE] Mixing: %d - Rounds: %d - xor_score: %f\n", MIX_ID+1, N_ROUNDS, x_max);
	
	cout << "[COMPARE] XOR differences computed" << endl;

	
	// |=====================|
	// |    CIRC analysis    |
	// |=====================|
	
	cout << "[COMPARE] Starting CIRC" << endl;
	// Every mixing layer
	cout << "[CIRC m" << MIX_ID+1 << " r" << N_ROUNDS << "] Mixing layer " << MIX_ID+1 << "/" << n_mix << " - " << N_ROUNDS << " rounds" << endl;
	double c_max = -1000000;

	// vector<int> in_diff = {0x7000}; // DEBUG

	// Check every possible input difference activating a single sbox
	#pragma omp parallel for reduction(max:c_max) // Compute the max for each thread and then the global one
	for(int e = 0;e < in_diff.size();e++){
		int D = in_diff[e];
		// printf("[CIRC m%d r%d] In diff: 0x%.4x (%d/%d)\n",MIX_ID+1,N_ROUNDS,D,e+1,(int)(in_diff.size()));
		
		// Possible outcomes
		vector<int> cnt_ddt((1<<ciph_size), 0.0);
		
		for(int k=0;k < n_keys;k++){ 
			for(int x=0;x < (1<<ciph_size);x++){
				int xD = circ(x,D,M);
				int fx = cipher(x,L,keys[k],N_ROUNDS);
				int fxd = cipher(xD,L,keys[k],N_ROUNDS);
				int out_D = circ(fx,fxd,M);
				assert(out_D != 0); // DEBUG
				cnt_ddt[out_D]++;
			}
		}

		// Averaging
		double tmp_best = -1000000;
		for(int i=0;i<(1<<ciph_size);i++){
			double avg = (double)cnt_ddt[i] / n_keys;
			double tmp_score = (double)log2((double)avg/(1<<ciph_size));
			c_max = max(c_max, tmp_score);
			tmp_best = max(tmp_best, tmp_score);
			// if (c_max == tmp_score){ // DEBUG
				// printf("[CIRC m%d r%d - NEW MAX] 0x%.4x --> 0x%.4x %f p=2^(%f) \n",MIX_ID+1,N_ROUNDS,D,i,avg,tmp_score);
			// }
			if(tmp_score>c_MAX){
				printf("[CIRC m%d r%d - HIGH] 0x%.4x --> 0x%.4x %f p=2^(%f) \n",MIX_ID+1,N_ROUNDS,D,i,avg,tmp_score);
			}

		}
		printf("[CIRC m%d r%d] In diff: 0x%.4x (%d/%d) - score: %f \n",MIX_ID+1,N_ROUNDS,D,e+1,(int)(in_diff.size()), tmp_best);
		// cout << "[CIRC m" << MIX_ID+1 << " r" << N_ROUNDS << "] cmax[" << e+1 << "]: " << c_max << endl;
	}

	printf("[SCORE] Mixing: %d - Rounds: %d - circ_score: %f\n", MIX_ID+1, N_ROUNDS, c_max);
	
	cout << "[COMPARE] Circ differences computed" << endl;

	// |========================|
	// |      Final Scores      |
	// |========================|

	printf("[SCORE] Mixing: %d - Rounds: %d - diff: %1.3f - xor: %1.3f - circ: %1.3f\n",MIX_ID+1,N_ROUNDS,c_max-x_max,x_max,c_max); 
	cout << endl;

	return 0;
}


// Encryption
int cipher (int x, vector<int>& L, vector<int>& key, int N_ROUNDS){
	int ciph_size = n_sbox * n_bits; // Total size of cipher
	int tmp, res;
	int mask = (1<<n_bits)-1; // 1111 - mask for a single sbox
	
	for(int round=0;round<N_ROUNDS;round++){
		tmp = 0;
		res = 0;
		// Sbox
		for(int i=0;i<n_sbox;i++){
			res ^= (Sbox[ (x >> (i*n_bits)) & mask ] << (i*n_bits));
		}
		// Mixing
		for(int i=0;i<ciph_size;i++){
			tmp ^= ((__builtin_popcount(res & L[i])%2) << (ciph_size-i-1));
		}
		// Key add
		x = tmp^key[round];
	}
	return x;
}

// Parallel sum
int circ(int x, int y, vector<vector<int> >& M){
	// FULL CACHE
	// pair<int,int> in_p(x,y);
	// if (full_cache.count(in_p)){
	// 	return full_cache[in_p];
	// }
	
	int res=0, tmp_x, tmp_y, tmp_sum;
	int mask = (1<<n_bits)-1;
	
	for(int i=0;i<n_sbox;i++){
		tmp_x = (x >> (i*n_bits))&mask;
		tmp_y = (y >> (i*n_bits))&mask;
		// tmp_sum = sum(tmp_x, tmp_y, M); // Using the sum function
		tmp_sum = block_cache[tmp_x][tmp_y]; // Cached

		res ^= (tmp_sum << (i*n_bits));
	}
	// full_cache[in_p] = res;	
	return res;
}

// Alternative sum
int sum(int x, int y, vector<vector<int> >& M){
	vector<int> A(n_bits, 0x0);
	int res;
	for(int i=0;i<n_bits;i++){
		if ((x>>i) & 0x1) { //i-th component of x
			for (int j=0;j<n_bits;j++){
				A[j] ^= M[n_bits-i-1][j];
			}
		}
	}
	res = 0;
	for(int i=0;i<n_bits;i++){
		res ^= (__builtin_popcount(y&A[i])%2) << (n_bits-i-1);
	}
	if((__builtin_popcount(x)%2) == 1){
		return (res^x);
	}else{
		return (res^x^y);
	}
}
