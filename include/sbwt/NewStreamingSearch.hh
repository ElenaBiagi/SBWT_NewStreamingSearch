#pragma once

#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include "SBWT.hh"
#include "variants.hh"
#include "sdsl/vectors.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/rmq_support.hpp"

#include "suffix_group_optimization.hh"

using namespace std;
using namespace sbwt;

int64_t get_char_idx(char c){
    switch(c){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

// Returns pairs (a_1, b_1), (a_2, b_2)..., such that 
// - a_i is the length of the longest match ending at query[i]
// - b_i is the colex rank of one arbitrary k-mer having the longest match.
vector<pair<int64_t,int64_t> > new_streaming_search(const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, const sdsl::rmq_succinct_sct<>& rmqLCS, const char* input, int64_t len){ //
    const sdsl::bit_vector& A_bits = sbwt.get_subset_rank_structure().A_bits;
    const sdsl::bit_vector& C_bits = sbwt.get_subset_rank_structure().C_bits;
    const sdsl::bit_vector& G_bits = sbwt.get_subset_rank_structure().G_bits;
    const sdsl::bit_vector& T_bits = sbwt.get_subset_rank_structure().T_bits;

    const sdsl::rank_support_v5<>& A_bits_rs = sbwt.get_subset_rank_structure().A_bits_rs;
    const sdsl::rank_support_v5<>& C_bits_rs = sbwt.get_subset_rank_structure().C_bits_rs;
    const sdsl::rank_support_v5<>& G_bits_rs = sbwt.get_subset_rank_structure().G_bits_rs;
    const sdsl::rank_support_v5<>& T_bits_rs = sbwt.get_subset_rank_structure().T_bits_rs;

//    cout << " Extracting k" << endl;
    const int64_t k = sbwt.get_k();

//    cout << " Extracting C array"<< endl;
    const vector<int64_t>& C = sbwt.get_C_array();

    vector<pair<int64_t,int64_t> > ans;
    const sdsl::bit_vector* DNA_bitvectors[4] = {&A_bits, &C_bits, &G_bits, &T_bits};
    const sdsl::rank_support_v5<>* DNA_rs[4] = {&A_bits_rs, &C_bits_rs, &G_bits_rs, &T_bits_rs};

    uint64_t  bit_len = A_bits.size(); // all the same

//    cout<< " Searching for the first char"<< endl;
    // Search the first char
    // Start from the middle of the string
    int64_t start = bit_len/2;
    int64_t curr_pos = start;
    int64_t match_len = 0;
    for(int64_t i = 0; i < len; i++){
//        for (auto res :ans){
//            cout<< res.first << "  and " << res.second << endl;
//        }
        char c = static_cast<char>(input[i] & ~32); // convert to uppercase using a bitwise operation //char c = toupper(input[i]);
        int64_t char_idx = get_char_idx(c);
        if(char_idx == -1) {// Invalid character
            match_len = 0;
            curr_pos = -1;
        } else {
            const sdsl::bit_vector &Bit_v = *(DNA_bitvectors[char_idx]);
            const sdsl::rank_support_v5<> &Bit_rs = *(DNA_rs[char_idx]);
            if (curr_pos == -1) {
                curr_pos = start; // restart from the middle of the string
            }
            // CASE 1
            if (Bit_v[curr_pos] == 1) { // I.second
                //match_len = ((match_len) < k ? match_len + 1 : k); //  no need to chek since the max length of LCS if k-1, no kmers are the same
                match_len += (match_len < k);
                //match_len += -(match_len < k) & 1;
                curr_pos = C[char_idx] + Bit_rs(curr_pos);
            } else {
                // CASE 2
                int64_t r = curr_pos - 1, l = curr_pos + 1;
                while (r >= 0 && !Bit_v[r]) {
                    r--;
                }
                while (l < bit_len && !Bit_v[l]) {
                    l++;
                }
                if (r >= 0){
                    if (l < bit_len) { // both ok
                        uint64_t r_lcs = LCS[rmqLCS(r + 1, curr_pos)];
                        uint64_t l_lcs = LCS[rmqLCS(curr_pos + 1, l)];
                        int64_t next_lcs = r_lcs;
                        next_lcs += -(r_lcs > l_lcs) & (l_lcs - r_lcs);

                        //pair<int64_t, int64_t> next = max(make_pair(r_lcs, r), make_pair(l_lcs, l));
                        //match_len = ((next.first <= match_len) ? next.first : match_len) + 1;

                        //if (a > b) { x += y;}
                        // Branchless
                        //x += -(a > b) & y;
                        match_len -= -(match_len > next_lcs) & (match_len-next_lcs); // add ++

                        //x = (a > b) ? z : y;
                        // Branchless
                        //x = y ^ ((z ^ y) & -(a > b));
                        int64_t next_pos = l ^ ((r ^ l) & -(r_lcs > l_lcs));
                        curr_pos = C[char_idx] + Bit_rs(next_pos);
                    } else { // r ok
                        uint64_t r_lcs = LCS[rmqLCS(r + 1, curr_pos)];
                        //match_len = (r_lcs <= match_len) ? r_lcs + 1 : match_len + 1;
                        match_len -= -(match_len > r_lcs) & (match_len-r_lcs); // add ++
                        curr_pos = C[char_idx] + Bit_rs(r);
                    }
                } else if (l < bit_len) { // l ok
                    uint64_t l_lcs = LCS[rmqLCS(curr_pos + 1, l)];
                    //match_len = (l_lcs <= match_len) ? l_lcs + 1 : match_len + 1;
                    match_len -= -(match_len > l_lcs) & (match_len - l_lcs); //1234 add ++
                    curr_pos = C[char_idx] + Bit_rs(l);
                } else { // if the char is a DNA char but is not in the input SBWT
                    match_len = -1;
                    curr_pos = -1;
                }
                match_len++;
            }
        }
        ans.push_back({match_len,curr_pos});
    }
    return ans;
}
//vector<pair<int64_t,int64_t> > new_streaming_search(const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, const string& query) {
//    return new_streaming_search(sbwt, LCS, query.c_str(), query.size());
//}
