#pragma once

#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include "SBWT.hh"
#include "variants.hh"
#include "sdsl/vectors.hpp"
#include "sdsl/int_vector.hpp"

#include "suffix_group_optimization.hh"

using namespace std;
using namespace sbwt;

// defined in NewStreamingSearchRMQ.hh
//int64_t get_char_idx(char c){
//    switch(c){
//        case 'A': return 0;
//        case 'C': return 1;
//        case 'G': return 2;
//        case 'T': return 3;
//        default: return -1;
//    }
//}

// Returns pairs (a_1, b_1), (a_2, b_2)..., such that 
// - a_i is the length of the longest match ending at query[i]
// - b_i is the colex rank of one arbitrary k-mer having the longest match.
vector<pair<int64_t,int64_t> > new_streaming_search_min(const sdsl::bit_vector** DNA_bitvectors, const sdsl::rank_support_v5<>** DNA_rs, const int64_t bit_len, const vector<int64_t>& C, const int64_t k, const sdsl::int_vector<>& LCS, const char* input, int64_t len){ //

    vector<pair<int64_t,int64_t> > ans;
    // Searching the first char start from the middle of the string
    int64_t start = bit_len/2;
    int64_t curr_pos = start;
    int64_t match_len = 0;
    for(int64_t i = 0; i < len; i++){
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
                match_len += (match_len < k);
                curr_pos = C[char_idx] + Bit_rs(curr_pos);
            } else {
                // CASE 2
                int64_t r = curr_pos - 1, l = curr_pos + 1;
                int64_t r_lcs = k,  l_lcs = k; // maximum possible match
                while(r >= 0 && !Bit_v[r]){
                    if (LCS[r-1]< r_lcs){r_lcs = LCS[r-1];}
                    r--;
                }
                while(l < bit_len && !Bit_v[l]){
                    if (LCS[l]< l_lcs){l_lcs = LCS[l];}
                    l++;
                }
                if (r >= 0) [[likely]]{
                    if (l < bit_len) [[likely]]{ // both ok
                        int64_t next_lcs = r_lcs;
                        next_lcs += -(r_lcs > l_lcs) & (l_lcs - r_lcs);
                        match_len -= -(match_len > next_lcs) & (match_len-next_lcs); // add ++
                        int64_t next_pos = l ^ ((r ^ l) & -(r_lcs > l_lcs));
                        curr_pos = C[char_idx] + Bit_rs(next_pos);
                    } else { // r ok
                        match_len -= -(match_len > r_lcs) & (match_len-r_lcs); // add ++
                        curr_pos = C[char_idx] + Bit_rs(r);
                    }
                } else if (l < bit_len) { // l ok
                    match_len -= -(match_len > l_lcs) & (match_len - l_lcs); // add ++
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
