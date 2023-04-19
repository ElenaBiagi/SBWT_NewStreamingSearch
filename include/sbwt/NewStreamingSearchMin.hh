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
vector<pair<int64_t,int64_t> > new_streaming_search_min(const sdsl::bit_vector** DNA_bitvectors, const sdsl::rank_support_v5<>** DNA_rs, const vector<int64_t>& C, const int64_t k, const sdsl::int_vector<>& LCS, const char* input, int64_t len){ //const plain_matrix_sbwt_t& sbwt

    vector<pair<int64_t,int64_t>> ans;

    // Search the first char
    // Start from the middle of the string
    const int64_t bit_len = DNA_bitvectors[0]->size();
    int64_t start = bit_len/2;
    int64_t curr_pos = start;
    uint64_t match_len = 0;
    for(int64_t i = 0; i < len; i++){

        char c = static_cast<char>(input[i] & ~32); // convert to uppercase using a bitwise operation //char c = toupper(input[i]);
        int64_t char_idx = get_char_idx(c);
        if(char_idx == -1) [[unlikely]]{// Invalid character
            match_len = 0;
            curr_pos = -1;
        } else {
            const sdsl::bit_vector& Bit_v = *(DNA_bitvectors[char_idx]);
            const sdsl::rank_support_v5<>& Bit_rs = *(DNA_rs[char_idx]);
            if (curr_pos == -1) [[unlikely]]{
                curr_pos = start; // restart from the middle of the string
            }
            // CASE 1
            if (Bit_v[curr_pos]) {
                match_len += (match_len < k);
                curr_pos = C[char_idx] + Bit_rs(curr_pos+1)-1;
            } else {
                // CASE 2
                int64_t r = curr_pos -1, l = curr_pos + 1 ;
                uint64_t r_lcs = LCS[curr_pos],  l_lcs = LCS[curr_pos + 1]; // maximum possible match
                while(r >= 0 && !Bit_v[r]){
                    r_lcs = min(r_lcs, LCS[r]);
                    //if (LCS[r]< r_lcs ){r_lcs = LCS[r];}
                    r--;
                }
                while(l < bit_len-1 && !Bit_v[l]){
                    l_lcs = min(l_lcs, LCS[l+1]);
                    //if (LCS[l+1]< l_lcs){l_lcs = LCS[l+1];}
                    l++;
                }
                if (r >= 0) [[likely]]{
                    if (l < bit_len) [[likely]]{ // both ok
                        pair<uint64_t, int64_t> next = max(make_pair(r_lcs, r), make_pair(l_lcs, l));
                        match_len = std::min(match_len, next.first) + 1; //match_len = ((next.first < match_len) ? next.first : match_len) + 1; //match_len -= -(match_len > next_lcs) & (match_len-next_lcs); // add ++
                        int64_t next_pos = next.second; //int64_t next_pos = l ^ ((r ^ l) & -(r_lcs > l_lcs));
                        curr_pos = C[char_idx] + Bit_rs(next_pos);
                    } else { // r ok
                        match_len = min(match_len, r_lcs) + 1;
                        curr_pos = C[char_idx] + Bit_rs(r);
                    }
                } else if (l < bit_len) { // l ok
                    match_len = min(match_len, l_lcs) + 1;
                    curr_pos = C[char_idx] + Bit_rs(l);
                } else [[unlikely]]{ // if the char is a DNA char but is not in the input SBWT
                    match_len = 0;
                    curr_pos = -1;
                }
            }
        }
        ans.push_back({match_len,curr_pos});
    }
    return ans;
}

